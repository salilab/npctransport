/**
 *  \file Scoring.cpp
 *  \brief description.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/Scoring.h>
#include <IMP/npctransport/SimulationData.h>
#include <IMP/npctransport/SitesPairScore.h>
#include <IMP/npctransport/SlabSingletonScore.h>
#ifdef IMP_NPC_GOOGLE
IMP_GCC_PUSH_POP(diagnostic push)
IMP_GCC_PRAGMA(diagnostic ignored "-Wsign-compare")
#include "third_party/npc/npctransport/data/npctransport.pb.h"
IMP_GCC_PUSH_POP(diagnostic pop)
#else
#include <IMP/npctransport/internal/npctransport.pb.h>
#endif
#include <IMP/npctransport/typedefs.h>
#include <IMP/npctransport/util.h>
#include <IMP/base_types.h>
#include <IMP/container_macros.h>
#include <IMP/algebra/vector_generators.h>
#include <IMP/atom/estimates.h>
#include <IMP/atom/distance.h>
#include <IMP/atom/Diffusion.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/Selection.h>
#include <IMP/base/exception.h>
#include <IMP/base/log.h>
#include <IMP/container/CloseBipartitePairContainer.h>
#include <IMP/container/ConsecutivePairContainer.h>
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/container/generic.h>
#include <IMP/core/BoundingBox3DSingletonScore.h>
#include <IMP/core/DistancePairScore.h>
#include <IMP/core/SphereDistancePairScore.h>
#include <IMP/core/HarmonicUpperBound.h>
#include <IMP/core/pair_predicates.h>
#include <IMP/core/RestraintsScoringFunction.h>
#include <IMP/core/XYZ.h>

#include <numeric>
#include <set>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
#define GET_ASSIGNMENT(name) name##_ = data.name().value()
#define GET_VALUE(name) name##_ = data.name()

Scoring::Scoring
(SimulationData* owner_sd,
 const ::npctransport_proto::Assignment &data)
: Object("Scoring%1%"),
  owner_sd_(owner_sd),
  //  is_updating_particles_(false),
  close_diffusers_container_(nullptr),
  otpp_(new core::OrderedTypePairPredicate()),
  scoring_function_(nullptr),
  predr_(nullptr),
  box_restraint_(nullptr),
  slab_restraint_(nullptr)
{
  IMP_ALWAYS_CHECK(owner_sd_ != nullptr,
                   "Must have non-null owner simulation data",
                   IMP::base::ValueException);
  GET_ASSIGNMENT(box_side);
  GET_ASSIGNMENT(tunnel_radius);
  GET_ASSIGNMENT(slab_thickness);
  GET_ASSIGNMENT(slab_is_on);
  GET_ASSIGNMENT(box_is_on);
  GET_ASSIGNMENT(interaction_k);
  GET_ASSIGNMENT(interaction_range);
  GET_ASSIGNMENT(backbone_k);
  GET_ASSIGNMENT(slack);
  GET_ASSIGNMENT(nonspecific_k);
  GET_ASSIGNMENT(nonspecific_range);
  GET_ASSIGNMENT(excluded_volume_k);
  GET_VALUE(range);
  //  update_particles();
}

IMP::ScoringFunction*
Scoring::get_scoring_function(bool update)
{
  if(!scoring_function_ || update){
    // set up the restraints for the BD simulation:
    RestraintsTemp rs =
      get_chain_restraints_on( get_sd()->get_diffusers() );
    if (box_is_on_) {
      rs.push_back(get_bounding_box_restraint(update));
    }
    if (slab_is_on_) {
      rs.push_back(get_slab_restraint(update));
    }
    rs.push_back(this->get_predicates_pair_restraint(update));
    //  is_updating_particles_ = false; // everything is supposed to be updated now
    //  IMP_NEW(core::RestraintsScoringFunction, rsf, (rs));
    scoring_function_  = new core::RestraintsScoringFunction(rs);
  }
  return scoring_function_;
}

IMP::ScoringFunction*
Scoring::get_custom_scoring_function
( const RestraintsTemp& extra_restraints,
  IMP::SingletonContainerAdaptor particles,
  IMP::SingletonContainerAdaptor optimizable_particles,
  bool is_attr_interactions_on) const
{
    // set up the restraints for the BD simulation:
  RestraintsTemp rs;
  rs += extra_restraints;
  rs += get_chain_restraints_on( particles );
  if (box_is_on_) {
    rs.push_back( create_bounding_box_restraint( particles ) );
  }
  if (slab_is_on_) {
    rs.push_back( create_slab_restraint( particles ) );
  }
  PairContainer* cpc = create_close_diffusers_container
    ( particles, optimizable_particles );
  rs.push_back( create_predicates_pair_restraint
                ( cpc, is_attr_interactions_on ) );
  //  is_updating_particles_ = false; // everything is supposed to be updated now
  IMP_NEW(core::RestraintsScoringFunction, rsf, (rs));
  return rsf.release();
}


IMP::PairContainer *
Scoring::get_close_diffusers_container(bool update)
{
  if(!close_diffusers_container_ || update){
    close_diffusers_container_ =
      create_close_diffusers_container
      ( get_sd()->get_diffusers(),
        get_sd()->get_optimizable_diffusers() );
      }
  return close_diffusers_container_;
}

IMP::container::PredicatePairsRestraint *
Scoring::get_predicates_pair_restraint
( bool update )
{
  if(!predr_ || update){
    predr_ = create_predicates_pair_restraint
      ( get_close_diffusers_container(), true /* attr restraints */ );
  }
  return predr_;
}

Restraint *
Scoring::get_bounding_box_restraint(bool update)
{
  IMP_USAGE_CHECK(box_is_on_, "box is not on - can't get restraint");
  if (update || !box_restraint_) {
    box_restraint_ =
      create_bounding_box_restraint( get_sd()->get_diffusers() );
  }
  return box_restraint_;
}

Restraint *
Scoring::get_slab_restraint(bool update)
{
  IMP_USAGE_CHECK(slab_is_on_, "slab is not on - can't get restraint");
  if (update || !slab_restraint_) {
    slab_restraint_ =
      create_slab_restraint( get_sd()->get_diffusers() );
  }
  return slab_restraint_;
}


/**
   add a pair score restraint that applies to particles of
   types t0 and t1 to the PredicatePairsRestraint object returned by
   get_predr().

   A SitesPairScore interaction means site-specific
   attractive forces between bidning sites on each particle,
   and non-specific attraction and repulsion (upon penetraion)
   between the particles themselves.

   \see SitesPairScore
   \see create_sites_pair_score
*/
void Scoring::add_interaction
( const ::npctransport_proto::Assignment_InteractionAssignment &idata)
{
  // extract interaction params
  core::ParticleType type0(idata.type0());
  core::ParticleType type1(idata.type1());
  double base_k = interaction_k_;
  if (idata.has_interaction_k()) {
    base_k = idata.interaction_k().value();
  }
  // no particles so drop it
  if (interaction_k_factors_.find(type0) == interaction_k_factors_.end() ||
      interaction_k_factors_.find(type1) == interaction_k_factors_.end()) {
    return;
  }
  double interaction_k = base_k * interaction_k_factors_.find(type0)
                                      ->second  // TODO: validate type exists
                         *
                         interaction_k_factors_.find(type1)->second;
  double base_range = interaction_range_;
  if (idata.has_interaction_range()) {
    base_range = idata.interaction_range().value();
  }
  double interaction_range =
      base_range * interaction_range_factors_.find(type0)->second *
      interaction_range_factors_.find(type1)->second;

  IMP_LOG(IMP::base::PROGRESS,
          "creating interaction "
            << idata.type0() << "," << idata.type1()
            << " effective_k = " << interaction_k
            << ", effective_range = " << interaction_range
            << ", nonspecific k = " << nonspecific_k_
            << ", nonspecific range = " << nonspecific_range_
            << ", excluded volume k = " << excluded_volume_k_
            << std::endl
          );

  // add the interaction restraint both for (t0,t1) and (t1,t0)  {
    // TODO: repulsion was also added in get_predr - do we double count here?
  {
    core::ParticleTypes pts1;
    pts1.push_back(type0);
    pts1.push_back(type1);
    int interaction_id1 = get_ordered_type_pair_predicate()->get_value(pts1);
    IMP::PairScore* ps1 = create_sites_pair_score
      ( interaction_range,                  // site-specific
        interaction_k, nonspecific_range_,  // non-specific = entire particle
        nonspecific_k_, excluded_volume_k_,
        get_sd()->get_sites(type0),
        get_sd()->get_sites(type1)
        );
   interaction_pair_scores_[interaction_id1] = ps1;
  }
  {
    core::ParticleTypes pts2;
    pts2.push_back(type1);
    pts2.push_back(type0);
    int interaction_id2 = get_ordered_type_pair_predicate()->get_value(pts2);
    IMP::PairScore* ps2 = create_sites_pair_score
      (  interaction_range,                  // site-specific
         interaction_k, nonspecific_range_,  // non-specific = entire particle
         nonspecific_k_, excluded_volume_k_,
         get_sd()->get_sites(type0),
         get_sd()->get_sites(type1)
         );
    interaction_pair_scores_[interaction_id2] = ps2;
  }
}



// returns max of site-site and non-specific interaction ranges
// for interaction type it, if applicable
double Scoring::get_interaction_range_for
  ( core::ParticleType t1, core::ParticleType t2,
     bool site_specific, bool non_specific) const
{
  // Get pair score
  core::ParticleTypes pts;
  pts.push_back(t1);  pts.push_back(t2);
  int iid = const_cast<core::OrderedTypePairPredicate*>
    (get_ordered_type_pair_predicate())->get_value(pts);
  t_map_pair_type_to_pair_score::const_iterator iter =
    interaction_pair_scores_.find(iid);
  IMP_USAGE_CHECK(iter != interaction_pair_scores_.end(),
                  "looking for a range for a non-existent interaction type");
  IMP::PairScore const* ps = iter->second;
  // Get ranges (assuming Sites or LinearInteraction pair scores)
  double range_ss = -1.0; // site-specific
  if(site_specific){
    SitesPairScore const* sps =
      dynamic_cast< SitesPairScore const* >( ps );
    if(sps) range_ss = sps->get_sites_range();
  }
  double range_ns = -1.0; // non-specific
  if(non_specific){
    LinearInteractionPairScore const* lips =
      dynamic_cast< LinearInteractionPairScore const* >( ps );
    if(lips) range_ns = lips->get_range_attraction();
  }
  return std::max(range_ss, range_ns);
}



/** add a restraint on particles in an FG nup chain */
Restraint*
Scoring::add_chain_restraint(atom::Hierarchy chain_root,
                             double rest_length,
                             std::string name)
{
  IMP_ALWAYS_CHECK( chain_root.get_number_of_children() > 0,
                    "No Particles passed.", IMP::base::ValueException );
  IMP_ALWAYS_CHECK( rest_length > 0, "negative rest length invalid",
                   IMP::base::ValueException );

  // Exclusive means that the particles will be in no other
  // ConsecutivePairContainer this assumption accelerates certain computations
  IMP_NEW(IMP::container::ExclusiveConsecutivePairContainer, xcpc,
          (chain_root.get_children(), name + " consecutive pairs"));
  // add chain restraint
  IMP_NEW(LinearWellPairScore, lwps,
          ( rest_length, get_backbone_k() ) );
  IMP::base::Pointer<Restraint> cr =
    IMP::container::create_restraint
    ( lwps.get(), xcpc.get(), "chain restraint %1%" );
  chain_scores_.push_back( lwps );
  ParticleIndex pi_chain_root = chain_root.get_particle_index();
  chain_restraints_map_[pi_chain_root] = cr;
  return cr;
}


Restraints
Scoring::get_chain_restraints_on
( IMP::SingletonContainerAdaptor particles ) const
{
  Restraints ret_rs;
  IMP_CONTAINER_FOREACH
    ( container::ListSingletonContainer,
      particles,
      {
        if(get_sd()->get_is_fg(_1)){
          IMP_USAGE_CHECK(atom::Hierarchy::get_is_setup(get_model(), _1),
                          "all particles in list must be hierarchy compatible");
          atom::Hierarchy h=atom::Hierarchy(get_model(), _1);
          ParticleIndex pi_dad = h.get_parent().get_particle_index();
          if( chain_restraints_map_.find(pi_dad) !=
              chain_restraints_map_.end() )
            {
              ret_rs.push_back
                ( chain_restraints_map_.find(pi_dad)->second );
            }
        }
      }
      );
  return ret_rs;
}


/***************************************************************/
/***************************** Creators ************************/
/***************************************************************/

// a close pair container for all pair of particles and
// optimizable particles
IMP::PairContainer*
Scoring::create_close_diffusers_container
( SingletonContainerAdaptor particles,
  SingletonContainerAdaptor optimizable_particles) const
{
  using namespace container;
  IMP_NEW(CloseBipartitePairContainer, cpc, // so range + 2*slack is what we get
          (particles, optimizable_particles, get_range(), slack_) );
  IMP_NEW( core::AllSamePairPredicate, aspp, () );
  cpc->add_pair_filter( aspp ); // only relevant for bipartite
  // IMP_NEW( container::ClosePairContainer, cpc,
  //            particles,
  //            get_range(), slack_) );
  IMP_NEW( ExclusiveConsecutivePairFilter, ecpf, () );
  cpc->add_pair_filter( ecpf );
  return cpc.release();
}


container::PredicatePairsRestraint
*Scoring::create_predicates_pair_restraint
(PairContainer* pair_container,
 bool is_attr_interactions_on) const
{
  // set linear repulsion upon penetration between all close pairs
  // returned by get_close_diffusers_container(), with different
  // scores for interactions between particles of different
  // (ordered) types
  IMP_NEW(container::PredicatePairsRestraint, predr,
          ( otpp_,
            pair_container ) );
  IMP_NEW(LinearSoftSpherePairScore, ssps, (get_excluded_volume_k()));
  predr->set_unknown_score(ssps.get());
  // add interaction scores with site-specific interactions
  if(is_attr_interactions_on) {
    t_map_pair_type_to_pair_score::const_iterator iter;
    for(iter = interaction_pair_scores_.begin();
        iter != interaction_pair_scores_.end();
        iter++)
      {
        // first = interaction_id; second = pair-score
        // TODO: using SitesPairScore as map type would be more efficient
        //       as set_score() method is templated?
        predr->set_score( iter->first, iter->second );
      }
  }
  return predr.release();
}



/**
   Creates bounding volume restraints such as box restraint and slab restraints,
   based on the box_size_, slab_height_, slab_radius_, etc. class variables
*/
Restraint* Scoring::create_bounding_box_restraint
( SingletonContainerAdaptor particles) const
{
  // Add bounding box restraint
  // TODO: what does backbone_spring_k_ has to do
  //       with bounding box constraint?
  IMP_NEW(core::HarmonicUpperBound, hub, (0, excluded_volume_k_));
  IMP_NEW(core::GenericBoundingBox3DSingletonScore<core::HarmonicUpperBound>,
          bbss, (hub.get(), get_sd()->get_box()));
  return
    container::create_restraint(bbss.get(),
                                particles.get(),
                                "bounding box");
}

Restraint * Scoring::create_slab_restraint
( SingletonContainerAdaptor particles)  const
{
  // Add cylinder restraint
  IMP_NEW(SlabSingletonScore, slab_score,
          (slab_thickness_ /* h */, tunnel_radius_ /*r*/, excluded_volume_k_));
  return
    container::create_restraint(slab_score.get(),
                                particles.get(),
                                "bounding slab");
}


/************************************************************/
/************* various simple getters and setters *******************/
/************************************************************/



/** returns the model associated with the owned SimulationData */
Model* Scoring::get_model() {
  return get_sd()->get_model();
}


/** returns the model associated with the owned SimulationData */
Model * Scoring::get_model() const {
  return get_sd()->get_model();
}


#undef GET_ASSIGNMENT
#undef GET_VALUE

IMPNPCTRANSPORT_END_NAMESPACE
