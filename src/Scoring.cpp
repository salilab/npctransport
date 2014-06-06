/**
 *  \file Scoring.cpp
 *  \brief description.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/Scoring.h>

#include <IMP/npctransport/FGChain.h>
#include <IMP/npctransport/SimulationData.h>
#include <IMP/npctransport/SitesPairScore.h>
#include <IMP/npctransport/SlabSingletonScore.h>
#include <IMP/npctransport/ZBiasSingletonScore.h>
#include <IMP/npctransport/internal/npctransport.pb.h>
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
#include <IMP/exception.h>
#include <IMP/log.h>
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
#include <limits>
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
  close_beads_container_(nullptr),
  otpp_(new core::OrderedTypePairPredicate()),
  scoring_function_(nullptr),
  predr_(nullptr),
  box_restraint_(nullptr),
  slab_restraint_(nullptr)
{
  IMP_ALWAYS_CHECK(owner_sd_ != nullptr,
                   "Must have non-null owner simulation data",
                   IMP::ValueException);
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
    ParticlesTemp beads = get_sd()->get_beads();
    RestraintsTemp rs =
      get_chain_restraints_on( beads );
    if (box_is_on_) {
      rs.push_back(get_bounding_box_restraint(update));
    }
    if (slab_is_on_) {
      rs.push_back(get_slab_restraint(update));
    }
    rs += get_z_bias_restraints();
    rs += get_custom_restraints();
    rs.push_back(this->get_predicates_pair_restraint(update));

    scoring_function_  = new core::RestraintsScoringFunction(rs);
  }
  return scoring_function_;
}

IMP::ScoringFunction*
Scoring::get_custom_scoring_function
( const RestraintsTemp& extra_restraints,
  IMP::SingletonContainerAdaptor beads,
  IMP::SingletonContainerAdaptor optimizable_beads,
  bool is_attr_interactions_on) const
{
  beads.set_name_if_default("NPCGetCustomScoringFunctionBeadsInput%1%");
  beads.set_name_if_default("NPCGetCustomScoringFunctionOptimiedBeadsInput%1%");
    // set up the restraints for the BD simulation:
  RestraintsTemp rs;
  rs += extra_restraints;
  rs += get_chain_restraints_on( beads );
  if (box_is_on_) {
    rs.push_back( create_bounding_box_restraint( beads ) );
  }
  if (slab_is_on_) {
    rs.push_back( create_slab_restraint( beads ) );
  }
  rs += get_custom_restraints(); // TODO: this is problematic cause not restricted to beads - need to decide
  PairContainer* cpc = create_close_beads_container
    ( beads, optimizable_beads );
  rs.push_back( create_predicates_pair_restraint
                ( cpc, is_attr_interactions_on ) );
  IMP_NEW(core::RestraintsScoringFunction, rsf, (rs));
  return rsf.release();
}


IMP::PairContainer *
Scoring::get_close_beads_container(bool update)
{
  if(!close_beads_container_ || update){
    close_beads_container_ =
      create_close_beads_container
      ( get_sd()->get_beads(),
        get_sd()->get_optimizable_beads() );
      }
  return close_beads_container_;
}

IMP::container::PredicatePairsRestraint *
Scoring::get_predicates_pair_restraint
( bool update )
{
  if(!predr_ || update){
    predr_ = create_predicates_pair_restraint
      ( get_close_beads_container(), true /* attr restraints */ );
  }
  return predr_;
}

Restraint *
Scoring::get_bounding_box_restraint(bool update)
{
  IMP_USAGE_CHECK(box_is_on_, "box is not on - can't get restraint");
  if (update || !box_restraint_) {
    box_restraint_ =
      create_bounding_box_restraint( get_sd()->get_beads() );
  }
  return box_restraint_;
}

Restraint *
Scoring::get_slab_restraint(bool update)
{
  IMP_USAGE_CHECK(slab_is_on_, "slab is not on - can't get restraint");
  if (update || !slab_restraint_) {
    slab_restraint_ =
      create_slab_restraint( get_sd()->get_beads() );
  }
  return slab_restraint_;
}


/**
   add a pair score restraint that applies to beads of
   types t0 and t1 to the PredicatePairsRestraint object returned by
   get_predr().

   A SitesPairScore interaction means site-specific
   attractive forces between bidning sites on each particle,
   and non-specific attraction and repulsion (upon penetraion)
   between the beads themselves.

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
  // no beads so drop it
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

  IMP_LOG(IMP::PROGRESS,
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



// add FG Nup chain with restraints
void
Scoring::add_chain_restraints(FGChain* chain)
{
  IMP_ALWAYS_CHECK( chain->get_rest_length_factor() > 0,
                    "non-positive rest length factor is not valid",
                    IMP::base::ValueException );
  if(chain->get_backbone_k() <= 0 ){
    chain->set_backbone_k( get_default_backbone_k() );
  }
  // save chain into map so its restraints can be accessed later
  chains_set_.insert( chain );
  // save mapping from all beads back to chains for easy access
  for(unsigned int i = 0; i < chain->get_number_of_beads(); i++){
    ParticleIndex bi = chain->get_bead_index(i);
    bead_to_chain_map_[bi] = chain;
  }
}


IMP::Restraints
Scoring::get_chain_restraints_on
( IMP::SingletonContainerAdaptor bead_particles ) const
{
  bead_particles.set_name_if_default("GetChainRestraintsOnInput%1%");
  IMP::Restraints R;
  FGChainsSet chains_found; // to prevent duplicates
  IMP_CONTAINER_FOREACH
    ( SingletonContainer,
      bead_particles,
      {
        if(get_sd()->get_is_fg(_1) &&
           bead_to_chain_map_.find(_1) != bead_to_chain_map_.end())
          {
            FGChain* chain =
              bead_to_chain_map_.find(_1)->second.get();
            IMP_USAGE_CHECK
              (chains_set_.find(chain) != chains_set_.end(),
               "For some reason a bead particle appears in bead to chain map"
               << " but its chain is missing from chains_set_");
            if( chains_found.find(chain) == chains_found.end() )
            {
              chains_found.insert(chain);
              R += chain->get_chain_restraints();
            }
          }
      }
      ); // IMP_CONTAINER_FOREACH
  return R;
}

Restraints
Scoring::get_all_chain_restraints() const
{
  Restraints R;
  for ( FGChainsSet::iterator it = chains_set_.begin();
        it != chains_set_.end(); it++)
    {
      FGChain* chain = it->get();
      R += chain->get_chain_restraints();
    }
  return R;
}


/***************************************************************/
/***************************** Creators ************************/
/***************************************************************/

// a close pair container for all pair of specified beads and
// optimizable beads
IMP::PairContainer*
Scoring::create_close_beads_container
( SingletonContainerAdaptor beads,
  SingletonContainerAdaptor optimizable_beads) const
{
  beads.set_name_if_default(
      "CreateCloseBeadsContainerBeadsInput%1%");
  optimizable_beads.set_name_if_default(
      "CreateCloseBeadsContainerOptimizableBeadsInput%1%");
  using namespace container;
  IMP_NEW(CloseBipartitePairContainer, cpc, // so range + 2*slack is what we get
          (beads, optimizable_beads, get_range(), slack_) );
  IMP_NEW( core::AllSamePairPredicate, aspp, () );
  cpc->add_pair_filter( aspp ); // only relevant for bipartite version
  // The following is for chains (TODO: strongly tied to
  // implementation but then again, it is possibly important to
  // exclude bonded interactions - see FGChain.cpp. A possible
  // solution would be to add a static interface method in FGChain to
  // return a filter for bonded interactions):
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
  // returned by get_close_beads_container(), with different
  // scores for interactions between beads of different
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
( SingletonContainerAdaptor beads) const
{
  beads.set_name_if_default("CreateBoundingBoxRestraintInput%1%");
  // Add bounding box restraint
  // TODO: what does backbone_spring_k_ has to do
  //       with bounding box constraint?
  IMP_NEW(core::HarmonicUpperBound, hub, (0, excluded_volume_k_));
  IMP_NEW(core::GenericBoundingBox3DSingletonScore<core::HarmonicUpperBound>,
          bbss, (hub.get(), get_sd()->get_box()));
  return
    container::create_restraint(bbss.get(),
                                beads.get(),
                                "bounding box");
}

Restraint * Scoring::create_slab_restraint
( SingletonContainerAdaptor particles)  const
{
  particles.set_name_if_default("CreateStabRestraintInput%1%");
  // Add cylinder restraint
  IMP_NEW(SlabSingletonScore, slab_score,
          (slab_thickness_ /* h */, tunnel_radius_ /*r*/, excluded_volume_k_));
  return
    container::create_restraint(slab_score.get(),
                                particles.get(),
                                "bounding slab");
}

void Scoring::add_z_bias_restraint
( SingletonContainerAdaptor ps, double k )
{
  ps.set_name_if_default("AddZBiasRestraintInput%1%");
  z_bias_restraints_.push_back
    ( create_z_bias_restraint( ps, k) );
}

void Scoring::add_z_bias_restraint(Particle* p, double k)
{
  ParticlesTemp ps(1, p);
  add_z_bias_restraint(ps, k);
}

IMP::Restraints
Scoring::get_z_bias_restraints()
{
  return z_bias_restraints_;
}

IMP::Restraint*
Scoring::create_z_bias_restraint(SingletonContainerAdaptor ps, double k)
const
{
  ps.set_name_if_default("CreateZBiasRestraintKInput%1%");
  // restrict restraint to be over tunnel if applicable
  double R = std::numeric_limits<double>::max();
  if (get_sd()->get_has_slab()) {
    R = get_sd()->get_tunnel_radius();
  }
  IMP_NEW(ZBiasSingletonScore, zbsc, (k, R) );
  return
    container::create_restraint(zbsc.get(),
                                ps.get(),
                                "z-bias potential");
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
