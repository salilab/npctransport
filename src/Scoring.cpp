/** *  \file Scoring.cpp
 *  \brief description.
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/Scoring.h>

#include <IMP/npctransport/FGChain.h>
#include <IMP/npctransport/linear_distance_pair_scores.h>
#include <IMP/npctransport/HarmonicSpringSingletonScore.h>
#include <IMP/npctransport/SimulationData.h>
#include <IMP/npctransport/SitesPairScore.h>
#include <IMP/npctransport/SlabWithCylindricalPorePairScore.h>
#include <IMP/npctransport/SlabWithToroidalPorePairScore.h>
#include <IMP/npctransport/AnchorToCylindricalPorePairScore.h>
#include <IMP/npctransport/PoreRadiusSingletonScore.h>
#include <IMP/npctransport/ZBiasSingletonScore.h>
#include <IMP/npctransport/internal/npctransport.pb.h>
#include <IMP/npctransport/typedefs.h>
#include <IMP/npctransport/util.h>

#include <IMP/base_types.h>
#include <IMP/container_macros.h>
#include <IMP/container/AllBipartitePairContainer.h>
#include <IMP/algebra/vector_generators.h>
#include <IMP/atom/estimates.h>
#include <IMP/atom/distance.h>
#include <IMP/atom/Diffusion.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/Selection.h>
#include <IMP/exception.h>
#include <IMP/log.h>
#include <IMP/container/CloseBipartitePairContainer.h>
#include <IMP/container/ClosePairContainer.h>
#include <IMP/container/PairContainerSet.h>
#include <IMP/container/ConsecutivePairContainer.h>
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/container/SingletonContainerSet.h>
#include <IMP/container/generic.h>
#include <IMP/core/BoundingBox3DSingletonScore.h>
#include <IMP/core/BoundingSphere3DSingletonScore.h>
#include <IMP/core/DistancePairScore.h>
#include <IMP/core/SphereDistancePairScore.h>
#include <IMP/core/HarmonicUpperBound.h>
#include <IMP/core/pair_predicates.h>
#include <IMP/core/RestraintsScoringFunction.h>
#include <IMP/core/XYZ.h>

#include <numeric>
#include <limits>
#include <set>
#include <string>
#include <boost/tuple/tuple.hpp>
#if defined(_MSC_VER)
#include <io.h>
#else
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

IMPNPCTRANSPORT_BEGIN_NAMESPACE
#define GET_ASSIGNMENT(name) name##_ = pb_assignment.name().value()
#define GET_VALUE(name) name##_ = pb_assignment.name()

Scoring::Scoring
(SimulationData* owner_sd,
 const ::npctransport_proto::Assignment &pb_assignment)
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
  GET_ASSIGNMENT(box_is_on);
  GET_ASSIGNMENT(interaction_k);
  GET_ASSIGNMENT(interaction_range);
  GET_ASSIGNMENT(backbone_k);
  GET_VALUE(is_backbone_harmonic);
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
    if (get_has_bounding_volume()) {
      rs.push_back(get_bounding_volume_restraint(update));
    }
    if (get_sd()->get_has_slab()) {
      rs.push_back(get_slab_restraint(update));
      if(get_sd()->get_is_pore_radius_dynamic()){
        rs += anchor_restraints_;
        rs.push_back(this->get_pore_radius_restraint());
      }
    }
    rs += get_z_bias_restraints();
    rs += get_custom_restraints();
    rs.push_back(this->get_predicates_pair_restraint(update));

    scoring_function_  = new core::RestraintsScoringFunction(rs);
    scoring_function_rs_ = rs;
  }
  return scoring_function_;
}

IMP::Restraints
Scoring::get_scoring_function_restraints(bool force_update)
{
  if(!scoring_function_ || force_update){
    get_scoring_function(force_update);
  }
  return scoring_function_rs_; // updated in get_scoring_function()
}

IMP::ScoringFunction*
Scoring::get_custom_scoring_function
( const RestraintsTemp& extra_restraints,
  IMP::SingletonContainerAdaptor non_optimizable_beads,
  IMP::SingletonContainerAdaptor optimizable_beads,
  bool is_attr_interactions_on) const
{
  non_optimizable_beads.set_name_if_default("NPCGetCustomScoringFunctionNonOptimizableBeadsInput%1%");
  optimizable_beads.set_name_if_default("NPCGetCustomScoringFunctionOptimizableBeadsInput%1%");
  IMP::SingletonContainers beads_list;
  beads_list.push_back(non_optimizable_beads);
  beads_list.push_back(optimizable_beads);
  IMP_NEW(IMP::container::SingletonContainerSet, beads, (beads_list));

    // set up the restraints for the BD simulation:
  RestraintsTemp rs;
  rs += extra_restraints;
  rs += get_chain_restraints_on( beads );
  if (get_has_bounding_volume()) {
    rs.push_back( create_bounding_volume_restraint( beads ) );
  }
  if (get_sd()->get_has_slab()) {
    rs.push_back( create_slab_restraint( beads ) );
    if(get_sd()->get_is_pore_radius_dynamic()){
      rs += anchor_restraints_; //TODO: is issue that may include unoptimized anchor beads?)
      rs.push_back(this->get_pore_radius_restraint());
    }
  }
  rs += get_custom_restraints(); // TODO: this is problematic cause not restricted to beads - need to decide
  PairContainer* cpc = create_close_beads_container
    ( non_optimizable_beads, optimizable_beads );
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
      ( get_particle_indexes(get_sd()->get_non_optimizable_beads()),
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
      ( get_close_beads_container(update), true /* attr restraints */ );
  }
  return predr_;
}

Restraint *
Scoring::get_bounding_volume_restraint(bool update)
{
  IMP_USAGE_CHECK(get_has_bounding_volume(),
                  "box is not on - can't get restraint");
  if (update || !box_restraint_) {
    box_restraint_ =
      create_bounding_volume_restraint( get_sd()->get_beads() );
  }
  return box_restraint_;
}

Restraint *
Scoring::get_slab_restraint(bool update)
{
  IMP_USAGE_CHECK(get_sd()->get_has_slab(),
                  "slab is not on - can't get restraint");
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
   attractive forces between binding sites on each particle,
   and non-specific attraction and repulsion (upon penetraion)
   between the beads themselves.

   \see SitesPairScore
   \see create_sites_pair_score
*/
void Scoring::add_interaction
( const ::npctransport_proto::Assignment_InteractionAssignment &idata)
{
  if(idata.is_on().value()==0){
    IMP_LOG(VERBOSE,
            "Skipping disabled interaction between " << idata.type0()
              << " and " << idata.type1() << std::endl);
    return;
  }
  // extract interaction params
  core::ParticleType type0(idata.type0());
  core::ParticleType type1(idata.type1());
  // no beads so drop it (using interaction k just as a proxy to test existence)
  if (interaction_k_factors_.find(type0) == interaction_k_factors_.end() ||
      interaction_k_factors_.find(type1) == interaction_k_factors_.end()) {
    return;
  }

  // compute factored k, range, and sigmas
  double base_k = interaction_k_;
  if (idata.has_interaction_k()) {
    base_k = idata.interaction_k().value();
  }
  double interaction_k = base_k * interaction_k_factors_.find(type0)->second
                                * interaction_k_factors_.find(type1)->second;
  double base_range = interaction_range_;
  if (idata.has_interaction_range()) {
    base_range= idata.interaction_range().value();
  }
  double interaction_range =
      base_range * interaction_range_factors_.find(type0)->second *
                   interaction_range_factors_.find(type1)->second;
  double sigma0=0.0, sigma1=0.0;
  if (idata.has_range_sigma0_deg() &&
      idata.has_range_sigma0_deg()) {
    sigma0 = idata.range_sigma0_deg().value();
    sigma1 = idata.range_sigma1_deg().value();
  }
  double nonspecific_k= nonspecific_k_;
  if (idata.has_nonspecific_k()){
    nonspecific_k= idata.nonspecific_k().value();
  }
  double nonspecific_range= nonspecific_range_;
  if (idata.has_nonspecific_range()){
    nonspecific_range= idata.nonspecific_range().value();
  }
  double excluded_volume_k= excluded_volume_k_;
  if (idata.has_excluded_volume_k()){
    excluded_volume_k= idata.excluded_volume_k().value();
  }
  IMP_LOG(IMP::PROGRESS,
          "creating interaction "
	  << idata.type0() << "," << idata.type1()
	  << " effective_k = " << interaction_k
	  << ", effective_range = " << interaction_range
	  << ", sigma0/1 = " << sigma0 << "," << sigma1
	  << ", nonspecific k = " << nonspecific_k
	  << ", nonspecific range = " << nonspecific_range
	  << ", excluded volume k = " << excluded_volume_k
	  << std::endl
          );

  // get active sites (or all if no active sites specified)
  algebra::Sphere3Ds sites0;
  algebra::Sphere3Ds sites1;
  {
    const algebra::Sphere3Ds sites0_all( get_sd()->get_sites(type0) );
    const algebra::Sphere3Ds sites1_all( get_sd()->get_sites(type1) );
    int n0= idata.active_sites0_size();
    int n1= idata.active_sites1_size();
    if(n0>0){
      for(int i=0; i<n0; i++){
        unsigned int site_id= idata.active_sites0(i);
        IMP_ALWAYS_CHECK(site_id<sites0_all.size(),
                         "Invalid active site id " << site_id
                         << " for " << type0 << " interacting with " << type1
                         << std::endl, IMP::ValueException);
        sites0.push_back(sites0_all[site_id]);
      }
    }else{
      sites0= sites0_all;
    }
    if(n1>0){
      for(int i=0; i<n1; i++){
        unsigned int site_id= idata.active_sites1(i);
        IMP_ALWAYS_CHECK(site_id<sites1_all.size(),
                         "Invalid active site id " << site_id
                         << " for " << type1 << " interacting with " << type0
                         << std::endl, IMP::ValueException);
        sites1.push_back(sites1_all[site_id]);
      }
    }else{
      sites1= sites1_all;
    }
  }

  // add the interaction restraint both for (t0,t1) and (t1,t0)
  {
    core::ParticleTypes pts1;
    pts1.push_back(type0);
    pts1.push_back(type1);
    int interaction_id1 = get_ordered_type_pair_predicate()->get_value(pts1);
    IMP_NEW(npctransport::SitesPairScore, ps1,
            (interaction_range,
             interaction_k,
	     sigma0, sigma1,
             nonspecific_range,
             nonspecific_k,
             excluded_volume_k,
             sites0,
             sites1)
            );
   interaction_pair_scores_[interaction_id1] = ps1;
  }
  {
    core::ParticleTypes pts2;
    pts2.push_back(type1);
    pts2.push_back(type0);
    int interaction_id2 = get_ordered_type_pair_predicate()->get_value(pts2);
    IMP_NEW(npctransport::SitesPairScore, ps2,
            (interaction_range,
             interaction_k,
	     sigma0, sigma1,
             nonspecific_range,
             nonspecific_k,
             excluded_volume_k,
             sites1,
             sites0)
            );
    interaction_pair_scores_[interaction_id2] = ps2;
  }
}

boost::tuple<IMP::Restraint*, Object*>
Scoring::create_backbone_restraint
(double rest_length_factor,
 double backbone_k,
 ParticlesTemp beads,
 std::string name) const
{
  Pointer<Restraint> r(nullptr);
  Pointer<Object> bonds_score(nullptr); // for backward compatibility
  if(is_backbone_harmonic_){
    // Harmonic with relaxing spring
    IMP_NEW( HarmonicSpringSingletonScore,
             hsss,
             (20.0 * backbone_k, // TODO: 20* is arbitrary - just to make the particles stick well to the ends of the spring. Need to think how to set a value that would also allow the particles to pull the spring and vice versa realistically
              backbone_k) );
    ParticleIndexes pis;
    for(unsigned int i= 0; i<beads.size()-1; i++){
      pis.push_back(beads[i]->get_index());
      IMP_USAGE_CHECK(RelaxingSpring::get_is_setup(beads[i]),
		      "if backbone is harmonic spring, all chain beads should be"
		      " decorated with RelaxingSpring except for the tail");
    }
    // TODO: update rest length and backbone k of springs or leave as
    // is? For now, this is handled in FGChain
    IMP_NEW(IMP::container::ListSingletonContainer,
            springs,
            (get_model(), pis));
    r= container::create_restraint( hsss.get(),
                                    springs.get(),
                                    "Spring bonds " + name );
    bonds_score= hsss;
  } else {
    // Linear
    IMP_NEW(LinearWellPairScore, lwps,
            (rest_length_factor, backbone_k));
    IMP_NEW(IMP::container::ExclusiveConsecutivePairContainer,
	    bead_pairs,
	    ( get_model(),
	      IMP::get_indexes(beads),
	      "Bonds %1% " + name + " consecutive pairs" )
	    );
    r=  container::create_restraint( lwps.get(),
                                     bead_pairs.get(),
                                     "Bonds " + name  );
    bonds_score= lwps;
  }
  IMP_USAGE_CHECK(r != nullptr && bonds_score != nullptr,
                  "Backbone restraint is expected to be constructed by now");
  return boost::make_tuple(r.release(), bonds_score.release());
}

IMP::PairScore const*
Scoring::get_predicate_pair_score
( core::ParticleType t1, core::ParticleType t2) const
{
  // Get pair score
  core::ParticleTypes pts;
  pts.push_back(t1);  pts.push_back(t2);
  int iid = const_cast<core::OrderedTypePairPredicate*>
    (get_ordered_type_pair_predicate())->get_value(pts);
  t_map_pair_type_to_pair_score::const_iterator iter =
    interaction_pair_scores_.find(iid);
  if(iter == interaction_pair_scores_.end()){
    return nullptr;
  }
  return iter->second;
 }

IMP::PairScore*
Scoring::get_predicate_pair_score
( core::ParticleType t1, core::ParticleType t2)
{
  PairScore const* ret_value= ((Scoring const*)(this))->get_predicate_pair_score(t1, t2);
  return const_cast<PairScore*>(ret_value);
 }


// returns max of site-site and non-specific interaction ranges
// for interaction type it, if applicable
double
Scoring::get_interaction_range_for
  ( core::ParticleType t1, core::ParticleType t2,
     bool site_specific, bool non_specific) const
{
  IMP::PairScore const* ps = get_predicate_pair_score(t1, t2);
  if(ps==nullptr){
    return 0; // when not defined then only repulsive force upon touching
  }
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

boost::tuple<unsigned int,
             std::vector<unsigned int>,
             std::vector<unsigned int>,
             bool >
Scoring::get_site_interactions_statistics
( ParticleIndex pi1, ParticleIndex pi2) const
{
  core::Typed t1(get_model(), pi1);
  core::Typed t2(get_model(), pi2);
  core::ParticleType pt1=t1.get_type();
  core::ParticleType pt2=t2.get_type();
  SitesPairScore const* ps =
    dynamic_cast<SitesPairScore const*>
    (get_predicate_pair_score(pt1, pt2));
  if(ps==nullptr){
    // when not defined or not sites pair score then only repulsive force
    // upon touching
    std::vector<unsigned int> empty;
    return boost::make_tuple(0, empty, empty, false);
  }
  boost::tuple<unsigned int,
               std::vector<unsigned int>,
               std::vector<unsigned int>,
               bool> ret_value;
  ps->evaluate_site_contributions(get_model(),
                                  ParticleIndexPair(pi1, pi2),
                                  nullptr,
                                  &ret_value);
  return ret_value;
}



// add FG Nup chain with restraints
void
Scoring::add_chain_restraints(FGChain* chain)
{
  IMP_ALWAYS_CHECK( chain->get_rest_length_factor() > 0,
                    "non-positive rest length factor is not valid",
                    IMP::ValueException );
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
        if(get_sd()->get_is_fg_bead(_1) &&
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
              R += chain->get_chain_restraints(this);
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
      R += chain->get_chain_restraints(this);
    }
  return R;
}


void
Scoring::remove_particle_type
(core::ParticleType pt)
{
  // I. Clean bead_to_chain_map_:
  for( t_particle_index_to_fg_chain_map::iterator iter= bead_to_chain_map_.begin();
       iter != bead_to_chain_map_.end(); ) {
    IMP_USAGE_CHECK(core::Typed::get_is_setup(get_model(), iter->first),
                    "all particles are expected to be typed");
    core::ParticleType cur_type= core::Typed(get_model(), iter->first).get_type();
    if(pt == cur_type){
      bead_to_chain_map_.erase(iter++);
    } else {
      iter++;
    }
  }
  // II. Clean chains_set_:
  for( FGChainsSet::iterator iter= chains_set_.begin();
       iter != chains_set_.end(); ) {
    atom::Hierarchy h_root= (*iter)->get_root();
    IMP_USAGE_CHECK(core::Typed::get_is_setup(get_model(),
                                              h_root.get_particle_index()),
                    "all particles are expected to be typed");
    core::ParticleType cur_type= core::Typed(get_model(),
                                       h_root.get_particle_index()).get_type();
    if(pt == cur_type){
      chains_set_.erase(iter++);
    }else{
      iter++;
    }
  }
  // III. TODO: Clean interaction_pair_scors - actually not essential
  // IV. TODO: remove also for z-bias partricles map
  // V. TODO: remove also restrained anchor beads (complicated!)
  // VI. force referesh of scoring function
  get_scoring_function(true);
}

/***************************************************************/
/***************************** Creators ************************/
/***************************************************************/

// a close pair container for all unordered non-self pairs of specified beads and
// optimizable beads that include at least one optimizable bead
IMP::PairContainer*
Scoring::create_close_beads_container
( SingletonContainerAdaptor non_optimizable_beads,
  SingletonContainerAdaptor optimizable_beads) const
{
  using namespace container;
  IMP_LOG(TERSE,
          "Creating a close pair container for " << non_optimizable_beads->get_contents().size()
            << " non-optimizable beads and " << optimizable_beads->get_contents().size()
            << " optimizable ones" << std::endl);
  non_optimizable_beads.set_name_if_default(
      "CreateCloseBeadsContainerNonOptimizableBeadsInput%1%");
  optimizable_beads.set_name_if_default(
      "CreateCloseBeadsContainerOptimizableBeadsInput%1%");

  PairContainers pc_list;
  {
    // filter for excluding self pairs
    IMP_NEW( core::AllSamePairPredicate, aspp, () );
    // Filter for excluding consecutive chain beads (TODO: strongly tied to
    // implementation but then again, it is possibly important to
    // exclude bonded interactions - see FGChain.cpp. A possible
    // solution would be to add a static interface method in FGChain to
    // return a filter for bonded interactions):
    IMP_NEW( ExclusiveConsecutivePairFilter, ecpf, () ); // TODO: bug - might affect non-bonded beads that appear consecutively too
    // Pairs of opptimizable with optimizable:
    IMP_NEW(ClosePairContainer, cpc, // so range + 2*slack is what we get
            (optimizable_beads, get_range(), slack_) );
    cpc->add_pair_filter( aspp );
    cpc->add_pair_filter( ecpf );
    pc_list.push_back(cpc);
    // Pairs of non-optiomizable with optimizable if applicable:
    IMP_NEW(CloseBipartitePairContainer, cbpc, // so range + 2*slack is what we get
            (non_optimizable_beads, optimizable_beads, get_range(), slack_) );
    cbpc->add_pair_filter( aspp );
    cbpc->add_pair_filter( ecpf );
    pc_list.push_back(cbpc);
  }

  IMP_NEW(PairContainerSet, pcs, (pc_list));
  return pcs.release();
}


container::PredicatePairsRestraint
*Scoring::create_predicates_pair_restraint
(PairContainer* pair_container,
 bool is_attr_interactions_on) const
{
  // set linear repulsion upon penetration among all close pairs
  // returned by get_close_beads_container(), with different
  // scores for interactions between beads of different
  // (ordered) types
  IMP_NEW(container::PredicatePairsRestraint, predr,
          ( otpp_,
            pair_container ) );
  predr->set_is_get_inputs_ignores_individual_scores(true);
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
   Create bounding volume restraints such as box restraint and slab restraints,
   based on the box_is_on_, box_size_, slab_height_, slab_radius_, etc. class variables
   box_is_on_ == 1 : Box
   box_is_on_ == 2 : Sphere of volume box_size_**3

   Returns nullptr if neither a box nor a sphere
*/
Restraint* Scoring::create_bounding_volume_restraint
( SingletonContainerAdaptor beads) const
{
  beads.set_name_if_default("CreateBoundingBoxRestraintInput%1%");
  // Add bounding box restraint
  // TODO: replace backbone_spring_k_ with designated bounding volume
  //       force coefficient
  IMP_NEW(core::HarmonicUpperBound, hub, (0, get_excluded_volume_k()));
  IMP::Pointer<SingletonScore> ss( nullptr );
  if (get_has_bounding_box()){
    ss= new core::GenericBoundingBox3DSingletonScore<core::HarmonicUpperBound>
      (hub.get(), get_sd()->get_bounding_box());
  }
  else if (get_has_bounding_volume()){
    ss= new core::BoundingSphere3DSingletonScore
      (hub.get(), get_sd()->get_bounding_sphere());
  }
  else {
    IMP_ALWAYS_CHECK(!get_has_bounding_volume(),
                     "Invalid bounding volume (check box_is_on_ value if using a protobuf file)",
                     IMP::ValueException);
    return nullptr;
  }
  return container::create_restraint(ss.get(),
                                     beads.get(),
                                     "bounding box");
}

Restraint * Scoring::create_slab_restraint
( SingletonContainerAdaptor particles)  const
{
  particles.set_name_if_default("CreateSlabRestraintInput%1%");
  IMP::Pointer<IMP::PairScore> slab_score;
  if (get_sd()->get_is_slab_with_cylindrical_pore()) {
    slab_score =new SlabWithCylindricalPorePairScore
      (get_excluded_volume_k());
  } else {
    slab_score =new SlabWithToroidalPorePairScore
      (get_excluded_volume_k());
  }
  IMP::ParticlesTemp slab_particles;
  slab_particles.push_back( get_sd()->get_slab_particle() );
  IMP_NEW(container::AllBipartitePairContainer, abpc,
	  (slab_particles,
           particles.get()) );
  return container::create_restraint(slab_score.get(),
				     abpc.get(),
                                     "membrane slab");
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
    R = get_sd()->get_pore_radius(); // TODO: incompatible with dynamic pore radius
  }
  IMP_NEW(ZBiasSingletonScore, zbsc, (k, R) );
  return
    container::create_restraint(zbsc.get(),
                                ps.get(),
                                "z-bias potential");
}

/**
   add anchor bead to scoring function
*/
void
Scoring::add_restrained_anchor_bead
(IMP::Particle* p)
{
  Model* m= get_model();
  ParticleIndex pi= p->get_index();
  IMP_ALWAYS_CHECK(IMP::core::XYZ::get_is_setup(m, pi),
                       "anchor bead must be an initialized XYZ particle",
                       IMP::ValueException);
  if(get_sd()->get_is_slab_with_cylindrical_pore()){
    ParticleIndex pi_slab= get_sd()->get_slab_particle()->get_index();
    IMP_ALWAYS_CHECK(SlabWithCylindricalPore::get_is_setup(m, pi_slab),
                     "slab expected to have been setup as cylindrical pore"
                     " at this point",
                     IMP::ValueException);
    SlabWithCylindricalPore slab(m, pi_slab);
    core::XYZ xyz(m, pi);
    IMP_NEW(AnchorToCylidnricalPorePairScore, atcp_ps,
            (slab, xyz.get_coordinates(), get_sd()->get_pore_anchored_beads_k()) );
    IMP_NEW(container::ListSingletonContainer,
            lsc1,
            (m, IMP::ParticleIndexes(1,pi_slab))
            );
    IMP_NEW(container::ListSingletonContainer,
            lsc2,
            (m, IMP::ParticleIndexes(1,pi))
            );
    IMP_NEW(container::AllBipartitePairContainer,
            abpc,
            (lsc1.get(), lsc2.get())
            );
    Pointer<Restraint> r= create_restraint(atcp_ps.get(),
                                           abpc.get());
    anchor_restraints_.push_back(r);
  }
  else
    {
      // TODO: support toroidal case; also no pore case?
      IMP_NOT_IMPLEMENTED;
    }
}

Restraint*
Scoring::get_pore_radius_restraint() const
{
  Model* m= get_model();
  ParticleIndex pi_slab= get_sd()->get_slab_particle()->get_index();
  IMP_ALWAYS_CHECK(get_sd()->get_has_slab() &&
                   SlabWithPore::get_is_setup(m, pi_slab),
                   "Pore radius restraint should only be created when there's a pore",
                   ValueException);
  IMP_NEW(PoreRadiusSingletonScore,
          prss,
          (get_sd()->get_pore_radius(),
           get_sd()->get_pore_radius_k())
          );
  IMP_NEW(container::ListSingletonContainer,
          lsc,
          (m, ParticleIndexes(1, pi_slab)) );
  return create_restraint(prss.get(),
                          lsc.get());//                          "pore radius");
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
#undef GET_ASSIGNMENT_DEF
#undef GET_VALUE
#undef GET_VALUE_DEF
IMPNPCTRANSPORT_END_NAMESPACE
