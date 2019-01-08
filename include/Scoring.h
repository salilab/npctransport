/**
 *  \file npctransport/Scoring.h
 *  \brief scoring associated with a SimulationData object
 *
 *  Copyright 2007-2019 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_SCORING_H
#define IMPNPCTRANSPORT_SCORING_H

#include "npctransport_config.h"

#include "FGChain.h"
#include "npctransport_proto.fwd.h"
#include "Parameter.h"
// #include "SimulationData.h"

#include <IMP/Model.h>
#include <IMP/PairContainer.h>
#include <IMP/ScoringFunction.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/Pointer.h>
#include <IMP/WeakPointer.h>
#include <boost/unordered_map.hpp>
#include <IMP/container/CloseBipartitePairContainer.h>
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/container/PairContainerSet.h>
#include <IMP/container/PredicatePairsRestraint.h>
#include <IMP/core/pair_predicates.h>
#include <IMP/core/BoundingBox3DSingletonScore.h>
#include <IMP/core/Typed.h>
#include <IMP/core/RestraintsScoringFunction.h>
#include "boost/tuple/tuple.hpp"

#include <boost/timer.hpp>
#include <string>


IMPNPCTRANSPORT_BEGIN_NAMESPACE

class SimulationData;
class LinearWellPairScore;
class HarmonicSpringSingletonScore;
//class FGChain;

//! Scoring associated with a SimulationData object
class IMPNPCTRANSPORTEXPORT Scoring: public Object
{
 private:
  Parameter<int> box_is_on_;
  Parameter<double> interaction_k_;
  Parameter<double> interaction_range_;
  Parameter<double> backbone_k_;
  Parameter<bool> is_backbone_harmonic_;
  Parameter<double> slack_;
  Parameter<double> nonspecific_k_;
  Parameter<double> nonspecific_range_;
  Parameter<double> excluded_volume_k_;
  Parameter<double> range_;
  // Per type scaling factors for the interaction parameters
  boost::unordered_map<core::ParticleType, double> interaction_range_factors_;
  boost::unordered_map<core::ParticleType, double> interaction_k_factors_;

 private:

  // the sim-data that uses this scoring object
  UncheckedWeakPointer<SimulationData> owner_sd_;

  /***************** Cache only variables ************/

  // see get_close_beads_container()
  PointerMember
    <IMP::PairContainer> close_beads_container_;

  // generates hash values ('predicates') for ordered types pairs
  // (e.g., pairs of ParticleTypes)
  PointerMember
    <IMP::core::OrderedTypePairPredicate> otpp_;

  PointerMember
    <IMP::core::RestraintsScoringFunction> scoring_function_;
  IMP::Restraints scoring_function_rs_; // the restraints associated with the last created scoring_function_, valid only when it is not null

  // contains all restraints between pairs of particle types
  PointerMember
    <container::PredicatePairsRestraint> predr_;

  typedef  boost::unordered_map< int, PointerMember<IMP::PairScore > >
    t_map_pair_type_to_pair_score;
  // scores between particle, mapped by their interaction id,
  // where interaction id is from OrderedTypePairPredicate otpp_
  // applied on each interacting particle type
  t_map_pair_type_to_pair_score
    interaction_pair_scores_;

  PointerMember
    <IMP::Restraint> box_restraint_;

  PointerMember
    <IMP::Restraint> slab_restraint_;

  // scores on anchor particles
  IMP::Restraints anchor_restraints_;

  // pore radius restraint
  PointerMember
    <IMP::Restraint> pore_radius_restraint_;

  // a map from particle indexes of chain identifiers to corresponding
  // chain backbone restraints (see chain_ids_map)
  typedef boost::unordered_set
    < PointerMember<FGChain> > FGChainsSet;
  FGChainsSet chains_set_;
  // a map from particle indexes of chain beads to the chain objects that they
  // represent (TODO: probably can be obliterated by adding appropriate
  // decorators)
  typedef boost::unordered_map
    < ParticleIndex, PointerMember<FGChain> >
    t_particle_index_to_fg_chain_map;
  t_particle_index_to_fg_chain_map bead_to_chain_map_;

  // particles to be z-biased on call to get_z_bias_restraints()
  // with key being the k value of each particle subset
  typedef boost::unordered_map<double, ParticlesTemp> t_z_bias_particles_map;
  t_z_bias_particles_map z_bias_particles_map_;

  // restraints for z-biasing particles (plural for different k force constants)
  IMP::Restraints z_bias_restraints_;

  // custom restraints that are added to the scoring function
  IMP::Restraints custom_restraints_;


 public:
  /**
     @param[in]  owner_sd npc simulation data object that owns this scoring
                          object
     @param[out] data protobuf data with parameters for
                      this assignment, used to parameterize
                      the scoring function
   */
  Scoring(SimulationData* owner_sd,
          const ::npctransport_proto::Assignment &data);


  /**
     returns a scoring function for the NPC simulation, based on:
     1) the chain restraints added by add_chain_restraints()
     2) default repulsive interactions on particles in get_beads()
     3) attractive interactions that depend on pair types, as added by
        add_interaction()
     4) the slab or bounding box restraints if relevant on all particles
        in get_beads()
     5) z-biasing potential if add_z_bias_restraint() was called for any
        particles

     @param if force_update is true, forces recreation of the scoring function,
                   o/w cached version may be retrieved if available, which may not
                   be completely up to date.

     @note if force_udpate=false, might fail to include e.g., new particles
           or new interactions that were added after the last call
   */
  IMP::ScoringFunction* get_scoring_function(bool force_update=false);

  /**
     returns the restrains associated with the scoring function return by
     get_scoring_function(force_update)
  */
  IMP::Restraints
    get_scoring_function_restraints(bool force_update=false);


  /**
     returns a custom scoring function for the NPC simulation, based on:
     1) the chain restraints added by add_chain_restraints()
        that overlap with 'particles'
     2) default repulsive interactions between particles in particles
        and optimizable particles
     3) attractive interactions that depend on pair types, as added by
        add_interaction()
     4) the slab or bounding box restraints if relevant on all particles
        in 'particles'

        @param extra_restraints  ad-hoc restraints to be added to scoring function
        @param non_optimizable_beads particles container on which to apply bounding volume
                         and pair constraints only if they also involve optimizable beads
        @param optimizable_beads particles container on which to apply bounding volume
                         and pair contstraints, either with other optimizable beads or with
                         non-optimizable ones.
        @param is_attr_interactions_on if false, only repulsive interactions will be
                                       computed between pairs of particles
        @note it is assumed the two sets of beads are disjoint
  */
  IMP::ScoringFunction*
    get_custom_scoring_function
    ( const RestraintsTemp& extra_restraints,
      SingletonContainerAdaptor non_optimizable_beads,
      SingletonContainerAdaptor optimizable_beads,
      bool is_attr_interactions_on = true ) const;

  //! same as other get_custom_scoring_function variant but
  //! based on a list of particle indexes and allowing an empty vector
  //! of non-optimizable beads as input
  //! @note  The optimizable beads are still
  //! treated dynamically (e.g. if they're associated with a dynamic container)
  IMP::ScoringFunction*
    get_custom_scoring_function
    ( const RestraintsTemp& extra_restraints,
      ParticleIndexes non_optimizable_beads,
      SingletonContainerAdaptor optimizable_beads,
      bool is_attr_interactions_on = true ) const
    {
      IMP_NEW(IMP::container::ListSingletonContainer,
              nob_lsc,
              (get_model(), non_optimizable_beads));
      return get_custom_scoring_function(extra_restraints,
                                         nob_lsc,
                                         optimizable_beads,
                                         is_attr_interactions_on);
    }


  /**
     a pair container that was used to define interaction particles
     that will be used in the next call to get_scoring_function(false)
     (so manipulating it might affect the scoring function)

     @param update if true, forces recreation of the cached container,
                   o/w cached version that was used in last call to
                   get_scoring_function() or to this function will
                   be used

     @return a container with all pairs of particles (a,b) such that:
     1) a is an optimizable bead particle and b is any bead
        particle (static or not).
     2) a and b are close (sphere surfaces within range get_range()+some slack)
     3) a and b do not appear consecutively within the model (e.g., bonded fg motifs)
     4) a and b are not the same particles
     5) b and a do not appear on the list (= unordered list)

     @note if udpate=false, might fail to include e.g., new particles
           or new interactions that were added after the last call
     @note TODO: right now will not return any consecutive particles - that's
                 erroneous (e.g., last particle of one chain and first particle
                 of next chain) though may be negligible in practice
     @note supposed to be robust to dynamic changes to the beads list,
           though need to double check (TODO)
  */
  IMP::PairContainer *get_close_beads_container(bool update=false);

  /**
     Returns the restraints over pairs of particles based on their type,
     which will be used in the next call to get_scoring_function(false)
     (so manipulating it might affect the scoring function). Also, if update
     is false, then it is guaranteed that these are the same restraints
     used by the last call to get_scoring_function()

     Different scores are used for particles of different (ordered)
     particle types.  When called for the first time, returns a new
     PredicatePairsRestraints over all diffusing particles and sets a
     default linear repulsion restraint between all pairs returned by
     get_close_beads_container()

     @param update if true, forces recreation of the cached container,
                   o/w cached version that was used in last call to
                   get_scoring_function(), if applicable

     @note if udpate=false, might fail to include e.g., new particles
           or new interactions that were added to the simulation
           after the last call
  */
  container::PredicatePairsRestraint *get_predicates_pair_restraint
    (bool update=false);


  /** returns the box/sphere restraint on >get_sd()->get_beads()
     which will be used in the next call to get_scoring_function(false)
     (so manipulating it might affect the scoring function). If update
     is false, then it is guaranteed that these are the same restraints
     used by the last call to get_scoring_function.

     @param update if true, forces recreation of the cached container,
                   o/w cached version that was used in last call to
                   get_scoring_function(), if applicable

     @note if udpate=false, might fail to include e.g., new particles
           that were added after the last call

     @see get_has_bounding_box(), get_has_bounding_sphere(),
          get_has_bounding_volume()
  */
  Restraint *get_bounding_volume_restraint(bool update=false);


  /** returns the slab restraint on >get_sd()->get_beads()
     which will be used in the next call to get_scoring_function(false)
     (so manipulating it might affect the scoring function). Also, if update
     is false, then it is guaranteed that these are the same restraints
     used by the last call to get_scoring_function()

     @param update if true, forces recreation of the cached container,
                   o/w cached version that was used in last call to
                   get_scoring_function(), if applicable

     @note if udpate=false, might fail to include e.g., new particles
           that were added after the last call
  */
  Restraint *get_slab_restraint(bool update=false);


  // swig doesn't equate the two protobuf types
  /**
     add a SitesPairScore restraint that applies to particles of
     types t0 and t1 (specified in idata) to the PredicatePairsRestraint
     object that is returned by get_predr().
     The restraint is symmetric, and applies to both (t0,t1) and (t1,t0).

     A SitesPairScore restraint means site-specific
     attractive forces between surface binding sites on each particle +
     non-specific attraction and soft repulsion between the entire particles.

     @param idata the protobuf data about the interaction (particle types,
            interaction coefficients, etc.)

     \see SitesPairScore
  */
  void add_interaction
    ( const ::npctransport_proto::Assignment_InteractionAssignment &idata);





  /***************************************************************/
  /***************************** Creators ************************/
  /***************************************************************/

  /**
     creates a new pair container that returns all pairs of particles (a,b) such that:
     1) a is any particle from 'non_optimizable_beads' or 'optimizable beads',
        and b is any particle in 'optimizable_beads'
     2) a and b are close (sphere surfaces within range get_range())
     3) a and b do not appear consecutively within the model (e.g., fg repeats)
     4) a is not b
     5) b and a do not appear in the list (effectively, a list of unordered pairs)

     @param non_optimizable_beads a container for particles that are assumed to be non-optimizable
                                  and therefore are not paired among themselves
     @param optimizable_beads A container for particles that are assumed to be optimizable, and may
                              be paired among themselves or with non_optimizable_beads
     @note it is assumed the two sets of beads are disjoint
  */
  IMP::PairContainer *create_close_beads_container
    ( SingletonContainerAdaptor non_optimiziable_beads,
      SingletonContainerAdaptor optimizable_beads)
    const;

  //! same as other create_close_beads_container variant but
  //! based on a list of particle indexes and allowing an empty vector
  //! of non-optimizable beads as input. The optimizable beads are still
  //! treated dynamically (e.g. if they're associated with a dynamic container)
  IMP::PairContainer *create_close_beads_container
    ( ParticleIndexes non_optimizable_beads,
      SingletonContainerAdaptor optimizable_beads)
    const
    {
      IMP_NEW(IMP::container::ListSingletonContainer,
              nob_lsc,
              (get_model(), non_optimizable_beads));
      return create_close_beads_container(nob_lsc,
                                          optimizable_beads);
    }


  /**
     Creates a new container for restraints over pairs of beads. Different
     scores are used for beads of different (ordered) particle types
     based on interaction_pair_scores, or just a default linear repulsive
     with k=get_excluded_volume_k() force on all other pairs in bead_pairs.

     @param bead_pairs a container for the pairs of beads to be restrained
     @param is_attr_interactions_on whether to include attractive interactions
                                    that were added by add_interaction()
  */
  container::PredicatePairsRestraint *create_predicates_pair_restraint
    ( PairContainer* bead_pairs,
      bool is_attr_interactions_on = true) const;

  /**
     Creates bounding volume (box or sphere) restraint based on the
     box_size_ class variable and get_has_bounding_sphere() and
     get_has_bounding_box() variables (typically based on box_is_on_),
     , and apply it to all beads returned by get_sd()->get_beads()

     @param beads beads on which to apply the constraint

     @note this restraint is (supposed to) guarantee to be always
           updated with the list inside the container

     @return a newly created box restraint
  */
  Restraint* create_bounding_volume_restraint
    ( SingletonContainerAdaptor beads ) const;

  /**
     Creates slab bounding volume restraint, based on the slab_thickness_
     and tunnel_radius_ class variables, and apply it to all beads
     returned by get_sd()->get_beads()

     @param beads beads on which to apply the constraint

     @note this restraint is (supposed to) guarantee to be always
           updated with the list inside the container

     @return a newly created slab restraint
   */
  Restraint* create_slab_restraint
    ( SingletonContainerAdaptor beads ) const;


  /************************************************************/
  /************* various simple getters and setters *******************/
  /************************************************************/

  /** returns the model associated with the owned SimulationData */
  Model* get_model();

#ifndef SWIG
  /** returns the model associated with the owned SimulationData */
  Model* get_model() const;
#endif


  /** return the SimulationData object that owns this ScoringFunction */
  SimulationData* get_sd() {
    return owner_sd_;
  }

#ifndef SWIG
  /** return the SimulationData object that owns this ScoringFunction */
  SimulationData const* get_sd() const{
    return owner_sd_;
  }
#endif

    //! returns true if simulation has a bounding simulation box
  bool get_has_bounding_box() const {
    return box_is_on_==1;
  }

  //! returns true if simulation has a bounding simulation sphere ('cell-like')
  bool get_has_bounding_sphere() const {
    return box_is_on_==2;
  }

  //! returns true if simulation has any bounding volume (box or sphere are supported)
  bool get_has_bounding_volume() const {
    return get_has_bounding_box() || get_has_bounding_sphere();
  }

  /** returns the default spring constant between consecutive beads
      in a chain of beads - this is used for chains whose k
      is non-positive
  */
  double get_default_backbone_k() const { return backbone_k_; }

  bool get_is_backbone_harmonic() const { return is_backbone_harmonic_; }

#ifndef SWIG
  //! Create a backbone bond restraint over beads according to class flags
  /** create a backbone bond restraint with specified rest_length_factor
      (relative to sum of spheres radii) and backbone k over beads. The restraint score
      type is according to the class flags - either a linear or
      harmonic potential, with k in units of kcal/mol/A or kcal/mol/A^2,
      respectively, with the harmonic potential corresponding to a relaxing
      spring.

      returns both the restraint and the linear well pair score or harmonic
      singleton score that is used in it (which one depends on the IMP version,
      assigns nullptr to the other), as a pointer to a generic object (= common
      ancestor class)
  */
  boost::tuple< Restraint*, Object*>
    create_backbone_restraint
    ( double rest_length_factor,
      double backbone_k,
      ParticlesTemp beads,
      std::string name) const;
#endif

  /** returns the constant force applied by overlapping spheres on each
      other
  */
  double get_excluded_volume_k() const { return excluded_volume_k_; }

  /** returns the default interaction constant, which is used
      only if interaction_k is not specified for a certain interaction */
  double get_interaction_k() const { return interaction_k_; }

  //! returns the pair score for particles of types t1 and t2   */
  IMP::PairScore const* get_predicate_pair_score
    ( core::ParticleType t1, core::ParticleType t2) const;

  //! returns the pair score for particles of types t1 and t2   */
  IMP::PairScore* get_predicate_pair_score
    ( core::ParticleType t1, core::ParticleType t2);


  /** @return maximal actual interaction range (normalized by
      range factors) for a given interaction type, as set by
      add_interaction(), over site-site and/or non-specific
      interactions if applicable.  If not applicable, returns -1.0

      @param t1,t2 the types of beads in this interaction
      @param site_specific include site-site ranges
      @param non_specific include non-specific sphere interaction ranges
      @note if both site_specific and non_specific are true, takes a
            maximum over both.
  */
  double get_interaction_range_for
    (core::ParticleType t1, core::ParticleType t2,
     bool site_specific = true,
     bool non_specific = false) const;

#ifndef SWIG
  //! return various statistics about site-site interactions
  //! between pi1 and pi2
  /**
     return a tuple with:
     a - number of site-site interaction that contribute
        to the score between pi1 and pi2 in the model, based on their
         current coordinates and their particle types
     b - vector of occupied pi1 sites counts (count number of interactions per site)
     c - vector of occupied pi2 sites counts (count number of interactions per site)
  */
  boost::tuple< unsigned int, std::vector<unsigned int>, std::vector<unsigned int> >
    get_site_interactions_statistics
    ( ParticleIndex pi1, ParticleIndex pi2) const;

  //! return various statistics about site-site interactions
  //! between p1 and p2
  /**
     return a tuple with:
     a) number of site-site interaction that contribute
        to the score between p1 and p2 in the model, based on their
         current coordinates and their particle types
     b) vector of occupied pi1 sites counts (count number of interactions per site)
     c) vector of occupied pi2 sites counts (count number of interactions per site)
  */
  boost::tuple< unsigned int, std::vector<unsigned int>, std::vector<unsigned int> >
    get_site_interactions_statistics
    ( Particle const* p1, Particle const* p2) const
  {
    return
      get_site_interactions_statistics
      (p1->get_index(), p2->get_index());
  }
#endif


  core::OrderedTypePairPredicate* get_ordered_type_pair_predicate()
  { return otpp_; }

#ifndef SWIG
  core::OrderedTypePairPredicate const* get_ordered_type_pair_predicate() const
  { return otpp_; }
#endif

  /** add FG chain to this scoring object, using the
      chain->get_chain_restraints() method to retrieve the restraints
      upon a call to get_scoring_function() or create_scoring_function(),
      etc.

      @param chain an FG chain associated with internal scoring

      @return pointer to the newly created restraint

      @note if the chain->get_backbone_k() is non-positive, it is reset to
            this->get_default_backbone_k()

      @note this method assumes that all such chains will be disjoint
      and so it is later possible to use
      container::ExclusiveConsecutivePairFilter to filter out all
      pairs of particles connected by such chain restraints.
  */
  void add_chain_restraints(FGChain* chain);


  /**
      returns restraints on all fg nup chains that were added by
      add_chain_restraints() so far and contain any bead in 'beads'.
      The restraints are updated directly from the original chain objects
      chain->get_chain_restraints() object.

      @param beads a least of chain fine beads over which restraints are
                   collected

      @return a list of all chain restraints that have been added so far
              and contain beads in 'beads'
  */
  Restraints get_chain_restraints_on
    ( SingletonContainerAdaptor beads ) const;

  /**
      @return all restraints of fg nup chains that were added by
      add_chain_restraints() so far
  */
  Restraints get_all_chain_restraints() const;

  //! remove all scoring function elements related to particle type pt
  //! NOTE: a few aspects of the scoring function do not support this function
  //!       as of 2018-06-19 - namely anchor restraints and z coordinate bias
  void remove_particle_type(core::ParticleType pt);

  /**
     adds a biasing potential for particles ps towards z, to the restraints
     returned by get_z_bias_restraints(), and for the scoring function
     generated by get_scoring_function()

     @param ps container for the particles to be biased
     @param k force constant

     \see get_z_bias_restraints()
     \see get_scoring_function()
  */
  void add_z_bias_restraint(SingletonContainerAdaptor ps, double k);

  /**
     adds a biasing potential for particle p towards z, to the restraints
     returned by get_z_bias_restraints(), and for the scoring function
     generated by get_scoring_function()


     @param p the particle to be biased
     @param k force constant
  */
  void add_z_bias_restraint(Particle* p, double k);

  /**
     return z_bias restraints for all particles that were added using
     add_z_bias_restraint()
  */
  IMP::Restraints get_z_bias_restraints();

  /**
      create a z bias restraint pulling along z-axis for the particles
      in ps. If a slab is present, the force acts only within tunnel radius
      distance from its central axis.

      @param ps the particles container
      @param k force constant for biasing ps
  */
  IMP::Restraint* create_z_bias_restraint(SingletonContainerAdaptor ps,
                                          double k) const;

  void add_custom_restraint(IMP::Restraint* r)
  { custom_restraints_.push_back(r); }

  void add_custom_restraints(IMP::Restraints R)
  { custom_restraints_ += R; }

  void clear_custom_restraints()
  { custom_restraints_.clear(); }

  IMP::Restraints get_custom_restraints() const
    { return custom_restraints_; }

  /**
     add anchor bead to scoring function, using current position of p
     as equilibrium position relative to pore

     @param p anchored particle

     @note currently only cylindrical pore case is supported
   */
  void add_restrained_anchor_bead(IMP::Particle* p);

  Restraint* get_pore_radius_restraint() const;

#ifndef SWIG
  /**
     returns pointers to the collection of the FG Chains stored in
     this scoring object (can be used to e.g. scale scoring
     information up or down during optimization)
   */
  FGChains get_fg_chains() {
    return FGChains(chains_set_.begin(), chains_set_.end());
  }
#endif

  double get_range() const { return range_; }

  /**
      sets the multiplicative scaling of the interaction range for
      pair interactions involving a particle of this type, for all
      future calls to add_interaction()
  */
  void set_interaction_range_factor
    ( IMP::core::ParticleType type, double value )
    {
      IMP_ALWAYS_CHECK(value>0, "interaction_range_factor must be positive",
                       IMP::ValueException);
      interaction_range_factors_[type] = value;
    }

  /**
      sets the multiplicative scaling of the interaction k constants for
      pair interactions involving a particle of this type, for all future
      calls to add_interaction()
  */
  void set_interaction_k_factor
    ( IMP::core::ParticleType type, double value )
    {
      IMP_ALWAYS_CHECK(value>0, "interaction_k_factor must be positive",
                       IMP::ValueException);
      interaction_k_factors_[type] = value;
    }

 private:

 public:
  IMP_OBJECT_METHODS(Scoring);
};




IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SCORING_H */
