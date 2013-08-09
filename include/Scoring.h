/**
 *  \file npctransport/Scoring.h
 *  \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_SCORING_H
#define IMPNPCTRANSPORT_SCORING_H

#include "npctransport_config.h"
#include <IMP/Model.h>
#include <IMP/PairContainer.h>
#include <IMP/ScoringFunction.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/base/Pointer.h>
#include <IMP/base/WeakPointer.h>
#include <IMP/base/map.h>
#include <IMP/container/CloseBipartitePairContainer.h>
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/container/PairContainerSet.h>
#include <IMP/container/PredicatePairsRestraint.h>
#include <IMP/core/pair_predicates.h>
#include <IMP/core/BoundingBox3DSingletonScore.h>
#include <IMP/core/Typed.h>
#include <IMP/core/RestraintsScoringFunction.h>
#include "linear_distance_pair_scores.h"
#include "Parameter.h"
//#include "SimulationData.h"
#include "SlabSingletonScore.h"

#include <boost/timer.hpp>
#include <string>

#ifndef SWIG
namespace npctransport_proto {
  class Assignment;
  class Assignment_FGAssignment;
  class Assignment_InteractionAssignment;
  class Assignment_FloaterAssignment;
  class Assignment_ObstacleAssignment;
}
#endif

IMPNPCTRANSPORT_BEGIN_NAMESPACE

class SimulationData;

class IMPNPCTRANSPORTEXPORT Scoring: public base::Object
{
 private:
  Parameter<double> box_side_;
  Parameter<double> tunnel_radius_;
  Parameter<double> slab_thickness_;
  Parameter<bool> box_is_on_;
  Parameter<bool> slab_is_on_;
  Parameter<double> interaction_k_;
  Parameter<double> interaction_range_;
  Parameter<double> backbone_k_;
  Parameter<double> slack_;
  Parameter<double> nonspecific_k_;
  Parameter<double> nonspecific_range_;
  Parameter<double> excluded_volume_k_;
  Parameter<double> range_;
  // Per type scaling factors for the interaction parameters
  base::map<core::ParticleType, double> interaction_range_factors_;
  base::map<core::ParticleType, double> interaction_k_factors_;

 private:

  base::UncheckedWeakPointer<SimulationData> owner_sd_;

  // true if the list of particles in the owner_sd has changed
  // and in the middle of changing it
  //  bool is_updating_particles_;

  /***************** Cache only variables ************/

  // see get_close_diffusers_container()
  base::PointerMember
    <IMP::PairContainer> close_diffusers_container_;

  // generates hash values ('predicates') for ordered types pairs
  // (e.g., pairs of ParticleTypes)
  base::PointerMember
    <IMP::core::OrderedTypePairPredicate> otpp_;

  base::PointerMember
    <IMP::core::RestraintsScoringFunction> scoring_function_;

  // contains all restraints between pairs of particle types
  base::PointerMember
    <container::PredicatePairsRestraint> predr_;

  typedef  IMP::base::map< int, base::PointerMember< IMP::PairScore > >
    map_pair_type_to_pair_score;
  // scores between particle, mapped by their interaction id,
  // where interaction id is from OrderedTypePairPredicate otpp_
  // applied on each interacting particle type
  map_pair_type_to_pair_score
    interaction_pair_scores_;

  base::PointerMember
    <IMP::Restraint> box_restraint_;

  base::PointerMember
    <IMP::Restraint> slab_restraint_;

  LinearWellPairScores chain_scores_;

  // a map from particle indexes of chain parents to corresponding
  // chain backbone restraints
  typedef IMP::base::map< ParticleIndex, IMP::base::PointerMember<Restraint> >
    t_chain_restraints_map;
  t_chain_restraints_map chain_restraints_map_;

 public:
  /**
     @param[in]  owner_sd npc simulation data object that owns this scoring
                          object
     @param[out] data protobuf data with parameters for
                      this assignment, used to parametrize
                      the scoring function
   */
  Scoring(SimulationData* owner_sd,
          const ::npctransport_proto::Assignment &data);


  /**
     returns a scoring function for the NPC simulation, based on:
     1) the chain restraints added by add_chain_restraint()
     2) default repulsive interactions on particles in get_diffusers()
     3) attractive interactions that depend on pair types, as added by
        add_interaction()
     4) the slab or bounding box restraints if relevant on all particles
        in get_diffusers()

     @param update if true, forces recreation of the scoring function,
                   o/w cached version may be retrieved, which may not
                   be completely up to date

     @note if udpate=false, might fail to include e.g., new particles
           or new interactions that were added after the last call
   */
  IMP::ScoringFunction* get_scoring_function(bool update=false);

  /**
     returns a custom scoring function for the NPC simulation, based on:
     1) the chain restraints added by add_chain_restraint()
        that overlap with 'particles'
     2) default repulsive interactions between particles in particles
        nad optimizable particles
     3) attractive interactions that depend on pair types, as added by
        add_interaction()
     4) the slab or bounding box restraints if relevant on all particles
        in 'particles'

        @param extra_restraints  ad-hoc restraints to be added to scoring function
        @param particles particles container on which to apply bounding volume
                         and pair constraints
        @param optimizable_particles interaction scores (both repulsive or attractive)
                                     are computed only for interactions that involve
                                     these optimizable particles.
        @param is_attr_interactions_on if false, only repulsive interactions will be
                                       computed between pairs of particles
  */
  IMP::ScoringFunction*
    get_custom_scoring_function
    ( const RestraintsTemp& extra_restraints,
      SingletonContainerAdaptor particles,
      SingletonContainerAdaptor optimizable_particles,
      bool is_attr_interactions_on = true ) const;

  /** Update the scoring function that the list of particles,
      to be retrieved from owning simulation data, has changed
   */
  //  void update_particles() const{
  //    is_updating_particles_ = true;
  //  }


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
  /**
     a pair container that was used to define interaction particles
     that will be used in the next call to get_scoring_function(false)
     (so manipulating it might affect the scoring function)

     @param update if true, forces recreation of the cached container,
                   o/w cached version that was used in last call to
                   get_scoring_function() or to this function will
                   be used

     @return a container with all pairs of particles (a,b) such that:
     1) a is an optimizable diffuser particle and b is any diffuser
        particle (static or not).
     2) a and b are close (sphere surfaces within range get_range())
     3) a and b do not appear consecutively within the model (e.g., fg repeats)

     @note if udpate=false, might fail to include e.g., new particles
           or new interactions that were added after the last call
     @note TODO: right now will not return any consecutive particles - that's
                 erroneous (e.g., last particle of one chain and first particle
                 of next chain) though may be negligible in practice
     @note supposed to be robust to dynamic changes to the diffusers list,
           though need to double check (TODO)
  */
  IMP::PairContainer *get_close_diffusers_container(bool update=false);

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
     get_close_diffusers_container()

     @param update if true, forces recreation of the cached container,
                   o/w cached version that was used in last call to
                   get_scoring_function(), if applicable

     @note if udpate=false, might fail to include e.g., new particles
           or new interactions that were added to the simulation
           after the last call
  */
  container::PredicatePairsRestraint *get_predicates_pair_restraint
    (bool update=false);


  /** returns the box restraint on >get_sd()->get_diffusers()
     which will be used in the next call to get_scoring_function(false)
     (so manipulating it might affect the scoring function). Also, if update
     is false, then it is guaranteed that these are the same restraints
     used by the last call to get_scoring_function

     @param update if true, forces recreation of the cached container,
                   o/w cached version that was used in last call to
                   get_scoring_function(), if applicable

     @note if udpate=false, might fail to include e.g., new particles
           that were added after the last call
  */
  Restraint *get_bounding_box_restraint(bool update=false);


  /** returns the slab restraint on >get_sd()->get_diffusers()
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
#ifndef SWIG
  /**
     add a SitesPairScore restraint that applies to particles of
     types t0 and t1 (specified in idata) to the PredicatePairsRestraint
     object that is returned by get_predr().
     The restraint is symmetric, and applies to both (t0,t1) and (t1,t0).

     A SitesPairScore restraint means site-specific
     attractive forces between surface bidning sites on each particle +
     non-specific attraction and soft repulsion between the entire particles.

     @param idata the protobuf data about the interaction (particle types,
            interaction coefficients, etc.)

     \see SitesPairScore
  */
  void add_interaction
    ( const ::npctransport_proto::Assignment_InteractionAssignment &idata);
#endif



  /***************************************************************/
  /***************************** Creators ************************/
  /***************************************************************/

  /**
     creates a new pair container that returns all pairs of particles (a,b) such that:
     1) a is any particle from 'particles' and b is any optimizable particle
        in 'optimizable_particles' (static or not).
     2) a and b are close (sphere surfaces within range get_range())
     3) a and b do not appear consecutively within the model (e.g., fg repeats)

     @param particles a container For particles over which the return value works
     @param optimizable_particles A container for particles that are also optimizable.

     @note TODO: right now will not return any consecutive particles - that's
                 erroneous (e.g., last particle of one chain and first particle
                 of next chain) though may be negligible in practice
  */
  IMP::PairContainer *create_close_diffusers_container
    ( SingletonContainerAdaptor particles,
      SingletonContainerAdaptor optimizable_particles)
    const;

  /**
     Creates a new container for restraints over pairs of particles. Different
     scores are used for particles of different (ordered) particle types
     based on interaction_pair_scores, or just a default linear repulsive
     with k=get_excluded_volume_k() force on all other particle pairs in particle_pairs.

     @param particle_pairs a container for the pairs of particles to be restrained
     @param is_attr_interactions_on whether to include attractive interactions
                                    that were added by add_interaction()
  */
  container::PredicatePairsRestraint *create_predicates_pair_restraint
    ( PairContainer* particle_pairs,
      bool is_attr_interactions_on = true) const;

  /**
     Creates bounding box restraint based on the box_size_
     class variable, and apply it to all diffusers returned
     by get_sd()->get_diffusers()

     @param particles particles on which to apply the constraint

     @note this restraint is (supposed to) gurantee to be always
           updated with the list inside the container

     @return a newly created box restraint
  */
  Restraint* create_bounding_box_restraint
    ( SingletonContainerAdaptor particles ) const;

  /**
     Creates slab bounding volume restraint, based on the slab_thickness_
     and tunnel_radius_ class variables, and apply it to all diffusers
     returned by get_sd()->get_diffusers()

     @param particles particles on which to apply the constraint

     @note this restraint is (supposed to) gurantee to be always
           updated with the list inside the container

     @return a newly created slab restraint
   */
  Restraint* create_slab_restraint
    ( SingletonContainerAdaptor particles ) const;


  /************************************************************/
  /************* various simple getters and setters *******************/
  /************************************************************/

  /** returns the model associated with the owned SimulationData */
  Model* get_model();

#ifndef SWIG
  /** returns the model associated with the owned SimulationData */
  Model* get_model() const;
#endif

  // returns true if a bounding box restraint is defined */
  bool get_has_bounding_box() const
  { return box_is_on_; }

  /** returns true if a slab restraint is defined */
  bool get_has_slab() const
  { return slab_is_on_; }

  double get_backbone_k() const { return backbone_k_; }

  double get_excluded_volume_k() const { return excluded_volume_k_; }

  double get_interaction_k() const { return interaction_k_; }

  core::OrderedTypePairPredicate* get_ordered_type_pair_predicate()
  { return otpp_; }

#ifndef SWIG
  core::OrderedTypePairPredicate const* get_ordered_type_pair_predicate() const
  { return otpp_; }
#endif

  /** create an fg chain restraint on consecutive chain particles
      and store it in the internal list of chain restraint, which
      will be used in future calls to get_scoring_function(true) or
      create_scoring_function()

      @note this method assumes that all such chains will be disjoint
      and so it is later possible to use
      container::ExclusiveConsecutivePairFilter to filter out all
      pairs of particles connected by such chain restraints.

      @param chain_root the root of the chain whose particles are restrained
      @param rest_length the rest length of consecutive chain beads
      @param name the name of the chain (to be used for naming the restraint

      @return pointer to the newly created restraint
  */
  Restraint* add_chain_restraint(IMP::atom::Hierarchy chain_root,
                                 double rest_length,
                                 std::string name);


  /**
      returns all restraints on fg nup chains that were added by
      add_chain_restraint() so far and pertain to particles in 'particles'

      @param particles a list of diffusing particles, such that every returned restraint
                       must pertain to a particle in the list

      @return pointers to all chain restraints that have been added so far
              and pertain to particle in 'particles' list
  */
  Restraints get_chain_restraints_on
    ( SingletonContainerAdaptor particles ) const;

  /**
      @return all restraints on fg nup chains that were added by
      add_chain_restraint() so far
  */
  Restraints get_all_chain_restraints() const {
    Restraints rs;
    std::transform( chain_restraints_map_.begin(),
                    chain_restraints_map_.end(),
                    std::back_inserter(rs),
                    boost::bind
                    ( &t_chain_restraints_map::value_type::second, _1)
                  );
    return rs;
  }


  /**
     returns a reference to the collection of score functions for FG backbones
     (can be used to e.g. scale them up or down during optimization)
   */
  LinearWellPairScores get_chain_scores() { return chain_scores_; }

  double get_slab_thickness() const { return slab_thickness_; }

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
                       IMP::base::ValueException);
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
                       IMP::base::ValueException);
      interaction_k_factors_[type] = value;
    }

 private:

 public:
  IMP_OBJECT_METHODS(Scoring);
};




IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SCORING_H */
