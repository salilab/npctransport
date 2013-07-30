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

  /* NOTE: the following ar mutable cause they do not affect the
     external logical state of the scoring functions
     (as they're all updated before being revealed externally
     based on the non-mutable parameters and owning sd
  */

  // true if the list of particles in the owner_sd has changed
  // and in the middle of changing it
  //  mutable bool is_updating_particles_;

  // see get_close_diffusers_container()
  mutable base::Pointer
    <IMP::PairContainer> close_diffusers_container_;

  // generates hash values ('predicates') for ordered types pairs
  // (e.g., pairs of ParticleTypes)
  mutable base::Pointer
    <IMP::core::OrderedTypePairPredicate> otpp_;

  mutable base::Pointer
    <IMP::core::RestraintsScoringFunction> scoring_function_;

  // contains all restraints between pairs of particle types
  mutable base::Pointer
    <container::PredicatePairsRestraint> predr_;

  // scores between particle, mapped by their interaction id,
  // where interaction id is from OrderedTypePairPredicate otpp_
  // applied on each interacting particle type
  mutable IMP::base::map< int, base::Pointer< IMP::PairScore > >
    sites_pair_scores_;

  mutable base::Pointer
    <IMP::Restraint> box_restraint_;

  mutable base::Pointer
    <IMP::Restraint> slab_restraint_;

  mutable Restraints chain_restraints_;

  mutable LinearWellPairScores chain_scores_;


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
     return a scoring function for the NPC simulation
   */
  IMP::ScoringFunction* get_scoring_function() const;

  /** update the scoring function that the list of particles,
      to be retrieved from owning simulation data, has changed
   */
  //  void update_particles() const{
  //    is_updating_particles_ = true;
  //  }


  /** return the SimulationData object that owns this ScoringFunction */
  SimulationData* get_sd() {
    return owner_sd_;
  }


  /**
     a pair container that returns all pairs of particles (a,b) such that:
     1) a is an optimizable diffuser particle and b is any diffuser
        particle (static or not).
     2) a and b are close (sphere surfaces within range get_range())
     3) a and b do not appear consecutively within the model (e.g., fg repeats)

     @note TODO: right now will not return any consecutive particles - that's
                 erroneous though may be negligible
  */
  IMP::PairContainer *get_close_diffusers_container() const;

  /**
     Returns the container for restraints over pairs of particles. Different
     scores
     are used for particles of different (ordered) particle types.
     When called for the first time, returns a new PredicatePairsRestraints
     over all diffusing particles and sets a default linear repulsion restraint
     between all pairs returned by get_close_diffusers_container()
  */
  container::PredicatePairsRestraint *get_predr() const;

  // returns true if a bounding box restraint has been defined */
  bool get_has_bounding_box() const
  { return box_restraint_ != nullptr && box_is_on_; }

  /** returns true if a slab restraint has been defined */
  bool get_has_slab() const
  { return slab_restraint_ != nullptr && slab_is_on_; }

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
  void add_interaction(
      const ::npctransport_proto::Assignment_InteractionAssignment &idata);
#endif

  double get_backbone_k() const { return backbone_k_; }

  double get_interaction_k() const { return interaction_k_; }

  /** create an fg chain restraint on consecutive chain particles

      @note this method assumes that all such chains will be disjoint
      and so it is later possible to use
      container::ExclusiveConsecutivePairFilter to filter out all
      pairs of particles connected by such chain restraints.

      @param particles the particles in the chain in consecutive
                       order
      @param rest_length the rest length of consecutive chain beads
      @param name the name of the chain (to be used for naming the restraint

      @return pointer to the newly created restraint
  */
  Restraint* add_chain_restraint(const ParticlesTemp &particles,
                                 double rest_length,
                                 std::string name);

  /** returns all restraints on fg nup chains */
  Restraints get_chain_restraints() const
  { return chain_restraints_; }

  /**
     returns a reference to the collection of score functions for FG backbones
     (can be used to e.g. scale them up or down during optimization)
   */
  LinearWellPairScores get_chain_scores() const { return chain_scores_; }


  Restraint *get_box_restraint() const {
    if ((/*is_updating_particles_ ||*/ !box_restraint_) && box_is_on_){
      const_cast<Scoring *>(this)->
        create_bounding_box_restraint_on_diffusers();
    }
    return box_restraint_;
  }

  Restraint *get_slab_restraint() const {
    if ((/*is_updating_particles_ ||*/ !slab_restraint_) && slab_is_on_){
      const_cast<Scoring *>(this)->
        create_slab_restraint_on_diffusers();
    }
    return slab_restraint_;
  }

  double get_slab_thickness() const { return slab_thickness_; }

  /**
     Creates bounding box restraint based on the box_size_
     class variable, and apply it to all diffusers returned
     by get_sd()->get_diffusers()
  */
  void create_bounding_box_restraint_on_diffusers();

  /**
     Creates slab bounding volume restraint, based on the slab_thickness_
     and tunnel_radius_ class variables, and apply it to all diffusers
     returned by get_sd()->get_diffusers()
   */
  void create_slab_restraint_on_diffusers();

  double get_range() const { return range_; }

  /**
      sets the multiplicative scaling of the interaction range for
      pair interactions involving a particle of this type
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
      pair interactions  involving a particle of this type
  */
  void set_interaction_k_factor
    ( IMP::core::ParticleType type, double value )
    {
      IMP_ALWAYS_CHECK(value>0, "interaction_k_factor must be positive",
                       IMP::base::ValueException);
      interaction_k_factors_[type] = value;
    }


  /**
     returns all diffusing particles that currently exist in
     this simulation data object get_sd() (or an empty list if none exists)
   */
  container::ListSingletonContainer const* get_diffusers() const;


 public:
  IMP_OBJECT_METHODS(Scoring);
};


IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SCORING_H */
