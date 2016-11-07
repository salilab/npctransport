/**
 *  \file initialize_positions_RAIIs.h
 *  \brief temporary RAII objects for setting various parameters during initializiation
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */


#ifndef IMPNPCTRANSPORT_INTERNAL_INITIALIZE_POSITIONS_RAIIS_H
#define IMPNPCTRANSPORT_INTERNAL_INITIALIZE_POSITIONS_RAIISH

#include "../npctransport_config.h"
#include "TAMDChain.h"
#include "../linear_distance_pair_scores.h"
#include "../FGChain.h"
#include <IMP/scoped.h>
#include <IMP/Optimizer.h>
#include <IMP/ScoringFunction.h>
#include <IMP/atom/BrownianDynamics.h>
#include <IMP/exception.h>
#include <IMP/Pointer.h>
#include <IMP/RAII.h>
#include <IMP/raii_macros.h>
#include <IMP/core/DistancePairScore.h>
#include <IMP/core/XYZ.h>
#include <cmath>

IMPNPCTRANSPORT_BEGIN_INTERNAL_NAMESPACE


/**
An RAII class for rescaling the rest length of an FGChain
    object by some factor f (upon construction or using set()).
    Restores the original value upon destruction
*/
class FGChainScaleRestLengthRAII
: public RAII
{
private:
  PointerMember< FGChain > chain_;
  double orig_;
  bool was_set_;

 public:
  IMP_RAII(
     FGChainScaleRestLengthRAII,
      (FGChain *chain, double f),
      { //Initialize
        was_set_ = false;
      },
      { // Set
        was_set_ = true;
        chain_ = chain;
        orig_ = chain_->get_rest_length_factor();
        chain_->set_rest_length_factor(orig_ * f);
      },
      { // Reset
        if(was_set_){
          chain_->set_rest_length_factor(orig_);
        }
      },
      { // Show });
      }
   );
};

/**
An RAII class for rescaling the bond force constant k of an FGChain
    object by some factor f (upon construction or using set()).
    Restores the original value upon destruction
*/
class FGChainScaleBackboneKRAII
: public RAII
{
private:
  PointerMember< FGChain > chain_;
  double orig_;
  bool was_set_;

 public:
  IMP_RAII(
     FGChainScaleBackboneKRAII,
      (FGChain *chain, double f),
      { //Initialize
        was_set_ = false;
      },
      { // Set
        was_set_ = true;
        chain_ = chain;
        orig_ = chain_->get_backbone_k();
        chain_->set_backbone_k(orig_ * f);
      },
      { // Reset
        if(was_set_){
          chain_->set_backbone_k(orig_);
        }
      },
      { // Show });
      }
   );
};


/**
    An RAII class for rescaling the k of all TAMD restraints in a
    TAMD chain by some factor f (upon construction or using
    set()). Restores the original value upon destruction
*/
class TAMDChainScaleKRAII
: public RAII
{
private:
  PointerMember< internal::TAMDChain > chain_;
  Floats origs_;
  bool was_set_;

 public:
  IMP_RAII
    (
     TAMDChainScaleKRAII,
     ( internal::TAMDChain *chain, double f),
     { //Initialize
       was_set_ = false;
     },
     { // Set
       was_set_ = true;
       chain_ = chain;
       core::HarmonicDistancePairScores& scores =
         chain_->get_tamd_springs_byref();
       origs_.resize( scores.size() );
       for(unsigned int i = 0; i < scores.size(); i++) {
         score_functor::Harmonic& h = scores[i]->get_score_functor();
         origs_[i] = h.get_k();
         h.set_k( origs_[i] * f);
       }
     },
     { // Reset
       if(was_set_){
         core::HarmonicDistancePairScores& scores =
           chain_->get_tamd_springs_byref();
         for(unsigned int i = 0; i < scores.size(); i++) {
           score_functor::Harmonic& h = scores[i]->get_score_functor();
           h.set_k( origs_[i] );
         }
       }
     },
     { // Show });
     }
     ); // close IMP_RAII
};


/** An RAII class for temporarily changing the scoring function
    of an optimizer
*/
class OptimizerSetTemporaryScoringFunctionRAII: public RAII {

  PointerMember< IMP::ScoringFunction > orig_sf_;
  PointerMember< IMP::Optimizer > o_;
  bool was_set_;

 public:
  IMP_RAII
    ( OptimizerSetTemporaryScoringFunctionRAII,
        (IMP::Optimizer* o, ScoringFunctionAdaptor sfa),
      { //Initialize
        was_set_ = false;
      },
      { // Set
        was_set_ = true;
        o_ = o;
        orig_sf_ = o_->get_scoring_function();
        o_->set_scoring_function(sfa);
      },
      { // Reset
        if(was_set_){
          o_->set_scoring_function(orig_sf_);
        }
      },
      { // Show
      } );
};

/** An RAII class for temporarily changing the temperature
    of a BrownianDynamics object
*/
  class BDSetTemporaryTemperatureRAII : public RAII {

  double orig_temp_;
    PointerMember
    < IMP::atom::BrownianDynamics > bd_;
  bool was_set_;

 public:
  IMP_RAII
  ( BDSetTemporaryTemperatureRAII,
      (IMP::atom::BrownianDynamics* bd, double temp),
      { //Initialize
        was_set_ = false;
      },
      { // Set
        was_set_ = true;
        bd_ = bd;
        orig_temp_ = bd_->get_temperature();
        bd_->set_temperature(temp);
      },
      { // Reset
        if(was_set_){
          bd_->set_temperature(orig_temp_);
        }
      },
      { // Show
      } );
};

/** An RAII class for temporarily changing the time step
    of a BrownianDynamics object
*/
  class BDSetTemporaryTimeStepRAII : public RAII {

  double orig_time_step_;
    PointerMember
    < IMP::atom::BrownianDynamics > bd_;
  bool was_set_;

 public:
  IMP_RAII
  ( BDSetTemporaryTimeStepRAII,
      (IMP::atom::BrownianDynamics* bd, double time_step),
      { //Initialize
        was_set_ = false;
      },
      { // Set
        was_set_ = true;
        bd_ = bd;
        orig_time_step_ = bd_->get_maximum_time_step();
        bd_->set_maximum_time_step(time_step);
      },
      { // Reset
        if(was_set_){
          bd_->set_maximum_time_step(orig_time_step_);
        }
      },
      { // Show
      } );
};


/** An RAII class for temporarily setting the optimization stsate of particles
*/
class TemporarySetOptimizationStateRAII: public RAII {
  bool orig_;
  IMP::WeakPointer<IMP::Model> m_;
  IMP::ParticleIndex pi_;
  bool was_set_;

 public:
  IMP_RAII
  ( TemporarySetOptimizationStateRAII,
      (IMP::ParticleAdaptor pa, bool is_optimized),
      { //Initialize
        was_set_ = false;
      },
      { // Set
        was_set_ = true;
        m_ = pa.get_model();
        pi_ = pa.get_particle_index();
        IMP_ALWAYS_CHECK(core::XYZ::get_is_setup( m_, pi_ ),
                         "p is not XYZ - can't set coordinates opt state",
                         IMP::ValueException);
        core::XYZ xyz( m_, pi_ );
        orig_ = xyz.get_coordinates_are_optimized();
        xyz.set_coordinates_are_optimized( is_optimized );
      },
      { // Reset
        if(was_set_){
          core::XYZ(m_, pi_).set_coordinates_are_optimized( orig_ );
        }
      },
      { // Show
      } );
};


IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE

#endif /* IMPNPCTRANSPORT_INTERNAL_INITIALIZE_POSITIONS_RAIIS_H */
