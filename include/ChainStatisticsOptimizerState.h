/**
 *  \file npctransport/ChainStatisticsOptimizerState.h
 *  \brief description
 *
 *  Copyright 2007-2018 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_CHAIN_STATISTICS_OPTIMIZER_STATE_H
#define IMPNPCTRANSPORT_CHAIN_STATISTICS_OPTIMIZER_STATE_H

#include "npctransport_config.h"
#include <IMP/Particle.h>
#include <IMP/OptimizerState.h>
//#include <IMP/optimizer_state_macros.h>
#include <IMP/core/PeriodicOptimizerState.h>
#include <IMP/npctransport/typedefs.h>
#include <deque>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/** Compute various statistics of a chain.*/
class IMPNPCTRANSPORTEXPORT ChainStatisticsOptimizerState
    : public core::PeriodicOptimizerState {
 private:
  typedef core::PeriodicOptimizerState P;

  // particles in the chain:
  ParticlesTemp ps_;

  // time series of the positions of particles in the chain:
  std::deque<algebra::Vector3Ds> positions_;

  // mean radius-of-gyration and its square since last reset
  double mean_rgyr_;
  double mean_rgyr2_;

  // mean end-to-end length and its square since last reset
  double mean_end_to_end_;
  double mean_end_to_end2_;

  // mean bond distance and its square since last reset
  double mean_bond_distance_;
  double mean_bond_distance2_;

  // number of samples over which means are computed
  int n_;

  double get_dt() const;

 public:
  /**
     @param ps the particles being wrapped
     @param periodicity frame interval for statistics, equiv. to set_period(1)
   */
  ChainStatisticsOptimizerState(const ParticlesTemp &ps,
                                unsigned int periodicity = 1);

  double get_correlation_time() const;

  //! returns a vector of diffusion coefficients
  //! for each particle in the chain, computed in the local
  //! reference frame of the chain (by locally aligning
  //! the chain)
  Floats get_local_diffusion_coefficients() const;

  //! get an estimate of the diffusion coefficient
  //! of the entire chain in A^2/fs units
  double get_diffusion_coefficient() const;

  //! returns the mean Rgyr of this chain
  double get_mean_radius_of_gyration() const{
    return mean_rgyr_;
  }

  //! returns the mean Rgyr^2, which can be used to compute
  //! std-dev (but recorded separately so it could be averaged with
  //! other chains)
  double get_mean_square_radius_of_gyration() const{
    return mean_rgyr2_;
  }

    //! returns the mean end-to-end distance of this chain
  double get_mean_end_to_end_distance() const{
    return mean_end_to_end_;
  }

  //! returns the mean square end-to-end distance, which can be used to compute
  //! std-dev (but recorded separately so it could be averaged with
  //! other chains)
  double get_mean_square_end_to_end_distance() const{
    return mean_end_to_end2_;
  }

    //! returns the mean bond distance of this chain
  double get_mean_bond_distance() const{
    return mean_bond_distance_;
  }

  //! returns the mean square bond distance, which can be used to compute
  //! std-dev (but recorded separately so it could be averaged with
  //! other chains)
  double get_mean_square_bond_distance() const{
    return mean_bond_distance2_;
  }

  /**
     Resets all the statistics about that chain
  */
  void reset();
  virtual void do_update(unsigned int call_num) IMP_OVERRIDE;
  IMP_OBJECT_METHODS(ChainStatisticsOptimizerState);
};
IMP_OBJECTS(ChainStatisticsOptimizerState, ChainStatisticsOptimizerStates);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_CHAIN_STATISTICS_OPTIMIZER_STATE_H */
