/**
 *  \file IMP/atom/BrownianDynamicsTAMDWithSlabSupport.h
 *  \brief Simple molecular dynamics optimizer.
 *
 *  Copyright 2007-2021 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPNPCTRANPORT_BROWNIAN_DYNAMICS_TAMD_WITH_SLAB_SUPPORT_H
#define IMPNPCTRANPORT_BROWNIAN_DYNAMICS_TAMD_WITH_SLAB_SUPPORT_H

#include "npctransport_config.h"
#include <IMP/atom/BrownianDynamicsTAMD.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

// for swig
class SimulationParameters;

//! Simple Brownian dynamics simulator.
/** This is an implementation of a Brownian Dynamics simulator.

    _Input particles and score_

    Each optimized particle must have x,y,z attributes
    that are optimizable. In addition, each optimized particle must be
    decorated with the Diffusion decorator. Optionally, the
    RigidBodyDiffusion decorator can be used to specify a rotational
    diffusion coefficient for core::RigidBody particles.  The
    optimizer assumes the scoring function to be energy in kcal/mol, and the xyz
    coordinates to be in angstroms and the diffusion coefficient of
    each particle be in \f$A^2/fs\f$ (or \f$Radian^2/fs\f$ for rotational
    diffusion coefficient).  Particles without optimized x,y,z
    and nonoptimized D are skipped.

    The optimizer can either automatically determine which particles
    to use from the model or be passed a SingletonContainer for the
    particles. If such a container is passed, particles added to it
    during optimization state updates are handled properly.

    _Simulation_

    At each simulation time step, each particle is translated in the
    direction of the sum of a random diffusion vector and the gradient
    of the scoring function (force field) at the particle
    coordinates. The translation is proportional to the particle
    diffusion coefficient, the time step size, and the inverse of kT.
    Note that particles masses are not considered, only their
    diffusion coefficients.

    Similarly, rigid bodies are rotated by the sum of a random torque and a
    force field torque, proportionally to the rotational diffusion
    coefficient, the time step size, and inversely proportional kT.

    If the skt (stochastic runge kutta) flag is true, the simulation is
    altered slightly to apply the SKT scheme.

    \see Diffusion
    \see RigidBodyDiffusion
  */
class IMPNPCTRANSPORTEXPORT BrownianDynamicsTAMDWithSlabSupport
: public IMP::atom::BrownianDynamicsTAMD {
 public:
  //! Create the optimizer
  /** If sc is not null, that container will be used to find particles
      to move, otherwise the model will be searched.
      @param m model associated with bd
      @param name name of bd object
      @param wave_factor for wave step function, see Simulator object,
                         if >1.001 or so, creates a wave of time steps
                         that are larger by up to wave_factor from
                         formal maximal time step

     @note wave_factor is an advanced feature - if you're not sure, just use
                       its default, see also Simulator::simulate_wave()
  */
  BrownianDynamicsTAMDWithSlabSupport(Model *m,
                       std::string name = "BrownianDynamicsTAMDWithSlabSupport%1%",
                                      double wave_factor = 1.0):
  BrownianDynamicsTAMD(m, name, wave_factor)
    {}

 protected:
  /** advances a chunk of ps from index begin to end

      @param dtfs time step in femtoseconds
      @param ikt inverse kT for current chunk step
      @param ps particle indexes to advance
      @param begin beginning index of chunk of ps
      @param end end index of chunk of ps
  */
  void do_advance_chunk(double dtfs, double ikt,
                        const ParticleIndexes &ps,
                        unsigned int begin, unsigned int end);
};

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_BROWNIAN_DYNAMICS_TAMD_WITH_SLAB_SUPPORT_H */
