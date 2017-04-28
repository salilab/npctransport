/**
 *  \file HarmonicSpringSingletonScore.h
 *  \brief Harmonic scores on the spring between a pair of particles.
 *         The spring resting length is dynamic
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_HARMONIC_SPRING_SINGLETON_SCORE_H
#define IMPNPCTRANSPORT_HARMONIC_SPRING_SINGLETON_SCORE_H

#include "npctransport_config.h"
#include "RelaxingSpring.h"
#include <IMP/compiler_macros.h>
#include <IMP/SingletonScore.h>
#include <IMP/pair_macros.h>
#include <IMP/core/XYZR.h>
#include <IMP/algebra/utility.h>

#include <boost/array.hpp>

IMPNPCTRANSPORT_BEGIN_NAMESPACE


/**
   A harmonic score between two bonded spheres
   connected by a spring whose resting length is dynamic
*/
class IMPNPCTRANSPORTEXPORT HarmonicSpringSingletonScore
: public SingletonScore
{
 private:
  double k1_;
  double k2_;

 public:
  /**
     a harmonic  well pair potential that keeps two particles connected
     by a spring particle whose resting length itself is dynamic

     @param k1 the force constant applied by the string on tethered
              particles in units of kcal/mol/A^2 - attractive if spring
              is longer than rest length, repulsive if shorter than rest
              length. The spring applies a force of k1*dX on each particle,
	      and a counter force of 2*k1*dX acts on the spring, where
	      dX=|distance(p1,p2)-(rest_length)|
     @param k2 the force constant for relaxation of the spring
              rest length to its equilibrium value in kcal/mol/A^2.
	      The force equals k2*dX where dX=|(rest_length)-(equilibrium length)|
     @param name the name of the score
   */
  HarmonicSpringSingletonScore(double k1, double k2,
                      std::string name = "HarmonicIDSingletonScore%1%");

  void set_k1(double k1)
  { k1_ = k1; }

  double get_k1() const
  { return k1_; }

  void set_k2(double k2)
  { k2_ = k2; }

  double get_k2() const
  { return k2_; }

  double evaluate_index(Model *m, const ParticleIndex pi,
                        DerivativeAccumulator *da) const IMP_OVERRIDE;

  ModelObjectsTemp do_get_inputs(Model *m, const ParticleIndexes &pis) const;

  IMP_SINGLETON_SCORE_METHODS(HarmonicSpringSingletonScore);

  IMP_OBJECT_METHODS(HarmonicSpringSingletonScore);
  ;
};

#ifndef IMP_DOXYGEN

inline double
HarmonicSpringSingletonScore
::evaluate_index
( Model *m, const ParticleIndex pi, DerivativeAccumulator *da) const
{
  IMP_OBJECT_LOG;

  IMP_USAGE_CHECK(RelaxingSpring::get_is_setup(m,pi),
		  "particle 0 is expected to be string in HarmonicSpringSingletonScore");
  RelaxingSpring s(m, pi);
  ParticleIndex pi0(s.get_bonded_particle_index_0());
  ParticleIndex pi1(s.get_bonded_particle_index_1());
  algebra::Sphere3D s0 = m->get_sphere(pi0);
  algebra::Sphere3D s1 = m->get_sphere(pi1);
  double rest_delta_length= s.get_rest_length();
  algebra::Vector3D delta_1_to_0 = s0.get_center() - s1.get_center();
  double delta_length_2 = delta_1_to_0.get_squared_magnitude();
  double delta_length = std::sqrt(delta_length_2);
  double dDelta = delta_length - rest_delta_length; // positive if spring is extended, negative if compressed
  double scoreDelta = 2 * 0.5 * k1_ * dDelta * dDelta; // x2 because each particle applies force on the spring independently in this model of a relaxing spring
  double dEq= rest_delta_length - s.get_equilibrium_rest_length(); // positive is rest length is stretched relative to equilibrium rest length
  double scoreEq = 0.5 * k2_ * dEq * dEq;
  double score= scoreDelta + scoreEq;
  IMP_LOG(TERSE, "dDelta: " << dDelta << " scoreDelta: " << scoreDelta
	  << " dEq: " << dEq << " scoreEq: " << scoreEq
	  << " total: " << score);

  // Compute derivatives:
  static const double MIN_DISTANCE = .00001;
  if (IMP_LIKELY( da && delta_length > MIN_DISTANCE )) { // Profiling note on use of likely(): in BD simulations, the simulation bottleneck is when da is true, and the spring is likely out of equilibrium
    double fParticles= k1_*dDelta; // force pulling particles closer together for a positive force (or apart if negative)
    double fSpring= k2_*dEq - 2*fParticles; // force pulling rest length down for a positive force (or up for negtive); the -2*f1 is the coutnerforce exerted on the spring by the two tethered particles
    s.add_to_rest_length_derivative(fSpring, *da);
    algebra::Vector3D deriv0( delta_1_to_0 * (fParticles / delta_length) );
    m->add_to_coordinate_derivatives(pi0, deriv0, *da);
    m->add_to_coordinate_derivatives(pi1, -deriv0, *da);
    IMP_LOG(TERSE, "\nderiv on pi0: " << deriv0);
    IMP_LOG(TERSE, "\nderiv on spring: " << -fSpring);
  }
  IMP_LOG(TERSE, std::endl);

  return score;
}

#endif

IMP_OBJECTS(HarmonicSpringSingletonScore, HarmonicSpringSingletonScores);


IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_HARMONIC_SPRING_SINGLETON_SCORE_H */
