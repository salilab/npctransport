/**
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#include <IMP/npctransport/SitesPairScore.h>
#include <IMP/npctransport/util.h>

#include <IMP/algebra/Transformation3D.h>
#include <IMP/algebra/Sphere3D.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/algebra/vector_generators.h>
#include <IMP/benchmark/benchmark_macros.h>
#include <IMP/benchmark/utility.h>
#include <IMP/container/AllBipartitePairContainer.h>
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/container/generic.h>
#include <IMP/container/ClosePairContainer.h>
#include <IMP/core/XYZR.h>
#include <IMP/core/MonteCarlo.h>
#include <IMP/core/SerialMover.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/core/ConjugateGradients.h>
#include <IMP/core/BallMover.h>
#include <IMP/core/SphereDistancePairScore.h>
#include <IMP/log_macros.h>
#include <IMP/flags.h>
#include <IMP/Model.h>
#include <IMP/Particle.h>
#include <IMP/Pointer.h>
#include <IMP/Restraint.h>
#include <IMP/scoped.h>
#include <IMP/PairPredicate.h>
#include <string>
#include <boost/ptr_container/ptr_vector.hpp>

using namespace IMP;
using namespace IMP::npctransport;
using namespace IMP::algebra;
using namespace IMP::benchmark;
using namespace IMP::core;
using namespace IMP::container;

namespace {
const double radius = 2;
#if IMP_BUILD < IMP_FAST
int number_of_particles = 5;
  //int step_size = 8;
#else
int number_of_particles = 40;
  //int step_size = 5;
#endif

// void debug_print_location(std::string context, bool reset = false) {
//   static int ncall = 0;
//   if (reset) ncall = 0;
//   std::cout << "Here " << ncall++ << " in context " << context << std::endl;
// }

/**
   sites 1 1 0.2,                8.85e-04,           -1.3e+04
   sites 1 1 0.4,                1.23e-03,           -2.0e+04
   sites 1 1 0.8,                9.03e-04,           -5.4e+04
   sites 6 1 0.2,                2.35e-03,           -3.3e+04
   sites 6 1 0.4,                2.40e-03,           -6.1e+04
   sites 6 1 0.8,                2.35e-03,           -1.3e+05
   sites 6 6 0.2,                1.98e-02,           -2.3e+04
   sites 6 6 0.4,                2.98e-02,           -3.1e+04
   sites 6 6 0.8,                1.14e-02,           -1.6e+05
   sites 11 1 0.2,               3.84e-03,           -3.6e+04
   sites 11 1 0.4,               3.83e-03,           -7.3e+04
   sites 11 1 0.8,               4.87e-03,           -1.2e+05
   sites 11 6 0.2,               5.31e-02,           -1.6e+04
   sites 11 6 0.4,               2.66e-02,           -6.2e+04
   sites 11 6 0.8,               1.93e-02,           -1.7e+05
   sites 11 11 0.2,              3.45e-02,           -4.4e+04
   sites 11 11 0.4,              3.51e-02,           -8.7e+04
   sites 11 11 0.8,              9.52e-02,           -6.6e+04
   sites 16 1 0.2,               1.02e-02,           -2.0e+04
   sites 16 1 0.4,               5.24e-03,           -7.7e+04
   sites 16 1 0.8,               5.62e-03,           -1.4e+05
   sites 16 6 0.2,               2.81e-02,           -4.3e+04
   sites 16 6 0.4,               2.93e-02,           -8.3e+04
   sites 16 6 0.8,               7.67e-02,           -6.3e+04
   sites 16 11 0.2,              9.07e-02,           -2.4e+04
   sites 16 11 0.4,              5.04e-02,           -8.7e+04
   sites 16 11 0.8,              5.02e-02,           -1.8e+05
   sites 16 16 0.2,              7.47e-02,           -4.4e+04
   sites 16 16 0.4,              1.42e-01,           -4.6e+04
   sites 16 16 0.8,              1.96e-01,           -6.7e+04
*/

/**
   create n optimizable rigid body particles for model m
   with random coordinates and orientation within bounding box bb
*/
ParticleIndexes create_particles(Model *m, const BoundingBox3D &bb, int n) {
  ParticleIndexes ret;
  for (int i = 0; i < n; ++i) {
    IMP_NEW(Particle, p, (m));
    ParticleIndex pi=p->get_index();
    XYZR d = XYZR::setup_particle(m,pi);
    d.set_radius(radius);
    d.set_coordinates(get_random_vector_in(bb));
    ret.push_back(pi);
    d.set_coordinates_are_optimized(true);
    Transformation3D tr(get_random_rotation_3d(), get_zero_vector_d<3>());
    ReferenceFrame3D rf(tr);
    RigidBody::setup_particle(m, pi, rf);
  }
  return ret;
}

  core::MonteCarloMover *create_serial_mover(IMP::Model* m, const ParticleIndexes &pis) {
  core::MonteCarloMovers movers;
  for (unsigned int i = 0; i < pis.size(); ++i) {
    double scale = core::XYZR(m, pis[i]).get_radius();
    movers.push_back(new core::BallMover(m, pis[i], scale * 2));
  }
  IMP_NEW(core::SerialMover, sm, (get_as<core::MonteCarloMoversTemp>(movers)));
  return sm.release();
}

/** Take a set of core::XYZR particles ps, and relax them relative to a set of
    restraints rs. Excluded volume is handle separately, so don't include it
    in 'rs', but only in 'excluded'.

    @param ps list of particles to be optimized
    @param rs list of restraints on particles (don't include excluded volume
   here)
    @param excluded XXXX (TODO: what is it exactly?) XXXX
    @param opt_states optimizer states to be added to the various optimizers
    (e.g. conjugate gradient, Monte-Carlo) during optimization

    \include optimize_balls.cpp
*/
  void optimize_balls(IMP::Model* m,
                      const ParticleIndexes &pis,
                      const RestraintsTemp &rs = RestraintsTemp(),
                      const PairPredicates &excluded = PairPredicates(),
                      const OptimizerStates &opt_states = OptimizerStates(),
                      LogLevel ll = DEFAULT) {
  // make sure that errors and log messages are marked as coming from this
  // function
  IMP_FUNCTION_LOG;
  SetLogState sls(ll);
  IMP_ALWAYS_CHECK(!pis.empty(), "No Particles passed.", ValueException);
  // double scale = core::XYZR(ps[0]).get_radius();
  Particles ps;
  for(unsigned int i=0; i<pis.size(); i++){
    ps.push_back(m->get_particle(pis[i]));
  }

  IMP_NEW(core::SoftSpherePairScore, ssps, (10));
  IMP_NEW(core::ConjugateGradients, cg, (m));
  cg->set_optimizer_states(opt_states);
  {
    // set up restraints for cg
    IMP_NEW(container::ListSingletonContainer, lsc, (m, pis));
    IMP_NEW(container::ClosePairContainer, cpc,
            (lsc, 0, core::XYZR(m, pis[0]).get_radius()));
    cpc->add_pair_filters(excluded);
    Pointer<Restraint> r =
        container::create_restraint(ssps.get(), cpc.get());
    cg->set_scoring_function(rs + RestraintsTemp(1, r.get()));
    cg->set_optimizer_states(opt_states);
  }
  IMP_NEW(core::MonteCarlo, mc, (m));
  mc->set_optimizer_states(opt_states);
  IMP_NEW(core::IncrementalScoringFunction, isf, (m, pis, rs));
  {
    // set up MC
    mc->add_mover(create_serial_mover(m, pis));
    // we are special casing the nbl term for montecarlo, but using all for CG
    mc->set_incremental_scoring_function(isf);
    // use special incremental support for the non-bonded part
    isf->add_close_pair_score(ssps, 0, ps, excluded);
    // make pointer vector
  }

  IMP_LOG_PROGRESS("Performing initial optimization" << std::endl);
  {
    boost::ptr_vector<ScopedSetFloatAttribute> attrs;
    for (unsigned int j = 0; j < attrs.size(); ++j) {
      attrs.push_back(
          new ScopedSetFloatAttribute(ps[j], core::XYZR::get_radius_key(), 0));
    }
    cg->optimize(1000);
  }
  // shrink each of the particles, relax the configuration, repeat
  for (int i = 0; i < 11; ++i) {
    boost::ptr_vector<ScopedSetFloatAttribute> attrs;
    double factor = .1 * i;
    IMP_LOG_PROGRESS("Optimizing with radii at " << factor << " of full"
                                                 << std::endl);
    for (unsigned int j = 0; j < pis.size(); ++j) {
      attrs.push_back(
          new ScopedSetFloatAttribute(ps[j], core::XYZR::get_radius_key(),
                                      core::XYZR(m,pis[j]).get_radius() * factor));
    }
    // changed all radii
    isf->set_moved_particles(isf->get_movable_indexes());
    for (int j = 0; j < 5; ++j) {
      mc->set_kt(100.0 / (3 * j + 1));
      mc->optimize(pis.size() * (j + 1) * 100);
      double e = cg->optimize(10);
      IMP_LOG_PROGRESS("Energy is " << e << std::endl);
      if (e < .000001) break;
    }
  }
}

/**
   create particles with sites on them, and with appropriate restraints,
   and benchmark whether they work well
 */
template <unsigned int NA, unsigned int NB, bool WHICH>
void test_one(double range) {
  BoundingBox3D bb = get_cube_d<3>(50);
  IMP_NEW(Model, m, ());
  ParticleIndexes psa = create_particles(m, bb, number_of_particles);
  ParticleIndexes psb = create_particles(m, bb, number_of_particles);
  Sphere3D s(get_zero_vector_d<3>(), radius);
  Vector3Ds sas = get_uniform_surface_cover(s, NA);
  Vector3Ds sbs = get_uniform_surface_cover(s, NB);
  optimize_balls(m, psa + psb);
  //  typedef TemplateSitesPairScore<NA, NB, WHICH> TSPS;
  //IMP_NEW(TSPS, tsps, (range, 1, 0, 0, 1, sas, sbs));
  IMP_NEW(SitesPairScore, sps, ( range, 1, // r, l
                                 0.0, 0.0, // no-skew version
                                 0, 0, 1, // non-specific r, k_attr, k_rep
                                 get_spheres_from_vectors(sas, 0.0),
                                 get_spheres_from_vectors(sbs, 0.0) )
          );
  IMP_NEW(ListSingletonContainer, lsca, (m, psa));
  IMP_NEW(ListSingletonContainer, lscb, (m, psb));
  IMP_NEW(AllBipartitePairContainer, abpc, (lsca, lscb));
  if (WHICH) {
    Pointer<Restraint> r = create_restraint(sps.get(), abpc.get());
    double scores = 0;
    double time = 0;
    IMP_TIME({
      scores += r->evaluate(true);
    },
             time);
    std::ostringstream oss;
    oss << "sites " << NA << " " << NB << " " << range;
    report(oss.str(), time, scores);
  }
  // {
  //   Pointer<Restraint> r = create_restraint(tsps.get(), abpc.get());
  //   double scores = 0;
  //   double time = 0;
  //   IMP_TIME({
  //     scores += r->evaluate(true);
  //   },
  //            time);
  //   std::ostringstream oss;
  //   oss << "template sites " << NA << " " << NB << " " << WHICH << " " << range;
  //   report(oss.str(), time, scores);
  // }
}

template <unsigned int NA, unsigned int NB, bool WHICH>
void test_ranges() {
  for (double range = .1 *radius; range < .5 *radius; range *= 2) {
    test_one<NA, NB, WHICH>(range);
  }
}
}

int main(int, char **) {
  if(IMP::get_check_level() >= IMP::USAGE){
    std::cout << "Skipping test because check level is greater or equal USAGE"
              << std::endl;
    return 0;
  }
  test_ranges<1, 1, false>();
  test_ranges<1, 1, true>();
  test_ranges<4, 1, false>();
  test_ranges<4, 1, true>();
  test_ranges<4, 4, false>();
  test_ranges<4, 4, true>();
  test_ranges<16, 1, false>();
  test_ranges<16, 1, true>();
  test_ranges<16, 4, false>();
  test_ranges<16, 4, true>();
  test_ranges<16, 8, false>();
  test_ranges<16, 8, true>();
  test_ranges<16, 16, false>();
  test_ranges<16, 16, true>();
  return IMP::benchmark::get_return_value();
}
