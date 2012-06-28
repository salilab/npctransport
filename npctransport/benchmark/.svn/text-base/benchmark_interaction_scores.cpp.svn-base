/**
 * Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#include <IMP.h>
#include <IMP/core.h>
#include <IMP/algebra.h>
#include <IMP/npctransport.h>
#include <IMP/benchmark.h>
#include <IMP/container.h>
#include <IMP/example/optimizing.h>

using namespace IMP;
using namespace IMP::npctransport;
using namespace IMP::algebra;
using namespace IMP::benchmark;
using namespace IMP::core;
using namespace IMP::container;
const double radius=2;
#if IMP_BUILD < IMP_FAST
int number_of_particles=5;
#else
int number_of_particles=40;
#endif


ParticlesTemp create_particles(Model *m,
                               const BoundingBox3D &bb,
                               int n) {
  ParticlesTemp ret;
  for ( int i=0; i< n; ++i) {
    IMP_NEW(Particle, p, (m));
    XYZR d= XYZR::setup_particle(p);
    d.set_radius(radius);
    d.set_coordinates(get_random_vector_in(bb));
    ret.push_back(p);
    d.set_coordinates_are_optimized(true);
  }
  return ret;
}

template <class Score>
void test_one(std::string name, Score *score,
              AllBipartitePairContainer *abpc) {
  Pointer<Restraint> r= create_restraint(score, abpc);
  double scores=0;
  double time=0;
  Pointer<ScoringFunction> sf=r->create_scoring_function();
  sf->evaluate(true);
  int n=0;
  IMP_TIME({
      scores+=sf->evaluate(true);
      ++n;
    }, time);
  report("linear interaction", name, time, scores/n);
}

int main(int argc, char**argv) {
  IMP_BENCHMARK();
  BoundingBox3D bb= get_cube_d<3>(10);
  IMP_NEW(Model, m, ());
  ParticlesTemp psa= create_particles(m, bb, number_of_particles);
  ParticlesTemp psb= create_particles(m, bb, number_of_particles);
  IMP_NEW(AllBipartitePairContainer, abpc,(psa, psb));
  test_one("functor", new FunctorLinearInteractionPairScore(1,2,.5), abpc);
  test_one("custom", new LinearInteractionPairScore(1,2,.5), abpc);
  return IMP::benchmark::get_return_value();
}
