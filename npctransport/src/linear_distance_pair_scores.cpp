/**
 *  \file DistancePairScore.cpp
 *  \brief A Score on the distance between a pair of particles.
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
 */

#include <IMP/npctransport/linear_distance_pair_scores.h>
IMPNPCTRANSPORT_BEGIN_NAMESPACE



LinearSoftSpherePairScore
::LinearSoftSpherePairScore(double k, std::string name):
    PairScore(name),
    k_(k) {
}
void LinearSoftSpherePairScore::do_show(std::ostream &) const
{
}
ParticlesTemp LinearSoftSpherePairScore::get_input_particles(Particle*p) const {
  return ParticlesTemp(1, p);
}
ContainersTemp
LinearSoftSpherePairScore::get_input_containers(Particle*) const {
  return ContainersTemp();
}



LinearInteractionPairScore
::LinearInteractionPairScore(double krep, double attr_range, double kattr,
                             std::string name):
    PairScore(name),
    attr_range_(attr_range), k_rep_(krep), k_attr_(kattr)  {
}

void LinearInteractionPairScore::do_show(std::ostream & o) const
{
  o << "k_repulsive = " << k_rep_ << " ; k_attraction = " << k_attr_
    << " ; attraction range = " << attr_range_ << std::endl;
}
ParticlesTemp LinearInteractionPairScore
::get_input_particles(Particle*p) const {
  return ParticlesTemp(1, p);
}
ContainersTemp
LinearInteractionPairScore::get_input_containers(Particle*) const {
  return ContainersTemp();
}




LinearWellPairScore
::LinearWellPairScore(double x0, double k,
                      std::string name):
    PairScore(name),
    x0_(x0), k_(k) {
}

void LinearWellPairScore::do_show(std::ostream & ) const
{
}
ParticlesTemp LinearWellPairScore
::get_input_particles(Particle*p) const {
  return ParticlesTemp(1, p);
}
ContainersTemp
LinearWellPairScore::get_input_containers(Particle*) const {
  return ContainersTemp();
}



IMPNPCTRANSPORT_END_NAMESPACE
