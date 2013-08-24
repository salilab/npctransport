/**
 *  \file BipartitePairsStatisticsOptimizerState.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/BipartitePairsStatisticsOptimizerState.h>
#include <IMP/container/ClosePairContainer.h>
#include <IMP/pair_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

namespace {
// update the average cur_avg of n-1 observation
// with a new observation new_val
inline Float update_average(Float old_avg, Float new_val, Float old_n) {
  return (new_val + old_avg * old_n) / (old_n + 1);
}
}

BipartitePairsStatisticsOptimizerState::BipartitePairsStatisticsOptimizerState(
    Model* m, InteractionType interaction_type, const ParticlesTemp& particlesI,
    const ParticlesTemp& particlesII, double contact_range, double slack)
  : P(m, "BipartitePairsStatisticsOptimizerState%1%"),
      updates_(0),
      interaction_type_(interaction_type),
      n_particles_I_(particlesI.size()),
      n_particles_II_(particlesII.size()) {
  reset();
  close_bipartite_pair_container_ =
      new IMP::container::CloseBipartitePairContainer(particlesI, particlesII,
                                                     contact_range, slack);
  range_ = contact_range;
  // TODO: do we want to add consecutive pair container for fg chains?
}

void BipartitePairsStatisticsOptimizerState::reset() {
  updates_ = 0;
  avg_ncontacts_ = 0;
  avg_pct_bound_particles_I_ = 0.0;
  avg_pct_bound_particles_II_ = 0.0;
}

void BipartitePairsStatisticsOptimizerState::do_update(unsigned int) {
  // count all the pairs that are currently in contact
  // and update average
  close_bipartite_pair_container_->before_evaluate();

  // update the rate of particles in contact with just anybody
  // from each group
  ParticleIndexes bounds_I;
  ParticleIndexes bounds_II;
  unsigned int ncontacts = 0;
  IMP_CONTAINER_FOREACH(IMP::container::CloseBipartitePairContainer,
                        close_bipartite_pair_container_,
                        {
                          // _1 = ParticleIndexPair
                          Model * m = get_model();
                          // filter out sphere pairs with distances in
                          // slack region
                          core::XYZR s0(m, _1[0]);
                          core::XYZR s1(m, _1[1]);
                          if(core::get_distance(s0, s1) < range_)
                            {
                              bounds_I.push_back(_1[0]);
                              bounds_II.push_back(_1[1]);
                              ncontacts++;
                            }
                        });
  avg_ncontacts_ = update_average(avg_ncontacts_, ncontacts, updates_ + 1);
  std::sort(bounds_I.begin(), bounds_I.end());
  std::sort(bounds_II.begin(), bounds_II.end());
  bounds_I.erase(std::unique(bounds_I.begin(), bounds_I.end()), bounds_I.end());
  bounds_II.erase(std::unique(bounds_II.begin(), bounds_II.end()),
                  bounds_II.end());
  pct_bound_particles_I_ = bounds_I.size() * 1.0 / n_particles_I_;
  pct_bound_particles_II_ = bounds_II.size() * 1.0 / n_particles_II_;
  avg_pct_bound_particles_I_ = update_average(
      avg_pct_bound_particles_I_, pct_bound_particles_I_, updates_ + 1);
  avg_pct_bound_particles_II_ = update_average(
      avg_pct_bound_particles_II_, pct_bound_particles_II_, updates_ + 1);
  ++updates_;
}

IMPNPCTRANSPORT_END_NAMESPACE
