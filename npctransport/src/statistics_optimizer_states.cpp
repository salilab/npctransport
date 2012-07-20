/**
 *  \file io.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/statistics_optimizer_states.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/algebra/geometric_alignment.h>
#include <IMP/atom/estimates.h>
#include <IMP/atom/Simulator.h>
#include <IMP/container/ClosePairContainer.h>
#include <IMP/container/ConsecutivePairContainer.h>
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/compatibility/set.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
BodyStatisticsOptimizerState
::BodyStatisticsOptimizerState(Particle*p):
  core::PeriodicOptimizerState("BodyStatisticsOptimizerState%1%"),
  p_(p){
}
void BodyStatisticsOptimizerState::reset() {
  positions_.clear();
  core::PeriodicOptimizerState::reset();
}
double BodyStatisticsOptimizerState::get_dt() const {
  return dynamic_cast<atom::Simulator*>(get_optimizer())
      ->get_maximum_time_step();
}

double BodyStatisticsOptimizerState::
get_correlation_time()const {
  double sum=0;
  int n=0;
  Floats angles;
  for (unsigned int i=0; i< positions_.size(); ++i) {
    double last=0;
    for (unsigned int j=i+1; j < positions_.size(); ++j) {
      algebra::Rotation3D rel= positions_[j].get_rotation()
          /positions_[i].get_rotation();
      double angle=algebra::get_axis_and_angle(rel).second;
      if (i==0) angles.push_back(angle);
      if (angle >1) {
        sum+=get_period()*get_dt()*(j-i-1 + (angle-last));
        ++n;
        break;
      }
      last=angle;
    }
  }
  /*std::cout << n << " events from " << angles
    << " with " << positions_.size() << " samples " << std::endl;*/
  if (n==0) {
    return std::numeric_limits<double>::infinity();
  }
  return sum/n;
}
double BodyStatisticsOptimizerState::get_diffusion_coefficient() const {
  algebra::Vector3Ds
    displacements(positions_.size()-1);
  for (unsigned int i=1; i< positions_.size(); ++i) {
    displacements[i-1]= positions_[i].get_translation()
        -positions_[i-1].get_translation();
  }
  return atom::get_diffusion_coefficient(displacements,
                                         get_period()*get_dt());
}


void BodyStatisticsOptimizerState
::do_update(unsigned int) {
  positions_.push_back(core::RigidBody(p_).get_reference_frame().
                          get_transformation_to());
  while (positions_.size() > 1000) {
    positions_.pop_front();
  }
}

void BodyStatisticsOptimizerState
::do_show(std::ostream&) const {
}


ChainStatisticsOptimizerState
::ChainStatisticsOptimizerState(const ParticlesTemp&p):
  core::PeriodicOptimizerState("ChainStatisticsOptimizerState%1%"),
  ps_(p) {
}
void ChainStatisticsOptimizerState::reset() {
  positions_.clear();
  core::PeriodicOptimizerState::reset();
}
double ChainStatisticsOptimizerState::get_dt() const {
  return dynamic_cast<atom::Simulator*>(get_optimizer())
      ->get_maximum_time_step();
}
double ChainStatisticsOptimizerState::
get_correlation_time()const {
  double sum=0;
  int n=0;
  for (unsigned int i=0; i< positions_.size(); ++i) {
    double last=0;
    for (unsigned int j=i+1; j < positions_.size(); ++j) {
      algebra::Transformation3D tr
          =algebra::get_transformation_aligning_first_to_second(positions_[i],
                                                                positions_[j]);
      algebra::Rotation3D rel = tr.get_rotation();
      double angle=algebra::get_axis_and_angle(rel).second;
      if (angle >1) {
        sum+=get_period()*(j-i-1 + (angle-last)) *get_dt();
        ++n;
        break;
      }
      last=angle;
    }
  }
  std::cout << n << " events" << std::endl;

  if (n==0) {
    return std::numeric_limits<double>::infinity();
  }
  return sum/n;
}

Floats ChainStatisticsOptimizerState::get_diffusion_coefficients() const {
  base::Vector<algebra::Vector3Ds >
    displacements(positions_[0].size(),
                  algebra::Vector3Ds( positions_.size()-1));
  for (unsigned int i=1; i< positions_.size(); ++i) {
    algebra::Transformation3D rel
        = algebra::get_transformation_aligning_first_to_second(positions_[i-1],
                                                               positions_[i]);
    for (unsigned int j=0; j < positions_[i].size(); ++j) {
      displacements[j][i-1]= rel.get_transformed(positions_[i-1][j])
        - positions_[i][j];
    }
  }
  Floats ret;
  for (unsigned int i=0; i < displacements.size(); ++i) {
    ret.push_back(atom::get_diffusion_coefficient(displacements[i],
                                                  get_period()*get_dt()));
  }
  return ret;
}

void ChainStatisticsOptimizerState
::do_update(unsigned int) {
  algebra::Vector3Ds vs;
  for (unsigned int i=0; i< ps_.size(); ++i) {
    vs.push_back(core::XYZ(ps_[i]).get_coordinates());
  }
  positions_.push_back(vs);
  while (positions_.size() > 1000) {
    positions_.pop_front();
  }
}

void ChainStatisticsOptimizerState
::do_show(std::ostream&) const {
}

double ChainStatisticsOptimizerState::get_diffusion_coefficient() const {
  algebra::Vector3Ds positions(positions_.size());
  for (unsigned int i=0; i< positions_.size(); ++i) {
    positions[i]= std::accumulate(positions_[i].begin(), positions_[i].end(),
                           algebra::get_zero_vector_d<3>())
        /positions_[i].size();
  }
  algebra::Vector3Ds
    displacements(positions.size()-1);
  for (unsigned int i=1; i< positions_.size(); ++i) {
    displacements[i-1]= positions[i] -positions[i-1];
  }
  return atom::get_diffusion_coefficient(displacements,
                                         get_period()*get_dt());
}

BipartitePairsStatisticsOptimizerState::BipartitePairsStatisticsOptimizerState
(Model* m,
 InteractionType interaction_type,
 const ParticlesTemp& particlesI,
 const ParticlesTemp& particlesII,
 double contact_range,
 double slack)
  :    core::PeriodicOptimizerState("BipartitePairsStatisticsOptimizerState%1%"),
       m_(m), updates_(0),
       interaction_type_(interaction_type),
       n_particles_I_( particlesI.size() ),
       n_particles_II_( particlesII.size() ) {
  reset();
  close_bipartite_pair_container_ =
    new IMP::container::CloseBipartitePairContainer
    (particlesI, particlesII,  contact_range, slack);
  close_bipartite_pair_container_->add_pair_filter
    (new container::ExclusiveConsecutivePairFilter());
  //TODO: verify all particles belong to the same model
  // for(int i = 0; i < = particlesI.size(); i++){
}

void BipartitePairsStatisticsOptimizerState::reset() {
  updates_=0;
  avg_ncontacts_=0;
  avg_pct_bound_particles_I_=0.0;
  avg_pct_bound_particles_II_=0.0;
}
namespace {
  // update the average cur_avg of n-1 observation
  // with a new observation new_val
  inline Float update_average( Float old_avg, Float new_val, Float old_n ) {
    return ( new_val + old_avg * old_n ) / (old_n + 1);
  }
}

void BipartitePairsStatisticsOptimizerState
::do_update(unsigned int) {
  // count all the pairs that are currently in contact
  // and update average
  close_bipartite_pair_container_->before_evaluate();
  unsigned int ncontacts =
    close_bipartite_pair_container_->get_number_of_particle_pairs ();
  avg_ncontacts_ =
    update_average( avg_ncontacts_, ncontacts, updates_ + 1 );

  // update the rate of particles in contact with just anybody
  // from each group
  ParticlesTemp bounds_I;
  ParticlesTemp bounds_II;
  for(unsigned int i = 0 ; i < ncontacts ; i++) {
    ParticlePair cur_pair =
      close_bipartite_pair_container_->get_particle_pair (i);
    bounds_I.push_back( cur_pair[0] );
    bounds_II.push_back( cur_pair[1] );
  }
  std::sort(bounds_I.begin(), bounds_I.end());
  std::sort(bounds_II.begin(), bounds_II.end());
  bounds_I.erase(std::unique(bounds_I.begin(), bounds_I.end()), bounds_I.end());
  bounds_II.erase(std::unique(bounds_II.begin(), bounds_II.end()),
                 bounds_II.end());
  Float pct_bound_particles_I = bounds_I.size()
    * 1.0 / n_particles_I_;
  Float pct_bound_particles_II = bounds_II.size()
    * 1.0 / n_particles_II_;
  avg_pct_bound_particles_I_ =
    update_average( avg_pct_bound_particles_I_,
                    pct_bound_particles_I,
                    updates_ + 1 );
  avg_pct_bound_particles_II_ =
    update_average( avg_pct_bound_particles_II_,
                    pct_bound_particles_II,
                    updates_ + 1);
  // std::cout <<
  //          "BipartitePairsStatisticsOptimizerState: "
  //          "Interaction type " << interaction_type_.first.get_string()
  //          << ", " << interaction_type_.second.get_string()
  //          <<" avg_pct_bound_I,II " << avg_pct_bound_particles_I_
  //          << ", " << avg_pct_bound_particles_II_
  //          << " avg_ncontacts " << avg_ncontacts_
  //           << " update # " << old_update_number + 1 << std::endl;
  ++updates_;
}

void BipartitePairsStatisticsOptimizerState
::do_show(std::ostream& o) const {
  o << "Average ncontacts " << avg_ncontacts_
    << "after " << updates_ << "updates" << std::endl;
}



IMPNPCTRANSPORT_END_NAMESPACE
