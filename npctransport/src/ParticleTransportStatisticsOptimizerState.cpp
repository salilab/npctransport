/**
 *  \file io.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/ParticleTransportStatisticsOptimizerState.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/compatibility/set.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE


ParticleTransportStatisticsOptimizerState::ParticleTransportStatisticsOptimizerState
(Particle *p, Float bottom_z, Float top_z) :
  core::PeriodicOptimizerState("ParticleTransportStatisticsOptimizerState%1%"),
  p_(p),
  bottom_z_(bottom_z), top_z_(top_z)
{
  this->reset();
  cur_z_ = core::RigidBody(p_).get_coordinates()[2];
}

void
ParticleTransportStatisticsOptimizerState::reset()
{
  n_transports_up_ = 0;
  n_transports_down_ = 0;
  n_entries_bottom_ = 0;
  n_entries_top_ = 0;
  is_reset_ = true;
  std::cout << "ParticleTransportStatistics - RESET" << std::endl;
}

void
ParticleTransportStatisticsOptimizerState
::do_update(unsigned int) {
  prev_z_ = cur_z_;
  cur_z_ = core::RigidBody(p_).get_coordinates()[2];
  if(is_reset_){ // ignore previous z if reset
    prev_z_ = cur_z_;
    is_reset_ = false;
  }
  // update transport events
  if(cur_z_ > top_z_ && prev_z_ <= top_z_
     && n_entries_bottom_ > 0 && !is_last_entry_from_top_) {
    n_transports_up_++;
    std::cout << "Particle EXIT UP " << std::endl
              << *p_ << std::endl
              << " n_transports_up = " << n_transports_up_
              << " prev_z = " << prev_z_
              << std::endl;
  }
  if(cur_z_ < bottom_z_ && prev_z_ >= bottom_z_
     && n_entries_top_ > 0 && is_last_entry_from_top_) {
    n_transports_down_++;
    std::cout << "Particle EXIT DOWN" << std::endl
              << *p_ << std::endl
              << " n_transports_down = " << n_transports_down_
              << " prev_z = " << prev_z_
              << std::endl;
  }
  // update channel entry directionality
  if(cur_z_ < top_z_ && prev_z_ >= top_z_) {
    is_last_entry_from_top_ = true;
    n_entries_top_++;
    std::cout << "Particle ENTRY TOP" << std::endl
              << *p_ << std::endl
              << " n_entries_top = " << n_entries_top_
              << " prev_z = " << prev_z_
              << std::endl;
  }
  if(cur_z_ > bottom_z_ && prev_z_ <= bottom_z_) {
    is_last_entry_from_top_ = false;
    n_entries_bottom_++;
    std::cout << "Particle ENTRY BOTTOM " << std::endl
              << *p_ << std::endl
              << " n_entries_bottom = " << n_entries_bottom_
              << " prev_z = " << prev_z_
              << std::endl;
  }
}

void ParticleTransportStatisticsOptimizerState
::do_show(std::ostream& o) const {
  o << "Particle "<< *p_
    << " Number of transports up: " << n_transports_up_
    << " Number of transports down: " << n_transports_down_
    << std::endl;
}
IMPNPCTRANSPORT_END_NAMESPACE
