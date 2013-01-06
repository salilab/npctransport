/**
 *  \file io.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/ParticleTransportStatisticsOptimizerState.h>
#include <IMP/npctransport/Transporting.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/compatibility/set.h>
#include <IMP/compatibility/nullptr.h>
#include <IMP/base/check_macros.h>
#include <IMP/base/exception.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE


ParticleTransportStatisticsOptimizerState::ParticleTransportStatisticsOptimizerState
(Particle *p,
 Float bottom_z,
 Float top_z,
 IMP::WeakPointer<IMP::atom::Simulator> owner)
  : core::PeriodicOptimizerState("ParticleTransportStatisticsOptimizerState%1%"),
    p_(p),
    bottom_z_(bottom_z), top_z_(top_z),
    owner_(owner)
{
  this->reset();
  IMP_ALWAYS_CHECK( ! Transporting::particle_is_instance( p_ ),
                    "Particle already defined as a transporting particle,"
                    " and cannot be tracked by this object",
                    IMP::base::ValueException);
  double cur_z = core::XYZ(p_).get_coordinates()[2];
  Transporting::setup_particle(p_, false); // initial value doesn't matter

}

void
ParticleTransportStatisticsOptimizerState::reset()
{
  n_transports_up_ = 0;
  n_transports_down_ = 0;
  n_entries_bottom_ = 0;
  n_entries_top_ = 0;
  transport_time_points_in_ns_.clear();
  is_reset_ = true;
  std::cout << "ParticleTransportStatistics - RESET" << std::endl;
}

void
ParticleTransportStatisticsOptimizerState
::do_update(unsigned int) {
  const double fs_in_ns = 1.0E+6;
  Transporting p_transporting( p_ );
  double prev_z = p_transporting.get_last_tracked_z();
  double cur_z = core::XYZ(p_).get_coordinates()[2];
  p_transporting.set_last_tracked_z( cur_z ); // save for RMF, etc.
  if(is_reset_){ // ignore previous z if reset
    prev_z = cur_z;
    is_reset_ = false;
  }
  // update transport events
  if(cur_z > top_z_ && prev_z <= top_z_
     && n_entries_bottom_ > 0 && ! p_transporting.get_is_last_entry_from_top() )
    {
      n_transports_up_++;
      double time_ns = 0.0;
      if( owner_ != nullptr ) {
        time_ns = owner_->get_current_time() / fs_in_ns;
      }
      transport_time_points_in_ns_.push_back( time_ns );
      std::cout << "EXIT UP " << p_->get_name()
                << " n_transports_up = " << n_transports_up_
                << " prev_z = " << prev_z
                << std::endl;
    }
  if(cur_z < bottom_z_ && prev_z >= bottom_z_
     && n_entries_top_ > 0 && p_transporting.get_is_last_entry_from_top() ) {
    n_transports_down_++;
    double time_ns = 0.0;
    if( owner_ != nullptr ) {
      time_ns = owner_->get_current_time() / fs_in_ns;
    }
    transport_time_points_in_ns_.push_back( time_ns );
    std::cout << "EXIT DOWN " << p_->get_name()
              << " n_transports_down = " << n_transports_down_
              << " prev_z = " << prev_z
              << std::endl;
  }
  // update channel entry directionality
  if(cur_z < top_z_ && prev_z >= top_z_) {
    p_transporting.set_is_last_entry_from_top( true );
    n_entries_top_++;
  }
  if(cur_z > bottom_z_ && prev_z <= bottom_z_) {
    p_transporting.set_is_last_entry_from_top( false );
    n_entries_bottom_++;
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
