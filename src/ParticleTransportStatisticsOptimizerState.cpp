/**
 *  \file io.cpp
 *  \brief description.
 *
 *  Copyright 2007-2020 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/ParticleTransportStatisticsOptimizerState.h>
#include <IMP/npctransport/Transporting.h>
#include <IMP/npctransport/Statistics.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/check_macros.h>
#include <IMP/exception.h>
#include <IMP/log.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

ParticleTransportStatisticsOptimizerState::
    ParticleTransportStatisticsOptimizerState(
        Particle *p, Float bottom_z, Float top_z,
        WeakPointer<IMP::npctransport::Statistics> statistics_manager,
        WeakPointer<IMP::atom::Simulator> owner)
      : P(p->get_model(), "ParticleTransportStatisticsOptimizerState%1%"),
        p_(p),
        bottom_z_(bottom_z),
        top_z_(top_z),
        statistics_manager_(statistics_manager),
        owner_(owner) {
  IMP_ALWAYS_CHECK(!Transporting::get_is_setup(p_),
                   "Particle already defined as a transporting particle,"
                   " and cannot be tracked by this object",
                   IMP::ValueException);
  Transporting::setup_particle(p_, false);  // initial value doesn't matter
  this->reset();
}

void ParticleTransportStatisticsOptimizerState::reset() {
  P::reset();
  n_transports_up_ = 0;
  n_transports_down_ = 0;
  transport_time_points_in_ns_.clear();
  Transporting(p_).set_n_entries_bottom(0);
  Transporting(p_).set_n_entries_top(0);
  is_reset_ = true;
  IMP_LOG(PROGRESS, "ParticleTransportStatistics - RESET" << std::endl);
}

void ParticleTransportStatisticsOptimizerState::do_update(unsigned int) {
  const double fs_in_ns = 1.0E+6;
  Transporting pt(p_);
  double prev_z = pt.get_last_tracked_z();
  double cur_z = core::XYZ(p_).get_coordinates()[2];
  pt.set_last_tracked_z(cur_z);  // save for RMF, etc.
  if (is_reset_) {               // ignore previous z if reset
    prev_z = cur_z;
    is_reset_ = false;
  }
  // update transport events
  if (cur_z > top_z_ && prev_z <= top_z_ && pt.get_n_entries_bottom() > 0 &&
      !pt.get_is_last_entry_from_top()) {
    n_transports_up_++;
    double time_ns = 0.0;
    if (owner_ != nullptr) {
      time_ns = owner_->get_current_time() / fs_in_ns;
    }
    transport_time_points_in_ns_.push_back(time_ns);
    IMP_LOG(PROGRESS, "EXIT UP " << p_->get_name() << " n_transports_up = "
            << n_transports_up_ << " prev_z = " << prev_z << std::endl);
  }
  if (cur_z < bottom_z_ && prev_z >= bottom_z_ && pt.get_n_entries_top() > 0 &&
      pt.get_is_last_entry_from_top()) {
    n_transports_down_++;
    double time_ns = 0.0;
    if (owner_ != nullptr) {
      time_ns = owner_->get_current_time() / fs_in_ns;
    }
    transport_time_points_in_ns_.push_back(time_ns);
    IMP_LOG(PROGRESS,"EXIT DOWN " << p_->get_name() << " n_transports_down = "
            << n_transports_down_ << " prev_z = " << prev_z << std::endl);
  }
  // update channel entry directionality
  if (cur_z < top_z_ && prev_z >= top_z_) {
    pt.set_is_last_entry_from_top(true);
    pt.set_n_entries_top(pt.get_n_entries_top() + 1);
  }
  if (cur_z > bottom_z_ && prev_z <= bottom_z_) {
    pt.set_is_last_entry_from_top(false);
    pt.set_n_entries_bottom(pt.get_n_entries_bottom() + 1);
  }
}

IMPNPCTRANSPORT_END_NAMESPACE
