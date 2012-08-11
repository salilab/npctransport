/**
 *  \file io.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/internal/main.h>
#include <IMP/base_types.h>
#include <boost/timer.hpp>
//#include <IMP/benchmark/Profiler.h>
#include <IMP/npctransport/initialize_positions.h>

IMPNPCTRANSPORT_BEGIN_INTERNAL_NAMESPACE
namespace {
  bool run_it(SimulationData *sd,
              unsigned int number_of_frames,
              boost::timer&total_time) {
    do {
      unsigned int cur_frames
          = std::min<unsigned int>(1000000,
                                   number_of_frames);
      std::cout << "Running..." << std::endl;
      sd->get_bd()->optimize(cur_frames);
      if (sd->get_maximum_number_of_minutes() > 0
          && total_time.elapsed()/60 > sd->get_maximum_number_of_minutes()) {
        sd->set_interrupted(true);
        std::cout << "Terminating..." << std::endl;
        return true;
      }
      number_of_frames-=cur_frames;
    } while (number_of_frames > 0);
    return false;
  }
}

void do_main_loop(SimulationData *sd, const ParticlePairsTemp &links,
                  bool quick, std::string final_config,
                  bool debug_initialize, std::string init_rmf) {
  using namespace IMP;
  base::Pointer<rmf::SaveOptimizerState> final_sos;
  if (!final_conformations.empty()) {
    final_sos = sd->create_rmf_writer(final_conformations);
  }
  sd->set_was_used(true);
  boost::timer total_time;
  for (unsigned int i=0; i< sd->get_number_of_trials(); ++i) {
    IMP::base::CreateLogContext clc("iteration");
    boost::timer timer;
    IMP::set_log_level(SILENT);
    if (!quick) sd->reset_rmf();
    std::cout<< "Initializing..." << std::endl;
    if (init_rmf == "")
      initialize_positions(sd, links, debug_initialize);
    else{
      sd->initialize_positions_from_rmf(init_rmf);
      std::cout << "Initializing positions from RMF file "
                << init_rmf << std::endl;
    }
    if (debug_initialize) break;
    sd->get_bd()->set_log_level(IMP::PROGRESS);
    if (final_sos) {
      final_sos->update_always();
    }
    /*IMP::benchmark::Profiler p;
    if(i == 0)
      p.set("profiling.pprof");*/
    sd->get_bd()->set_current_time(0);
    std::cout << "Equilibrating..." << std::endl;
    if (run_it(sd, sd->get_number_of_frames()
               * sd->get_statistics_fraction(), total_time)) {
      return;
    }
    sd->reset_statistics_optimizer_states();
    std::cout << "Running..." << std::endl;
    // now run the rest of the sim
    sd->get_bd()->optimize(sd->get_number_of_frames()
                           *(1.0- sd->get_statistics_fraction()));
    bool abort=run_it(sd, sd->get_number_of_frames()
                      * (1.0- sd->get_statistics_fraction()), total_time);
    //p.reset();
    sd->update_statistics(timer);
    if (final_sos) {
      final_sos->update_always();
    }
    if (abort) break;
  }
}
IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE
