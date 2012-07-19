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
void do_main_loop(SimulationData *sd, const ParticlePairsTemp &links,
                  bool quick, std::string final_config) {
  using namespace IMP;
  for (unsigned int i=0; i< sd->get_number_of_trials(); ++i) {
    IMP::base::CreateLogContext clc("iteration");
    boost::timer timer;
    IMP::set_log_level(SILENT);
    if (!quick) sd->reset_rmf();
    std::cout<< "Initializing..." << std::endl;
    initialize_positions(sd, links);
    sd->get_bd()->set_log_level(SILENT);
    sd->get_bd()->set_log_level(IMP::PROGRESS);
    /*IMP::benchmark::Profiler p;
    if(i == 0)
      p.set("profiling.pprof");*/
    sd->get_bd()->set_current_time(0);
    std::cout << "Equilibrating..." << std::endl;
    sd->get_bd()->optimize(sd->get_number_of_frames()
                           * sd->get_statistics_fraction());
    sd->reset_statistics_optimizer_states();
    std::cout << "Running..." << std::endl;
    // now run the rest of the sim
    sd->get_bd()->optimize(sd->get_number_of_frames()
                           *(1.0- sd->get_statistics_fraction()));
    //p.reset();
    sd->update_statistics(timer);
    std::cout << "Writing..." << std::endl;
    sd->write_geometry(final_config);
  }
}
IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE
