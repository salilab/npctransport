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
#include <IMP/log.h>
//#include <IMP/benchmark/Profiler.h>
#include <IMP/npctransport/initialize_positions.h>
#include <IMP/rmf/frames.h>
#include <IMP/ScoringFunction.h>
#include <IMP/Model.h>
#include <IMP/core/rigid_bodies.h>

IMPNPCTRANSPORT_BEGIN_INTERNAL_NAMESPACE
namespace {
  // TODO: move to H file?
  //! print this score for the current state of sd
  void print_score(SimulationData *sd, std::string header = "Score = ");

  void print_score(SimulationData *sd, std::string header) {
          std::cout << header
                    << sd->get_bd()->get_scoring_function()->evaluate(false)
                    << " ; PredicatePairsRestraint score = "
                    << sd->get_predr()->evaluate(false)
                    <<  std::endl;
  }

  /**
     Run simulation <sd> for <number_of_frames> frames, in chunks of
     optimization that last <max_frames_per_chunks> frames each.
     Statistics are updated after each chunk of optimization, using
     <timer> to time the current simulation.

     @param sd simulation data used for optimization
     @param number_of_frames total number of simulation frames requires
     @param timer a timer that was reset before this simulation trial was
                  initialized, to be used for tracking statistics
     @param total_time the total time that the simulation has spent.
                       The simulation would terminate at the end of an
                       optimization chunk in which this time has
                       elapsed, that is sd->get_maximum_number_of_minutes()
     @param silent_statistics if true, do not update statistics file
                              (e.g., during equilibration)
     @param max_frames_per_chunk maximal number of frames to be simulated
                                 in a single optimization chunk

     @return true if succesful, false if terminated abnormally
  */
  bool run_it(SimulationData *sd,
              unsigned int number_of_frames,
              boost::timer& timer,
              boost::timer& total_time,
              bool silent_statistics = false,
              unsigned int max_frames_per_chunk = 50000) {
    // TODO: next line is a temporary hack - needed for some reason to
    // force the pair predicates to evaluate predicate pairs restraints
    sd->get_m()->update();
    do {
      unsigned int cur_nframes
        = std::min<unsigned int>(max_frames_per_chunk,
                                 number_of_frames);
#pragma omp parallel num_threads(3)
      {
#pragma omp single
        {
          std::cout << "Optimizing for " << cur_nframes
                    << " frames in this iteration" << std::endl;
          sd->get_bd()->optimize(cur_nframes);
          print_score(sd);
          if(! silent_statistics) {
            sd->update_statistics(timer, cur_nframes);
          }
        }
      }
      if (sd->get_maximum_number_of_minutes() > 0
          && total_time.elapsed()/60 > sd->get_maximum_number_of_minutes()) {
        sd->set_interrupted(true);
        std::cout << "Terminating..." << std::endl;
        return false;
      }
      number_of_frames-=cur_nframes;
    } while (number_of_frames > 0);
    return true;
  }
}

void do_main_loop(SimulationData *sd,
                  const RestraintsTemp &init_restraints,
                  bool quick, bool init_only, std::string final_conformations,
                  bool debug_initialize, std::string init_rmf) {
  using namespace IMP;
  IMP::Pointer<rmf::SaveOptimizerState> conformations_rmf_sos
    = sd->get_rmf_sos_writer();
  RMF::FileHandle final_rmf_fh;
  if(!final_conformations.empty()){
    final_rmf_fh=RMF::create_rmf_file(final_conformations);
    sd->link_rmf_file_handle(final_rmf_fh);
  }
  boost::timer total_time;
  for (unsigned int i=0; i< sd->get_number_of_trials(); ++i) {
    IMP::base::CreateLogContext clc("iteration");
    boost::timer timer;
    IMP::set_log_level(SILENT);
    std::cout << "Simulation trial " << i << " out of "
              << sd->get_number_of_trials() << std::endl;
    if (!quick)
      sd->reset_rmf();
    std::cout<< "Initializing...";
    if (init_rmf == "") {
 	  initialize_positions(sd, init_restraints, debug_initialize);
          std::cout << " from scratch" << std::endl;
    }
    else{
      sd->initialize_positions_from_rmf(init_rmf, -1);
      std::cout << " from last frame of existing RMF file " << init_rmf << std::endl;
    }
    print_score( sd, "Initial score = " );
    if (debug_initialize) break;
    sd->get_bd()->set_log_level(IMP::PROGRESS);
    if (conformations_rmf_sos) {
      conformations_rmf_sos->update_always("After initialization");
    }
    {
      ParticlesTemp ps = sd->get_diffusers()->get_particles();
      for(unsigned int i = 0; i < ps.size(); i++){
        std::cout << ps[i] << ", "
          //        <<  IMP::core::XYZ(ps[i]).get_coordinates()
                  << IMP::core::RigidBody(ps[i]).get_reference_frame()
                  << std::endl;
      }
    }
    /*IMP::benchmark::Profiler p;
    if(i == 0)
      p.set("profiling.pprof");*/
    sd->get_bd()->set_current_time(0);
    {
      double equilibrate_fraction =  1.0 - sd->get_statistics_fraction() ;
      unsigned int nframes_equilibrate = (unsigned int)
        (sd->get_number_of_frames() * equilibrate_fraction );
      std::cout << "Equilibrating for " << nframes_equilibrate
                << " frames..." << std::endl;
      bool ok = run_it
        (sd, nframes_equilibrate, timer, total_time,
         true /* silent stats */);
      if(! ok)
        return;
    }
    std::cout << "Equilibration finished succesfully" << std::endl;
    if (init_only) {
      continue; // skip optimization
    }
    {
      timer.restart();
      sd->reset_statistics_optimizer_states();
      unsigned int nframes_run = (unsigned int)
        ( sd->get_number_of_frames() * sd->get_statistics_fraction() );
      std::cout << "Running for " << nframes_run << " frames..." << std::endl;
      if (conformations_rmf_sos) {
        conformations_rmf_sos->update_always("Before running (post equilibration)");
      }
      // now run the rest of the sim
      bool ok = run_it(sd, nframes_run, timer, total_time,
                       false /* silent stats */);
      if (! ok){
        return;
      }
      if (conformations_rmf_sos) {
        conformations_rmf_sos->update_always("Final frame");
      }
      if( !final_conformations.empty() ) {
        std::cout << "Printing last frame to "
                  << final_conformations << std::endl;
        IMP::rmf::save_frame
          ( final_rmf_fh, final_rmf_fh.get_number_of_frames() );
      }
      std::cout << "Run trial #" << i << " finished succesfully" << std::endl;
    }
  }
  std::cout << "Entire run finished" << std::endl;
  print_score( sd, "Final score = " );
  {
    ParticlesTemp ps = sd->get_diffusers()->get_particles();
    for(unsigned int i = 0; i < ps.size(); i++){
      std::cout << ps[i] << ", "
        //      << IMP::core::XYZ(ps[i]).get_coordinates()
                << IMP::core::RigidBody(ps[i]).get_reference_frame()
                << std::endl;
    }
  }
}
IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE
