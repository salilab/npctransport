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
  void print_score_and_positions(SimulationData *sd,
                                 bool print_positions = false,
                                 std::string header = "Score = ");

  void print_score_and_positions(SimulationData *sd,
                                 bool print_positions,
                                 std::string header) {
    std::cout << header
              << sd->get_bd()->get_scoring_function()->evaluate(false)
              << " ; PredicatePairsRestraint score = "
              << sd->get_predr()->evaluate(false)
              <<  std::endl;
    if(print_positions){
      ParticlesTemp ps = sd->get_diffusers()->get_particles();
      for(unsigned int i = 0; i < ps.size(); i++){
        std::cout << ps[i] << ", "
                  << IMP::core::RigidBody(ps[i]).get_reference_frame()
                  << std::endl;
      }
    }
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
     @param first_only some debug mode parameter for having a real short
                        simulation

     @return true if succesful, false if terminated abnormally
  */
  bool run_it(SimulationData *sd,
              unsigned int number_of_frames,
              boost::timer& timer,
              boost::timer& total_time,
              bool silent_statistics = false,
              unsigned int max_frames_per_chunk = 50000,
              bool first_only=false) {
    // TODO: next line is a temporary hack - needed for some reason to
    // force the pair predicates to evaluate predicate pairs restraints
    sd->get_m()->update();
    do {
      unsigned int cur_nframes
        = std::min<unsigned int>(first_only
                                 ?max_frames_per_chunk/10:max_frames_per_chunk,
                                 number_of_frames);
      //IMP_THREADS((sd, silent_statistics, cur_nframes),{
          std::cout << "Optimizing for " << cur_nframes
                    << " frames in this iteration" << std::endl;
          sd->get_bd()->optimize(cur_nframes);
          print_score_and_positions(sd);
          if(! silent_statistics) {
            sd->update_statistics(timer, cur_nframes);
          }
          std::cout << "Done" << std::endl;
          //});
      if (sd->get_maximum_number_of_minutes() > 0
          && total_time.elapsed()/60 > sd->get_maximum_number_of_minutes()) {
        sd->set_interrupted(true);
        std::cout << "Terminating..." << std::endl;
        return false;
      }
      number_of_frames-=cur_nframes;
    } while (number_of_frames > 0 && !first_only);
    return true;
  }
}

/**
   is_initial_optimization - whether to do initial optimization of positions
   is_full_run - whetehr to run the full BD simulation
 */
void do_main_loop(SimulationData *sd,
                  const RestraintsTemp &init_restraints,
                  bool quick,
                  bool is_initial_optimization,
                  bool is_equilibration,
                  bool is_full_run,
                  std::string final_conformations,
                  bool debug,
                  bool first_only) {
  using namespace IMP;
  const int max_frames_per_chunk=50000;
  try {
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
      if (is_initial_optimization) {
        std::cout<< "Doing initial coordinates optimization..." << std::endl;
        initialize_positions(sd, init_restraints, debug);
        sd->get_bd()->set_current_time( 0.0 );
      }
      print_score_and_positions( sd, debug, "Score right before BD = " );
      if (conformations_rmf_sos) {
        conformations_rmf_sos->update_always("Right before BD");
      }
      /*IMP::benchmark::Profiler p;
        if(i == 0)
        p.set("profiling.pprof");*/
      sd->get_bd()->set_log_level(IMP::PROGRESS);
      unsigned int nframes_run = (unsigned int)
        ( sd->get_number_of_frames() * sd->get_statistics_fraction() );
      unsigned int nframes_equilibrate =
        sd->get_number_of_frames() - nframes_run;
      if(is_equilibration){
        std::cout << "Equilibrating for " << nframes_equilibrate
                  << " frames..." << std::endl;
        bool ok = run_it
          (sd, nframes_equilibrate, timer, total_time,
           true /* silent stats */, max_frames_per_chunk, first_only);
        if(! ok || first_only)
          return;
        //if(nframes_equilibrate > 0) {
        // if equilibrated, ignore equilibration stats
        // TODO: removed for now since this may be incosistent with
        //       consecutive runs
        //sd->reset_statistics_optimizer_states();
        // }
        sd->get_bd()->set_current_time( 0.0 );
        std::cout << "Equilibration finished succesfully" << std::endl;
      }
      if(is_full_run) {
        timer.restart();
        std::cout << "Running for " << nframes_run << " frames..." << std::endl;
        if (conformations_rmf_sos) {
          conformations_rmf_sos->update_always
            ("Before running (post equilibration)");
        }
        // now run the rest of the sim
        bool ok = run_it(sd, nframes_run, timer, total_time,
                         false /* silent stats */, max_frames_per_chunk, first_only);
        if (! ok){
          return;
        }
        std::cout << "Run trial #" << i << " finished succesfully" << std::endl;
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
    }
    std::cout << "Entire run finished" << std::endl;
    print_score_and_positions( sd, debug, "Final score = " );
  }
  catch (IMP::base::Exception& e){
    std::cout << "Internal do_main_loop() function aborted with exception:\n\t"
              << e.what() << std::endl;
  }
}
IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE
