/**
 *  \file io.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

// TODO: cleanup all headers - probably most are not needed

#include <IMP/npctransport/main.h>
#include <IMP/base_types.h>
#include <boost/timer.hpp>
#include <IMP/log.h>
//#include <IMP/benchmark/Profiler.h>
#include <IMP/npctransport/initialize_positions.h>
#include <IMP/rmf/frames.h>
#include <IMP/ScoringFunction.h>
#include <IMP/Model.h>
#include <IMP/core/rigid_bodies.h>

#include <IMP/algebra.h>
#include <IMP/core.h>
#include <IMP/atom.h>
#include <IMP/display.h>
#include <IMP/rmf.h>
#include <IMP.h>
#include <IMP/container.h>
#include <IMP/base/CreateLogContext.h>
#include <IMP/base/random.h>
#ifdef IMP_NPC_GOOGLE
IMP_GCC_PUSH_POP(diagnostic push)
IMP_GCC_PRAGMA(diagnostic ignored "-Wsign-compare")
#include "third_party/npc/npctransport/data/npctransport.pb.h"
IMP_GCC_PUSH_POP(diagnostic pop)
#include <IMP/npctransport/internal/google_main.h>
#else
#include <IMP/npctransport/internal/npctransport.pb.h>
#include <IMP/npctransport/internal/boost_main.h>
#endif
#include <IMP/npctransport.h>
//#include <IMP/benchmark/Profiler.h>

// use the example code for now to work bugs out of it
#include <IMP/example/creating_restraints.h>
#include <IMP/example/randomizing.h>
#include <IMP/example/counting.h>
#include <IMP/example/optimizing.h>

#include <IMP/base/nullptr.h>
#include <IMP/base/nullptr_macros.h>
#include <IMP/base/check_macros.h>
#include <IMP/base/exception.h>
#include <boost/cstdint.hpp>
#include <numeric>
#include <cmath>
#include <iostream>
#include <ctime>

#include <IMP/npctransport/protobuf.h>


IMP_NPC_PARAMETER_INT64(work_unit, -1, "The work unit");
IMP_NPC_PARAMETER_INT64(log_level, 0, "The log level to use");
IMP_NPC_PARAMETER_STRING(configuration, "configuration.pb",
                         "input configuration file in protobuf format"
                         " [default: %default]");
IMP_NPC_PARAMETER_STRING(output, "output.pb",
                         "output assignments and statistics file in protobuf format,"
                         " recording the assignment being executed"
                         " [default: %default]");
IMP_NPC_PARAMETER_STRING(restart, "",
                         "output file of a previous run, from which to restart"
                         " this run (also initializing the final coordinates"
                         " from this previous run, if they exist in the output"
                         " file)"
                         " [default: %default]");
IMP_NPC_PARAMETER_STRING(conformations, "conformations.rmf",
                         "RMF file for recording the conforomations along the "
                         " simulation [default: %default]");
IMP_NPC_PARAMETER_STRING(init_rmffile, "",
                         "[OBSOLETE]"
                         " RMF file for initializing the simulation with its"
                         " last frame (to continue a previous run). Note that"
                         " this option overrides the coordinates specified in"
                         " an older output file when using the --restart flag");
IMP_NPC_PARAMETER_STRING(final_conformations, "final_conformations.rmf",
                         "RMF file for recording the initial and final conformations "
                         " [default: %default]");
IMP_NPC_PARAMETER_BOOL(verbose, false,
                       "Print more info during run");
IMP_NPC_PARAMETER_BOOL(first_only, false,
                       "Only do the first simulation block");
IMP_NPC_PARAMETER_BOOL(initialize_only, false,
                       "Run the initialization and then stop");
IMP_NPC_PARAMETER_BOOL(quick, false,
                       "Reduce all steps to the minimum");
IMP_NPC_PARAMETER_BOOL(show_steps, false,
                       "Show the steps for each modified variable");
IMP_NPC_PARAMETER_BOOL(show_number_of_work_units, false,
                       "Show the number of work units");


IMPNPCTRANSPORT_BEGIN_NAMESPACE

/*********************************** internal functions *********************************/

// TODO: move to H file?
//! print this score for the current state of sd
void print_score_and_positions(SimulationData *sd,
                                 bool print_positions = false,
                                 std::string header = "Score = ");

void print_score_and_positions(SimulationData *sd,
                               bool print_positions,
                               std::string header) {
#pragma omp critical
  std::cout << header
            << sd->get_bd()->get_scoring_function()->evaluate(false)
            << " ; PredicatePairsRestraint score = "
            << sd->get_predr()->evaluate(false)
            <<  std::endl;
  if(print_positions){
    ParticlesTemp ps = sd->get_diffusers()->get_particles();
      for(unsigned int i = 0; i < ps.size(); i++){
#pragma omp critical
        std::cout << ps[i] << ", "
                  << IMP::core::RigidBody(ps[i]).get_reference_frame()
                  << std::endl;
      }
    }
}

/** writes the output assignment file based on the configuration parameters
    either assign a new work unit from a configuration file, or restart
    from a previous output file

    @param actual_seed the actual random seed used in the simulation,
                       to be saved in the output file
*/
inline void write_output_based_on_flags( boost::uint64_t actual_seed ) {
  if (FLAGS_restart.empty()) {
    int num=IMP::npctransport::assign_ranges
      (FLAGS_configuration, FLAGS_output,
       FLAGS_work_unit, FLAGS_show_steps, actual_seed);
    if (FLAGS_show_number_of_work_units) {
#pragma omp critical
      std::cout << "work units " << num << std::endl;
    }
  } else { // FLAGS_resart.empty()
#pragma omp critical
    std::cout << "Restart simulation from " << FLAGS_restart << std::endl;
      ::npctransport_proto::Output prev_output;
    // copy to new file to avoid modifying input file
    std::ifstream file(FLAGS_restart.c_str(), std::ios::binary);
    bool read=prev_output.ParseFromIstream(&file);
    IMP_ALWAYS_CHECK(read, "Couldn't read restart file " << FLAGS_restart,
                     IMP::base::ValueException);
    prev_output.mutable_assignment()->set_random_seed( actual_seed );
    std::ofstream outf(FLAGS_output.c_str(), std::ios::binary);
    prev_output.SerializeToOstream(&outf);
  }
}

/********************* public functions **********************/

// initialize and return a simulation data object based on
// program command line parameters
IMP::npctransport::SimulationData *startup(int argc, char *argv[]) {
  IMP_NPC_PARSE_OPTIONS(argc, argv);
#pragma omp critical
  std::cout << "Random seed is " << IMP::base::get_random_seed() << std::endl;
  set_log_level(IMP::base::LogLevel(FLAGS_log_level));
  IMP::base::Pointer<IMP::npctransport::SimulationData> sd;
  write_output_based_on_flags( IMP::base::get_random_seed() );
  sd= new IMP::npctransport::SimulationData(FLAGS_output,
                                            FLAGS_quick);
    if (!FLAGS_conformations.empty()) {
      sd->set_rmf_file_name(FLAGS_conformations);
    }
    if(!FLAGS_init_rmffile.empty()) {
      sd->initialize_positions_from_rmf(RMF::open_rmf_file_read_only(FLAGS_init_rmffile),
                                        -1);
#pragma omp critical
      std::cout << "Initialize coordinates from last frame of an existing RMF file "
                << FLAGS_init_rmffile << std::endl;
    }
    return sd.release();
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
   @note if FLAGS_first_only==true, make a really short simulation

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
      = std::min<unsigned int>(FLAGS_first_only
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
  } while (number_of_frames > 0 && !FLAGS_first_only);
  return true;
}

//  Run simulation using preconstructed SimulationData object sd,
//  with ad-hoc init restratins init_restraints
void do_main_loop(SimulationData *sd,
                  const RestraintsTemp &init_restraints)
{
  using namespace IMP;
  const int max_frames_per_chunk=50000;
  /** initial optimization and equilibration needed unless starting
      from another output file or rmf file */
  bool is_initial_optimization =
      FLAGS_restart.empty() && FLAGS_init_rmffile.empty();
  bool is_BD_equilibration = is_initial_optimization;
  bool is_BD_full_run = !FLAGS_initialize_only;

  IMP::Pointer<rmf::SaveOptimizerState> conformations_rmf_sos
    = sd->get_rmf_sos_writer();
  RMF::FileHandle final_rmf_fh;
  if(!FLAGS_final_conformations.empty()){
    final_rmf_fh=RMF::create_rmf_file(FLAGS_final_conformations);
    sd->link_rmf_file_handle(final_rmf_fh);
  }
  boost::timer total_time;
  for (unsigned int i=0; i< sd->get_number_of_trials(); ++i) {
    IMP::base::CreateLogContext clc("iteration");
    boost::timer timer;
    IMP::set_log_level(SILENT);
    std::cout << "Simulation trial " << i << " out of "
              << sd->get_number_of_trials() << std::endl;
    if (is_initial_optimization) {
      std::cout<< "Doing initial coordinates optimization..." << std::endl;
      initialize_positions(sd, init_restraints, FLAGS_verbose);
      sd->get_bd()->set_current_time( 0.0 );
    }
    print_score_and_positions( sd, FLAGS_verbose, "Score right before BD = " );
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
    if(is_BD_equilibration){
      std::cout << "Equilibrating for " << nframes_equilibrate
                << " frames..." << std::endl;
      bool ok = run_it
        (sd, nframes_equilibrate, timer, total_time,
         true /* silent stats */, max_frames_per_chunk);
      if(! ok || FLAGS_first_only)
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
    if(is_BD_full_run) {
      timer.restart();
      std::cout << "Running for " << nframes_run << " frames..." << std::endl;
      if (conformations_rmf_sos) {
        conformations_rmf_sos->update_always
          ("Before running (post equilibration)");
      }
      // now run the rest of the sim
      bool ok = run_it(sd, nframes_run, timer, total_time,
                       false /* silent stats */, max_frames_per_chunk);
      if (! ok){
        return;
      }
      std::cout << "Run trial #" << i << " finished succesfully" << std::endl;
    }
    if (conformations_rmf_sos) {
      conformations_rmf_sos->update_always("Final frame");
    }
    if( !FLAGS_final_conformations.empty() ) {
      std::cout << "Printing last frame to "
                << FLAGS_final_conformations << std::endl;
      IMP::rmf::save_frame
        ( final_rmf_fh, final_rmf_fh.get_number_of_frames() );
    }
  }
  std::cout << "Entire run finished" << std::endl;
  print_score_and_positions( sd, FLAGS_verbose, "Final score = " );
}

IMPNPCTRANSPORT_END_NAMESPACE
