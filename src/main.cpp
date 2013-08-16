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
#include <IMP/base/log.h>
#include <IMP/base/exception.h>
#include <IMP/base/check_macros.h>
//#include <IMP/benchmark/Profiler.h>
#include <IMP/npctransport/initialize_positions.h>
#include <IMP/rmf/frames.h>
#include <IMP/ScoringFunction.h>
#include <IMP/Model.h>
#include <IMP/core/rigid_bodies.h>

//#include <IMP/algebra.h>
//#include <IMP/core.h>
//#include <IMP/atom.h>
//#include <IMP/display.h>
//#include <IMP/rmf.h>
//#include <IMP.h>
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
//#include <IMP/example/creating_restraints.h>
//#include <IMP/example/randomizing.h>
//#include <IMP/example/counting.h>
//#include <IMP/example/optimizing.h>

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

IMPNPCTRANSPORT_BEGIN_NAMESPACE
boost::int64_t work_unit = -1;
IMP::base::AddIntFlag work_unitadder
( "work_unit", "The work unit",
  &work_unit);
std::string configuration = "configuration.pb";
IMP::base::AddStringFlag configuration_adder
( "configuration", "input configuration file in protobuf format"
  " [default: %default]",
  &configuration);
std::string output = "output.pb";
IMP::base::AddStringFlag output_adder
( "output", "output assignments and statistics file in protobuf format,"
  " recording the assignment being executed"
  " [default: %default]",
  &output);
std::string restart = "";
IMP::base::AddStringFlag restart_adder
( "restart", "output file of a previous run, from which to restart"
  " this run (also initializing the final coordinates"
  " from this previous run, if they exist in the output"
  " file)"
  " [default: %default]",
  &restart);
std::string conformations = "conformations.rmf";
IMP::base::AddStringFlag conformations_adder
( "conformations", "RMF file for recording the conforomations along the "
  " simulation [default: %default]",
  &conformations);
std::string init_rmffile = "";
IMP::base::AddStringFlag init_rmf_adder
( "init_rmffile", "[OBSOLETE]"
  " RMF file for initializing the simulation with its"
  " last frame (to continue a previous run). Note that"
  " this option overrides the coordinates specified in"
  " an older output file when using the --restart flag",
  &init_rmffile);
std::string final_conformations = "final_conformations.rmf";
IMP::base::AddStringFlag final_conformations_adder
( "final_conformations",
  "RMF file for recording the initial and final conformations "
  " [default: %default]",
  &final_conformations);
bool verbose = false;
IMP::base::AddBoolFlag verbose_adder
( "verbose", "Print more info during run",
  &verbose);
bool first_only = false;
IMP::base::AddBoolFlag first_only_adder
( "first_only",
  "Only do the first simulation block",
  &first_only);
bool initialize_only = false;
IMP::base::AddBoolFlag initialize_only_adder
( "initialize_only", "Run the initialization and then stop",
  &initialize_only);
bool show_steps = false;
IMP::base::AddBoolFlag show_steps_adder
( "show_steps", "Show the steps for each modified variable",
  &show_steps);
bool show_number_of_work_units = false;
IMP::base::AddBoolFlag show_work_units_adder
( "show_number_of_work_units",
  "Show the number of work units" ,
  &show_number_of_work_units);
double short_init_factor = 1.0;
base::AddFloatFlag short_init_adder
( "short_init_factor",
  "Run an abbreviated version of system initialization, which takes"
  " a fraction of a full initialization, in the range (0.0..1.0]"
  " [default=1.0]",
  &short_init_factor);
double short_sim_factor = 1.0;
base::AddFloatFlag short_sim_adder
( "short_sim_factor",
  "Run an abbreviated version of the simulation, which takes"
  " a fraction of a full simulation (or more if >1.0)"
  " [default=1.0]",
  &short_sim_factor);


namespace {
  /*********************************** internal functions
   * *********************************/

  // TODO: move to H file?
  //! print this score for the current state of sd
  void print_score_and_positions(SimulationData *sd, bool print_positions = false,
                                 std::string header = "Score = ");

  void print_score_and_positions(SimulationData *sd, bool print_positions,
                                 std::string header) {
    IMP_OMP_PRAGMA(critical)
      std::cout << header
                << sd->get_bd()->get_scoring_function()->evaluate(false)
                << " ; PredicatePairsRestraint score = "
                << sd->get_scoring()
                    ->get_predicates_pair_restraint()->evaluate(false)
                << std::endl;
    if (print_positions) {
      ParticlesTemp ps = sd->get_diffusers()->get_particles();
      for (unsigned int i = 0; i < ps.size(); i++) {
        IMP_OMP_PRAGMA(critical)
          std::cout << ps[i] << ", " << IMP::core::RigidBody(ps[i])
          .get_reference_frame() << std::endl;
      }
    }
  }

  /** writes the output assignment file based on the configuration parameters
      either assign a new work unit from a configuration file, or restart
      from a previous output file

      @param actual_seed the actual random seed used in the simulation,
      to be saved in the output file
  */
  inline void write_output_based_on_flags(boost::uint64_t actual_seed) {
    if (restart.empty()) {
      int num = IMP::npctransport::assign_ranges(configuration, output, work_unit,
                                               show_steps, actual_seed);
      if (show_number_of_work_units) {
      IMP_OMP_PRAGMA(critical)
        std::cout << "work units " << num << std::endl;
    }
    } else {  // resart.empty()
    IMP_OMP_PRAGMA(critical)
      std::cout << "Restart simulation from " << restart << std::endl;
    ::npctransport_proto::Output prev_output;
    // copy to new file to avoid modifying input file
    std::ifstream file(restart.c_str(), std::ios::binary);
    bool read = prev_output.ParseFromIstream(&file);
    IMP_ALWAYS_CHECK(read, "Couldn't read restart file " << restart,
                     IMP::base::ValueException);
    prev_output.mutable_assignment()->set_random_seed(actual_seed);
    std::ofstream outf(output.c_str(), std::ios::binary);
    prev_output.SerializeToOstream(&outf);
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
     @param silent_statistics if true, do not update statistics file (e.g.,
                              during equilibration)
     @param max_frames_per_chunk maximal number of frames to be simulated
                                 in a single optimization chunk

     @note if first_only==true, make a really short simulation

     @return true if succesful, false if terminated abnormally
  */
  bool run_it(SimulationData *sd, unsigned int number_of_frames,
              boost::timer &timer, boost::timer &total_time,
              bool silent_statistics = false,
              unsigned int max_frames_per_chunk = 50000) {
    // TODO: next line is a temporary hack - needed for some reason to
    // force the pair predicates to evaluate predicate pairs restraints
    sd->get_model()->update();
    do {
      unsigned int cur_nframes = std::min<unsigned int>
        ( first_only ? max_frames_per_chunk / 10 : max_frames_per_chunk,
          number_of_frames);
      // IMP_THREADS((sd, silent_statistics, cur_nframes),{
      std::cout << "Optimizing for " << cur_nframes << " frames in this iteration"
                << std::endl;
      sd->get_bd()->optimize(cur_nframes);
      print_score_and_positions(sd);
      if (!silent_statistics) {
        sd->get_statistics()->update(timer, cur_nframes);
      }
      std::cout << "Done" << std::endl;
      //});
      if (sd->get_maximum_number_of_minutes() > 0 &&
          total_time.elapsed() / 60 > sd->get_maximum_number_of_minutes()) {
        sd->get_statistics()->set_interrupted(true);
        std::cout << "Terminating..." << std::endl;
        return false;
      }
      number_of_frames -= cur_nframes;
    } while (number_of_frames > 0 && !first_only);
    return true;
  }

} // anonymous namespace

/********************* public functions **********************/

// initialize and return a simulation data object based on
// program command line parameters
IMP::npctransport::SimulationData *startup(int argc, char *argv[]) {
  IMP_NPC_PARSE_OPTIONS(argc, argv);
  IMP_ALWAYS_CHECK( short_init_factor <= 1.0 && short_init_factor > 0,
                    "short_init_factor must be in the range (0..1.0]",
                    IMP::base::ValueException );
  IMP_OMP_PRAGMA(critical)
  std::cout << "Random seed is " << IMP::base::get_random_seed() << std::endl;
  IMP::base::Pointer<IMP::npctransport::SimulationData> sd;
  write_output_based_on_flags(IMP::base::get_random_seed());
  sd = new IMP::npctransport::SimulationData(output, IMP::base::run_quick_test);
  if (!conformations.empty()) {
    sd->set_rmf_file_name(conformations);
  }
  if (!init_rmffile.empty()) {
    sd->initialize_positions_from_rmf(
        RMF::open_rmf_file_read_only(init_rmffile), -1);
    IMP_OMP_PRAGMA(critical)
    std::cout
        << "Initialize coordinates from last frame of an existing RMF file "
        << init_rmffile << std::endl;
  }
  return sd.release();
}

//  Run simulation using preconstructed SimulationData object sd,
//  with ad-hoc init restratins init_restraints
void do_main_loop(SimulationData *sd, const RestraintsTemp &init_restraints) {
  using namespace IMP;
  sd->set_was_used( true );
  const int max_frames_per_chunk = 50000;
  /** initial optimization and equilibration needed unless starting
      from another output file or rmf file */
  bool is_initial_optimization = restart.empty() && init_rmffile.empty();
  bool is_BD_equilibration = is_initial_optimization;
  bool is_BD_full_run = !initialize_only;

  base::Pointer<rmf::SaveOptimizerState> conformations_rmf_sos =
      sd->get_rmf_sos_writer();
  RMF::FileHandle final_rmf_fh;
  if (!final_conformations.empty()) {
    final_rmf_fh = RMF::create_rmf_file(final_conformations);
    sd->link_rmf_file_handle(final_rmf_fh);
  }
  boost::timer total_time;
  for (unsigned int i = 0; i < sd->get_number_of_trials(); ++i) {
    IMP::base::CreateLogContext clc("iteration");
    boost::timer timer;
    //    IMP::base::set_log_level(IMP::base::PROGRESS);
    std::cout << "Simulation trial " << i << " out of "
              << sd->get_number_of_trials() << std::endl;
    if (is_initial_optimization) {
      sd->switch_suspend_rmf(true);
      std::cout << "Doing initial coordinates optimization..." << std::endl;
      initialize_positions(sd, init_restraints, verbose, short_init_factor);
      sd->get_bd()->set_current_time(0.0);
      sd->get_statistics()->reset_statistics_optimizer_states();
      sd->switch_suspend_rmf(false);
      //      conformation_rmf_sos->reset();
    }
    print_score_and_positions(sd, verbose, "Score right before BD = ");
    if (conformations_rmf_sos) {
      conformations_rmf_sos->update_always("Right before BD");
    }
    /*IMP::benchmark::Profiler p;
      if(i == 0)
      p.set("profiling.pprof");*/
    // sd->get_bd()->set_log_level(IMP::PROGRESS);
    unsigned int actual_nframes =
      (unsigned int)std::ceil(sd->get_number_of_frames() * short_sim_factor);
    unsigned int nframes_run = (unsigned int)(actual_nframes *
                                              sd->get_statistics_fraction());
    unsigned int nframes_equilibrate = actual_nframes - nframes_run;
    if (is_BD_equilibration) {
      std::cout << "Equilibrating for " << nframes_equilibrate << " frames..."
                << std::endl;
      bool ok = run_it(sd, nframes_equilibrate, timer, total_time,
                       true /* silent stats */, max_frames_per_chunk);
      if (!ok || first_only) return;
      // if(nframes_equilibrate > 0) {
      // if equilibrated, ignore equilibration stats
      // TODO: removed for now since this may be incosistent with
      //       consecutive runs
      // sd->reset_statistics_optimizer_states();
      // }
      sd->get_bd()->set_current_time(0.0);
      std::cout << "Equilibration finished succesfully" << std::endl;
    }
    if (is_BD_full_run) {
      timer.restart();
      std::cout << "Running for " << nframes_run << " frames..." << std::endl;
      if (conformations_rmf_sos) {
        conformations_rmf_sos->update_always(
            "Before running (post equilibration)");
      }
      // now run the rest of the sim
      bool ok = run_it(sd, nframes_run, timer, total_time,
                       false /* silent stats */, max_frames_per_chunk);
      if (!ok) {
        return;
      }
      std::cout << "Run trial #" << i << " finished succesfully" << std::endl;
    }
    if (conformations_rmf_sos) {
      conformations_rmf_sos->update_always("Final frame");
    }
    if (!final_conformations.empty()) {
      std::cout << "Printing last frame to " << final_conformations
                << std::endl;
      IMP::rmf::save_frame(final_rmf_fh, final_rmf_fh.get_number_of_frames());
    }
  }
  std::cout << "Entire run finished" << std::endl;
  print_score_and_positions(sd, verbose, "Final score = ");
}

IMPNPCTRANSPORT_END_NAMESPACE
