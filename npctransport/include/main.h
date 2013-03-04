/**
 *  \file main.h
 *  \brief Helper functions for executable .cpp files
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_MAIN_H
#define IMPNPCTRANSPORT_MAIN_H


#ifdef IMP_NPC_MAIN
#include "npctransport_config.h"
#include <IMP/algebra.h>
#include <IMP/core.h>
#include <IMP/atom.h>
#include <IMP/display.h>
#include <IMP/rmf.h>
#include <IMP.h>
#include <IMP/container.h>
#include <IMP/base/CreateLogContext.h>
#include <IMP/base/random.h>
#include "SimulationData.h"
#ifdef IMP_NPC_GOOGLE
IMP_GCC_PUSH_POP(diagnostic push)
IMP_GCC_PRAGMA(diagnostic ignored "-Wsign-compare")
#include "third_party/npc/npctransport/data/npctransport.pb.h"
IMP_GCC_PUSH_POP(diagnostic pop)
#else
#include <IMP/npctransport/internal/npctransport.pb.h>
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


#include "internal/main.h"
#include "protobuf.h"

#ifdef IMP_NPC_GOOGLE
#include "internal/google_main.h"
#else
#include "internal/boost_main.h"
#endif


/**
   reads all command line parameters defined in the previous invokations of
   IMP_NPC_PARAMETERS_XXX macros. It then initializes sim_data as a SimulationData
   objects, based on these command line paramateres.
 */
#define IMP_NPC_STARTUP(sim_data)                                       \
  IMP::base::Pointer<IMP::npctransport::SimulationData> sim_data        \
  = startup(argc, argv);                                                \
  if (!sim_data) return 1




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
IMP_NPC_PARAMETER_UINT64(random_seed, 0,
                         "an unsigned integer to be used as a random seed"
                         " in the IMP random numbers generator. If unspecified"
                         " or zero, IMP will use system time as a seed."
                         " (Note: if in doubt, use a 32 bit unsigned integer"
                         " seed in the range 0 to 4,294,967,295)");

/** Run simulation using preconstructed SimulationData object sim_data.
    init_restraints are used ad-hoc during initialization only,
    and unless initialized from an RMF file
*/
#define IMP_NPC_LOOP(sim_data, init_restraints)                         \
  {                                                                     \
  /** initial optimization and equilibration needed unless starting
      from another output file or rmf file */                           \
    bool is_initial_optimization =                                      \
      FLAGS_restart.empty() && FLAGS_init_rmffile.empty();              \
    bool is_BD_equilibration = is_initial_optimization;                 \
    bool is_BD_full_run = !FLAGS_initialize_only;                       \
    IMP::npctransport::internal::do_main_loop(sim_data,                 \
                                              init_restraints,          \
                                              FLAGS_quick,              \
                                              is_initial_optimization,  \
                                              is_BD_equilibration,      \
                                              is_BD_full_run,           \
                                              FLAGS_final_conformations, \
                                              FLAGS_verbose,            \
                                              FLAGS_first_only);        \
  }

namespace {
//! seeds the random number generator of IMP with seed
//! (or time if it is zero)
/**
   returns the actual seed used to initialize the generator
 */
inline std::size_t seed_randn_generator(std::size_t seed)
{
  if(seed == 0){
    seed =  static_cast<std::size_t> (std::time(0)); //IMP::nullptr)) ;
  }
  IMP::base::random_number_generator.seed( seed );
  return seed;
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

/**
   initialize and return a simulation data object based on
   program command line parameters
 */
inline IMP::npctransport::SimulationData *startup(int argc, char *argv[]) {
  IMP_NPC_PARSE_OPTIONS(argc, argv);
  boost::uint64_t actual_seed =
    seed_randn_generator( FLAGS_random_seed );
#pragma omp critical
  std::cout << "Random seed is " << actual_seed << std::endl;
  set_log_level(IMP::base::LogLevel(FLAGS_log_level));
  IMP::base::Pointer<IMP::npctransport::SimulationData> sd;
  try {
    write_output_based_on_flags( actual_seed );
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
  } catch (const IMP::base::IOException &e) {
#pragma omp critical
    std::cerr << "Error: " << e.what() << std::endl;
    return IMP_NULLPTR;
  }
  return sd.release();
}
}
#endif // IMP_NPC_MAIN

#endif /* IMPNPCTRANSPORT_MAIN_H */
