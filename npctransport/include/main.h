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
#include <IMP/npctransport.h>
//#include <IMP/benchmark/Profiler.h>

// use the example code for now to work bugs out of it
#include <IMP/example/creating_restraints.h>
#include <IMP/example/randomizing.h>
#include <IMP/example/counting.h>
#include <IMP/example/optimizing.h>

#include <IMP/compatibility/nullptr.h>
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




IMP_NPC_PARAMETER_INT(work_unit, -1, "The work unit");
IMP_NPC_PARAMETER_INT(log_level, 0, "The log level to use");
IMP_NPC_PARAMETER_STRING(configuration, "configuration.pb",
                         "input configuration file in protobuf format"
                         " [default: %default]");
IMP_NPC_PARAMETER_STRING(output, "output.pb",
                         "output assignments and statistics file in protobuf format,"
                         " recording the assignment being executed"
                         " [default: %default]");
IMP_NPC_PARAMETER_STRING(conformations, "conformations.rmf",
                         "RMF file for recording the conforomations along the "
                         " simulation [default: %default]");
IMP_NPC_PARAMETER_STRING(init_rmffile, "",
                         "RMF file for initializing the simulation with its"
                         " last frame (to continue a previous run). If not"
                         " specified, initialization is through"
                         " pre-optimization");
IMP_NPC_PARAMETER_STRING(final_conformations, "final_conformations.rmf",
                 "RMF file for recording the initial and final conformations "
                         " [default: %default]");
#ifdef IMP_NPC_GOOGLE
IMP_NPC_PARAMETER_BOOL(profile, false,
                       "Whether to turn on profiling for the first run");
#endif
IMP_NPC_PARAMETER_BOOL(debug_initialization, false,
                       "Print more info about initialization");
IMP_NPC_PARAMETER_BOOL(initialize_only, false,
                       "Run the initialization and then stop");
IMP_NPC_PARAMETER_BOOL(quick, false,
                       "Reduce all steps to the minimum");
IMP_NPC_PARAMETER_BOOL(show_steps, false,
                       "Show the steps for each modified variable");
IMP_NPC_PARAMETER_BOOL(show_number_of_work_units, false,
                       "Show the number of work units");
IMP_NPC_PARAMETER_INT(seed_offset, 0,
                      "a number to add to the random seed of IMP random numbers"
                      " generator" );
IMP_NPC_PARAMETER_BOOL(const_seed, 0,
                      "a number to be used as a constant seed"
                       " (instead of system time)"
                       " in the IMP random numbers generator."
                       " If unspecified or zero, use system time" );


#ifdef IMP_BENCHMARK_USE_GOOGLE_PERFTOOLS_PROFILE
#define IMP_NPC_SET_PROF(p, tf) if (FLAGS_profile && i==0) {          \
  p.set("profiling.pprof");                                           \
  }
#else
#define IMP_NPC_SET_PROF(p, tf)
#endif

/** Run simulation using preconstructed SimulationData (sim_data) object.
    init_restraints are used ad-hoc during initialization only,
    and unless initialized from an RMF file
*/
#define IMP_NPC_LOOP(sim_data, init_restraints)                        \
  IMP::npctransport::internal::do_main_loop(sim_data,                  \
                                            init_restraints,           \
                                            FLAGS_quick,               \
                                            FLAGS_initialize_only,     \
                                            FLAGS_final_conformations,\
                                            FLAGS_debug_initialization, \
                                            FLAGS_init_rmffile)

//! seeds the random number generator of IMP with const_seed
//! (or time if it is zero), offsetted by seed_offset
void seed_randn_generator(IntArg const_seed, IntArg seed_offset)
{
  IntArg seed = const_seed;
  if(seed == 0)
    seed = static_cast<IntArg> (std::time(IMP::nullptr));
  IMP::base::random_number_generator.seed( seed + seed_offset );
}

inline IMP::npctransport::SimulationData *startup(int argc, char *argv[]) {
  IMP_NPC_PARSE_OPTIONS(argc, argv);
  IMP_NPC_PRINTHELP;
  seed_randn_generator(FLAGS_const_seed, FLAGS_seed_offset);
  set_log_level(IMP::base::LogLevel(FLAGS_log_level));
  try {
    int num=IMP::npctransport::assign_ranges
        (FLAGS_configuration, FLAGS_output,
         FLAGS_work_unit, FLAGS_show_steps);
    if (FLAGS_show_number_of_work_units) {
#pragma omp critical
      std::cout << "work units " << num << std::endl;
    }
    set_log_level(IMP::base::LogLevel(FLAGS_log_level));
  } catch (const IMP::base::IOException &e) {
#pragma omp critical
   std::cerr << "Error: " << e.what() << std::endl;
    return IMP_NULLPTR;
  }
  IMP_NEW(IMP::npctransport::SimulationData, sd,(FLAGS_output,
                              FLAGS_quick));
  try {
    if (!FLAGS_conformations.empty()) {
      sd->set_rmf_file_name(FLAGS_conformations);
    }
  } catch(const RMF::Exception &e) {
#pragma omp critical
    std::cerr << "Error: " << e.what() << std::endl;
    return IMP_NULLPTR;
  }
  return sd.release();
}

#endif // IMP_NPC_MAIN

#endif /* IMPNPCTRANSPORT_MAIN_H */
