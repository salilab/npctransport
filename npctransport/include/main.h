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
#include <IMP/npctransport.h>
//#include <IMP/benchmark/Profiler.h>
#include "internal/main.h"
#include "protobuf.h"
#include <numeric>
#include <cmath>
#include <iostream>

// use the example code for now to work bugs out of it
#include <IMP/example/creating_restraints.h>
#include <IMP/example/randomizing.h>
#include <IMP/example/counting.h>
#include <IMP/example/optimizing.h>

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
  IMP_NPC_START_INT;                                                    \
  IMP_NPC_PRINTHELP;                                                    \
  set_log_level(IMP::base::LogLevel(FLAGS_log_level));                  \
  try {                                                                 \
    int num=IMP::npctransport::assign_ranges                            \
        (FLAGS_configuration, FLAGS_assignments,                        \
         FLAGS_work_unit, FLAGS_show_steps);                            \
    if (FLAGS_show_number_of_work_units) {                              \
      std::cout << "work units " << num << std::endl;                   \
    }                                                                   \
    set_log_level(IMP::base::LogLevel(FLAGS_log_level));                \
  } catch (const IMP::base::IOException &e) {                           \
    std::cerr << "Error: " << e.what() << std::endl;                    \
    return 1;                                                           \
  }                                                                     \
  IMP_NEW(SimulationData, sim_data,(FLAGS_assignments, FLAGS_statistics,\
                              FLAGS_quick));                            \
  try {                                                                 \
    if (!FLAGS_conformations.empty()) {                                 \
      sd->set_rmf_file_name(FLAGS_conformations);                       \
    }                                                                   \
  } catch(const RMF::Exception &e) {                                    \
    std::cerr << "Error: " << e.what() << std::endl;                    \
    return 1;                                                           \
  }                                                                     \




IMP_NPC_PARAMETER_INT(work_unit, -1, "The work unit");
IMP_NPC_PARAMETER_INT(log_level, 0, "The log level to use");
IMP_NPC_PARAMETER_STRING(configuration, "configuration.pb",
                         "input configuration file in protobuf format"
                         " [default: %default]");
IMP_NPC_PARAMETER_STRING(assignments, "assignments.pb",
                         "output assignments file in protobuf format,"
                         " recording the assignment being executed"
                         " [default: %default]");
IMP_NPC_PARAMETER_STRING(statistics, "statistics.pb",
                         "output statistics file in protobuf format"
                         " [default: %default]");
IMP_NPC_PARAMETER_STRING(final_configuration, "final.pym",
                         "output final configuration file"
                         " [default: %default]");
IMP_NPC_PARAMETER_STRING(conformations, "conformations.rmf",
                         "RMF file for recording the conforomations along the "
                         " simulation [default: %default]");
IMP_NPC_PARAMETER_STRING(init_rmffile, "",
                         "RMF file for initializing the simulation with its"
                         " last frame (to continue a previous run). If not"
                         " specified, initialization is through"
                         " pre-optimization");

#ifdef IMP_NPC_GOOGLE
IMP_NPC_PARAMETER_BOOL(profile, false,
                       "Whether to turn on profiling for the first run");
#endif
IMP_NPC_PARAMETER_BOOL(debug_initialization, false,
                       "Print more info about initialization");
IMP_NPC_PARAMETER_BOOL(quick, false,
                       "Reduce all steps to the minimum");
IMP_NPC_PARAMETER_BOOL(show_steps, false,
                       "Show the steps for each modified variable");
IMP_NPC_PARAMETER_BOOL(show_number_of_work_units, false,
                       "Show the number of work units");

#ifdef IMP_BENCHMARK_USE_GOOGLE_PERFTOOLS_PROFILE
#define IMP_NPC_SET_PROF(p, tf) if (FLAGS_profile && i==0) {          \
  p.set("profiling.pprof");                                           \
  }
#else
#define IMP_NPC_SET_PROF(p, tf)
#endif

/** TODO: what is the meaning of links? */
#define IMP_NPC_LOOP(sim_data, links)                                   \
  IMP::npctransport::internal::do_main_loop(sim_data, links, FLAGS_quick,     \
                                            FLAGS_final_configuration,\
                                            FLAGS_debug_initialization, \
                                            FLAGS_init_rmffile)

using namespace IMP::npctransport;

#endif // IMP_NPC_MAIN

#endif /* IMPNPCTRANSPORT_MAIN_H */
