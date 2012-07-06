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
#include <IMP/benchmark/Profiler.h>
#include "internal/main.h"
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


#define IMP_NPC_STARTUP                                                 \
  IMP_NPC_START_INT;                                                    \
  IMP_NPC_PRINTHELP;                                                    \
  int num=assign_ranges(FLAGS_configuration, FLAGS_assignments,         \
                FLAGS_work_unit, FLAGS_show_steps);                     \
  if (FLAGS_show_number_of_work_units) {                                \
    std::cout << "work units " << num << std::endl;                     \
  }                                                                     \
  set_log_level(LogLevel(FLAGS_log_level));                             \
  IMP_NEW(SimulationData, sd,(FLAGS_assignments, FLAGS_statistics,      \
                              FLAGS_quick));                            \
  if (!FLAGS_conformations.empty()) {                                   \
    sd->set_rmf_file_name(FLAGS_conformations);                         \
  }



IMP_NPC_PARAMETER_INT(work_unit, -1, "The work unit");
IMP_NPC_PARAMETER_INT(log_level, 0, "The log level to use");
IMP_NPC_PARAMETER_STRING(configuration, "configuration.pb",
                         "Configuration file");
IMP_NPC_PARAMETER_STRING(assignments, "assignments.pb", "Assignments file");
IMP_NPC_PARAMETER_STRING(statistics, "statistics.pb", "Statistics file");
IMP_NPC_PARAMETER_STRING(final_configuration, "final.pym",
                         "Where to write the final config");
IMP_NPC_PARAMETER_STRING(configurations, "conformations.rmf",
                         "Where to write the conformations.");
IMP_NPC_PARAMETER_BOOL(profile, false,
                         "Whether to turn on profiling for the first run");
IMP_NPC_PARAMETER_BOOL(quick, false, "Reduce all steps to the minimum");
IMP_NPC_PARAMETER_BOOL(show_steps, false,
                       "Show the steps for each modified variable");
IMP_NPC_PARAMETER_BOOL(show_number_of_work_units, false,
                       "Show the number of work units");

#ifdef IMP_BENCHMARK_USE_GOOGLE_PERFTOOLS_PROFILE
#define IMP_NPC_SET_PROF(p, tf) if (FLAGS_profile && i==0) {             \
  p.set("profiling.pprof");                                           \
  }
#else
#define IMP_NPC_SET_PROF(p, tf)
#endif

#define IMP_NPC_LOOP(links)                                     \
  IMP::npctransport::internal::do_main_loop(sd, links, FLAGS_quick,     \
                                            FLAGS_final_configuration)

using namespace IMP;
using namespace IMP::npctransport;
using namespace IMP::atom;
using namespace IMP::core;
using namespace IMP::container;
using namespace IMP::display;
using namespace IMP::algebra;
using namespace RMF;

#endif // IMP_NPC_MAIN

#endif /* IMPNPCTRANSPORT_MAIN_H */
