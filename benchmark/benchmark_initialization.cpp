/**
# * \file fg_simulation.cpp
# * \brief Simulate an fg and a kap interacting
#
# * Copyright 2007-2012 IMP Inventors. All rights reserved.
# */

#include <IMP/npctransport/main.h>
#include <IMP/npctransport/npctransport_config.h>
#include <IMP/npctransport/particle_types.h>
#include <IMP/algebra.h>
#include <IMP/core.h>
#include <IMP/atom.h>
#include <IMP/display.h>
#include <IMP/rmf.h>
#include <IMP.h>
#include <RMF/utility.h>
#include <IMP/container/SingletonsRestraint.h>
#include <IMP/base/CreateLogContext.h>
#include <IMP/base/exception.h>
#include <IMP/npctransport.h>
#include <IMP/ParticleTuple.h>
#include <IMP/base_types.h>
#include <IMP/base/Pointer.h>
#include <IMP/Restraint.h>
#include <IMP/SingletonScore.h>
#include <IMP/core/Typed.h>
#include <numeric>
#include <cmath>
#include <iostream>
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



int main(int argc, char *argv[])
{
#ifdef IMP_NPC_GOOGLE
    std::string input = "third_party/npc/npctransport/data/benchmark_initialize.pb";
#else
    std::string input = IMP::npctransport::get_data_path("benchmark_initialize.pb");
#endif
    IMP::base::AddStringFlag add_input("input", "Input file", &input);
  // preparation::
  try {
    IMP_NPC_PARSE_OPTIONS(argc, argv);
    std::string output = IMP::base::create_temporary_file_name("output", ".pb");


    int num=IMP::npctransport::assign_ranges
        (input, output, 100, false, IMP::base::get_random_seed() );

    IMP::base::Pointer<IMP::npctransport::SimulationData> sd
        = new IMP::npctransport::SimulationData(output, true);

    IMP::npctransport::initialize_positions(sd, IMP::RestraintsTemp(), false);
  }
  catch (const IMP::base::Exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return -1;
  }
  return 0;
 }
