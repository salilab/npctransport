/**
# * \file fg_simulation.cpp
# * \brief Simulate an fg and a kap interacting
#
# * Copyright 2007-2012 IMP Inventors. All rights reserved.
# */

#include <IMP/npctransport/main.h>
#include <IMP/npctransport/npctransport_config.h>
#include <IMP/npctransport/particle_types.h>
#include <IMP/npctransport/util.h>
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
#include <iostream>

#ifdef IMP_NPC_GOOGLE
#include <IMP/npctransport/internal/google_main.h>
#else
#include <IMP/npctransport/internal/boost_main.h>
#endif



int main(int argc, char *argv[]) {
#ifdef IMP_NPC_GOOGLE
  std::string config_txt =
    "third_party/npc/npctransport/data/benchmark_initialize.txt";
#else
  std::string config_txt =
    IMP::npctransport::get_data_path("benchmark_initialize.txt");
#endif
  IMP::base::AddStringFlag add_input("input", "Input file", &config_txt);
  // preparation::
  try {
     IMP_NPC_PARSE_OPTIONS(argc, argv);
     std::string config_pb = IMP::base::create_temporary_file_name
       ("benchmark_initalize.pb", ".pb");
     IMP::npctransport::configuration_txt2pb(config_txt, config_pb);
     std::string output = IMP::base::create_temporary_file_name("output", ".pb");

    int num = IMP::npctransport::assign_ranges(config_pb, output, 100, false,
                                               IMP::base::get_random_seed());

    IMP::base::Pointer<IMP::npctransport::SimulationData> sd =
        new IMP::npctransport::SimulationData(output, true);

    bool verbose = false;
    unsigned int acceleration_factor = 100;
    double short_init_factor = 1.0 / acceleration_factor;
    IMP::npctransport::initialize_positions(sd, IMP::RestraintsTemp(),
                                            verbose, short_init_factor);
    std::cout << "Energy after benchmark initialization "
              << "(acclerated by x" << acceleration_factor << "): "
              << sd->get_bd()->get_scoring_function()->evaluate(false)
              << std::endl;
  }
  catch (const IMP::base::Exception & e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return -1;
  }
  return 0;
}
