/**
# * \file fg_simulation.cpp
# * \brief Simulate an fg and a kap interacting
#
# * Copyright 2007-2012 IMP Inventors. All rights reserved.
# */

#include <IMP/npctransport/npctransport_config.h>
#include <IMP/npctransport/initialize_positions.h>
#include <IMP/npctransport/main.h>
#include <IMP/npctransport/protobuf.h>
#include <IMP/npctransport/util.h>

#include <IMP/benchmark/utility.h>
#include <IMP/benchmark/benchmark_macros.h>

//#include <IMP/ParticleTuple.h>
#include <IMP/base_types.h>
#include <IMP/base/CreateLogContext.h>
#include <IMP/base/exception.h>
#include <IMP/base/flags.h>
#include <IMP/base/Pointer.h>
#include <IMP/Restraint.h>
//#include <RMF/utility.h>

#include <numeric>
#include <iostream>

#include <IMP/npctransport/internal/boost_main.h>

std::string config_txt =
  IMP::npctransport::get_data_path("benchmark_initialize.txt");

IMP::base::AddStringFlag add_input("input",
                                   "Input text file that is "
                                   "converted to a config file", &config_txt);


int main(int argc, char *argv[]) {
  try {
    // preparation::
     IMP_NPC_PARSE_OPTIONS(argc, argv);
     bool verbose = false;
     unsigned int acceleration_factor = 200;
    double short_init_factor = 1.0 / acceleration_factor;
    std::string config_pb = IMP::base::create_temporary_file_name
      ("benchmark_initalize.pb", ".pb");
    std::string output =
      IMP::base::create_temporary_file_name("output", ".pb");
    /*     std::string rmf_file =
           IMP::base::create_temporary_file_name("init", ".rmf");
           std::cout << "Output: " << output << "  rmf-file: " << rmf_file
           << std::endl;
    */

    // assign and run::
    IMP::npctransport::configuration_txt2pb(config_txt, config_pb);
    IMP::npctransport::assign_ranges(config_pb, output, 100, false,
                                     IMP::base::get_random_seed());

    IMP::base::Pointer<IMP::npctransport::SimulationData> sd =
      new IMP::npctransport::SimulationData(output, true); //, rmf_file);

    if(IMP::base::get_check_level() >= IMP::base::USAGE ) {
       //      || IMP::base::get_is_quick_test) {
      std::cout << "skipping actual call to initialize_positions"
                << " when check level is larger than USAGE" << std::endl;
      return 0;
    }

    // Test and report:
    std::ostringstream algo_oss;
    algo_oss << "Init accelerate " << acceleration_factor;
#ifdef _OMP
    algo_oss << " (openmp)";
#else
    algo_oss << " (serial)";
#endif
    double timev = 0.0;
    IMP_TIME
      ( IMP::npctransport::initialize_positions(sd, IMP::RestraintsTemp(),
                                                verbose, short_init_factor),
        timev);
    double score = sd->get_bd()->get_scoring_function()->evaluate(false);
    std::cout << "Energy after benchmark initialization "
              << algo_oss.str() << " - " << score << std::endl;
    IMP::benchmark::report("init npc", algo_oss.str(), timev, score);
  }
  catch (const IMP::base::Exception & e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return -1;
  }
  return 0;
}
