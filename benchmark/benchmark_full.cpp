/**
# * \file fg_simulation.cpp
# * \brief Simulate an fg and a kap interacting
#
# * Copyright 2007-2019 IMP Inventors. All rights reserved.
# */

#include <IMP/npctransport/main.h>
#include <IMP/npctransport/npctransport_config.h>
#include <IMP/algebra.h>
#include <IMP/core.h>
#include <IMP/rmf.h>
#include <IMP.h>
#include <RMF/utility.h>
#include <IMP/container/SingletonsRestraint.h>
#include <IMP/CreateLogContext.h>
#include <IMP/exception.h>
#include <IMP/npctransport.h>
#include <IMP/base_types.h>
#include <IMP/enums.h>
#include <IMP/flags.h>
#include <IMP/Pointer.h>
#include <IMP/Restraint.h>
#include <IMP/SingletonScore.h>
#include <IMP/benchmark/utility.h>
#include <IMP/benchmark/benchmark_macros.h>
#include <IMP/core/Typed.h>
#include <numeric>
#include <cmath>
#include <iostream>
#include <IMP/npctransport/protobuf.h>
#include <IMP/npctransport/internal/boost_main.h>
#include <IMP/npctransport/internal/npctransport.pb.h>
#include <IMP/npctransport/internal/boost_main.h>


int main(int argc, char *argv[]) {
  try {
    IMP::Pointer<IMP::npctransport::SimulationData> sd;
    IMP_NPC_PARSE_OPTIONS(argc, argv);
    std::string input =
        IMP::npctransport::get_data_path("benchmark_full_restart.pb");
    std::string output = IMP::create_temporary_file_name("output", ".pb");
    {
      IMP_OMP_PRAGMA(critical);
      ::npctransport_proto::Output prev_output;
      std::ifstream file(input.c_str(), std::ios::binary);
      bool read = prev_output.ParseFromIstream(&file);
      IMP_ALWAYS_CHECK(read, "Couldn't read restart file " << input,
                       IMP::ValueException);
      prev_output.mutable_assignment()
          ->set_random_seed(IMP::get_random_seed());
      std::ofstream outf(output.c_str(), std::ios::binary);
      prev_output.SerializeToOstream(&outf);
    }
    IMP::set_log_level(IMP::SILENT);
    sd = new IMP::npctransport::SimulationData(output, false);
    sd->activate_statistics();

    sd->get_model()->update();
    double timev, score = 0;
    if(IMP::get_check_level() >= IMP::USAGE ) {
       //       || IMP::get_is_quick_test) {
      IMP_TIME(score += sd->get_bd()->optimize(1), timev);
    } else {
      IMP_TIME(score += sd->get_bd()->optimize(10000), timev);
    }
#ifdef _OMP
    std::string algo = "openmp";
#else
    std::string algo = "serial";
#endif
    IMP::benchmark::report("full npc", algo, timev, score);
  }
  catch (const IMP::Exception & e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
