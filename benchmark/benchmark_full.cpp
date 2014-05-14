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
#include <IMP/base/enums.h>
#include <IMP/base/flags.h>
#include <IMP/base/Pointer.h>
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
    IMP::base::Pointer<IMP::npctransport::SimulationData> sd;
    IMP_NPC_PARSE_OPTIONS(argc, argv);
    std::string input =
        IMP::npctransport::get_data_path("benchmark_full_restart.pb");
    std::string output = IMP::base::create_temporary_file_name("output", ".pb");
    {
      IMP_OMP_PRAGMA(critical);
      ::npctransport_proto::Output prev_output;
      std::ifstream file(input.c_str(), std::ios::binary);
      bool read = prev_output.ParseFromIstream(&file);
      IMP_ALWAYS_CHECK(read, "Couldn't read restart file " << input,
                       IMP::base::ValueException);
      prev_output.mutable_assignment()
          ->set_random_seed(IMP::base::get_random_seed());
      std::ofstream outf(output.c_str(), std::ios::binary);
      prev_output.SerializeToOstream(&outf);
    }
    IMP::base::set_log_level(IMP::base::SILENT);
    sd = new IMP::npctransport::SimulationData(output, false);

    sd->get_model()->update();
    double timev, score = 0;
    if(IMP::base::get_check_level() >= IMP::base::USAGE ) {
       //       || IMP::base::get_is_quick_test) {
      IMP_TIME(score += sd->get_bd()->optimize(1), timev);
    } else {
      IMP_TIME(score += sd->get_bd()->optimize(1000), timev);
    }
#ifdef _OMP
    std::string algo = "openmp";
#else
    std::string algo = "serial";
#endif
    IMP::benchmark::report("full npc", algo, timev, score);
  }
  catch (const IMP::base::Exception & e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
