/**
 * \file test_main.cpp
 * \brief Test the macros and loading and all that for mains
 * Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#define IMP_NPC_MAIN
#include <IMP/npctransport/main.h>
#include <RMF/utility.h>

int main(int , char *[]) {
  try {
    RMF::set_show_hdf5_errors(true);
    std::string config
        = IMP::npctransport::get_data_path("quick.pb");
    std::string assignment
        = IMP::base::create_temporary_file_name("assignment", ".pb");
    std::string statistics
        = IMP::base::create_temporary_file_name("statistics", ".pb");
    std::string output
        = IMP::base::create_temporary_file_name("output", ".rmf");
    set_log_level(LogLevel(IMP::base::VERBOSE));
    int num=assign_ranges(config, assignment,
                          0, true);
    std::cout << "num ranges " << num << std::endl;
    IMP_NEW(SimulationData, sd,(assignment,statistics,
                                true));
    sd->set_rmf_file_name(output);
    sd->get_m()->set_log_level(SILENT);
    std::cout << "Files are " << assignment << " and " << statistics
              << " and " << output
              << std::endl;
    IMP_NPC_LOOP(ParticlePairsTemp());
  } catch (IMP::Exception e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
