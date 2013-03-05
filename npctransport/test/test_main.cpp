/**
 * \file test_main.cpp
 * \brief Test the macros and loading and all that for mains
 * Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#define IMP_NPC_MAIN
#include <IMP/npctransport/main.h>
#include <IMP/base/random.h>

int main(int , char *[]) {
  // TODO: emulate real runtime parameters by initializig fake argc / argv?
  try {
    std::string config
        = IMP::npctransport::get_data_path("quick.pb");
    std::string assignment
        = IMP::base::create_temporary_file_name("output", ".pb");
    std::string output
        = IMP::base::create_temporary_file_name("output", ".rmf");
    set_log_level(IMP::base::LogLevel(IMP::base::SILENT));
    FLAGS_quick=true;
    int num=IMP::npctransport::assign_ranges
      (config, assignment, 0, true, IMP::base::get_random_seed());
    std::cout << "num ranges " << num << std::endl;
    IMP_NEW(IMP::npctransport::SimulationData, sd,(assignment, true));
    sd->set_rmf_file_name(output);
    sd->get_m()->set_log_level(IMP::base::SILENT);
    std::cout << "Files are " << assignment
              << " and " << output
              << std::endl;
    IMP_NPC_LOOP(sd, IMP::RestraintsTemp());
  } catch (IMP::Exception e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
