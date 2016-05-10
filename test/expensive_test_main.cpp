/**
 * \file test_main.cpp
 * \brief Test the macros and loading and all that for mains
 * Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#define IMP_NPC_MAIN
#include <IMP/npctransport/Scoring.h>
#include <IMP/npctransport.h>
#include <IMP/npctransport/main.h>
#include <IMP/exception.h>
#include <IMP/flags.h>
#include <IMP/types.h>
#include <vector>
#include <string>

int main(int argc, char *argv[]) {
  // TODO: emulate real runtime parameters by initializig fake argc / argv?
  try {
    // add a quick init flag (to make test quicker)
    IMP::Strings new_argv;
    for(int i=0; i < argc; i++){
      new_argv.push_back( argv[i] );
    }
    new_argv.push_back("--short_init_factor");
    new_argv.push_back("0.1");
    IMP::setup_from_argv(new_argv,
                               "Test of the main loop for npc transport");
    // prepare assignment and simulation data
    std::string config = IMP::npctransport::get_data_path("quick.pb");
    std::string assignment =
        IMP::create_temporary_file_name("output", ".pb");
    std::string output =
        IMP::create_temporary_file_name("output", ".rmf");
    IMP::set_log_level(IMP::LogLevel(IMP::PROGRESS));
    int num = IMP::npctransport::assign_ranges(config, assignment, 0, true,
                                               IMP::get_random_seed());
    std::cout << "num ranges " << num << std::endl;
    IMP_NEW(IMP::npctransport::SimulationData, sd,
            (assignment, true /* quick */));
    sd->activate_statistics();
    sd->set_rmf_file(output);
    sd->get_model()->set_log_level(IMP::PROGRESS);
    std::cout << "Files are " << assignment << " and " << output << std::endl;
    // simulate
    IMP::npctransport::do_main_loop(sd, IMP::RestraintsTemp());
  }
  catch (IMP::Exception e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
