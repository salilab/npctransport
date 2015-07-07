/**
 * \file fg_kap.cpp
 * \brief Simulate an fg and a kap interacting
 * Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#define IMP_NPC_MAIN
#include <IMP/npctransport/main.h>
#include <RMF/utility.h>

int main(int argc, char *argv[]) {
  using namespace IMP::npctransport;
  IMP::Pointer<SimulationData> sd = startup(argc, argv);
  do_main_loop(sd, IMP::RestraintsTemp());
  return 0;
}
