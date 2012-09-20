/**
 * \file fg_kap.cpp
 * \brief Simulate an fg and a kap interacting
 * Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#define IMP_NPC_MAIN
#include <IMP/npctransport/main.h>
#include <IMP/ParticleTuple.h>
#include <RMF/utility.h>

int main(int argc, char *argv[]) {
  RMF::set_show_hdf5_errors(true);
  IMP_NPC_STARTUP(sd);
  IMP_NPC_LOOP(sd, IMP::ParticlePairsTemp());
  return 0;
}
