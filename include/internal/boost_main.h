/**
 *  \file internal/boost_main.h
 *  \brief Macros for using npc executable with boost program options
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 */
#ifndef IMPNPCTRANSPORT_BOOST_MAIN_H
#define IMPNPCTRANSPORT_BOOST_MAIN_H

#include <IMP/flags.h>

#define IMP_NPC_PARSE_OPTIONS(argc, argv) \
  IMP::setup_from_argv(argc, argv, "Something npcish");

#endif /* IMPNPCTRANSPORT_BOOST_MAIN_H */
