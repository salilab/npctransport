/**
 *  \file boost_map.h
 *  \brief Macros for using npc executable with boost program options
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */
#ifndef IMPNPCTRANSPORT_BOOST_MAIN_H
#define IMPNPCTRANSPORT_BOOST_MAIN_H

#include <IMP/base/flags.h>

#define IMP_NPC_PARSE_OPTIONS(argc, argv) \
  IMP::base::setup_from_argv(argc, argv, "Something npcish");

#endif /* IMPNPCTRANSPORT_BOOST_MAIN_H */
