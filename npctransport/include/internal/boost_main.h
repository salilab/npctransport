/**
 *  \file boost_map.h
 *  \brief Macros for using npc executable with boost program options
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */
#ifndef IMPNPCTRANSPORT_BOOST_MAIN_H
#define IMPNPCTRANSPORT_BOOST_MAIN_H

#ifdef IMP_NPC_MAIN
#include <IMP/base/flags.h>

#define IMP_NPC_PARAMETER_BOOL(name, def, description)                  \
  bool FLAGS_##name=def;                                                \
  IMP::base::AddBoolFlag name##adder(#name, description, &FLAGS_##name)

#define IMP_NPC_PARAMETER_INT64(name, def, description)                 \
  boost::int64_t FLAGS_##name=def;                                      \
  IMP::base::AddIntFlag name##adder(#name, description, &FLAGS_##name)

#define IMP_NPC_PARAMETER_UINT64(name, def, description)                \
  boost::int64_t FLAGS_##name=def;                                      \
  IMP::base::AddIntFlag name##adder(#name, description, &FLAGS_##name)


#define IMP_NPC_PARAMETER_STRING(name, def, description)                \
  std::string FLAGS_##name=def;                                         \
  IMP::base::AddStringFlag name##adder(#name, description, &FLAGS_##name)

#define IMP_NPC_PARSE_OPTIONS(argc, argv)                               \
  IMP::base::setup_from_argv(argc, argv, "Something npcish");

#endif




#endif /* IMPNPCTRANSPORT_BOOST_MAIN_H */
