/**
 *  \file google_map.h
 *  \brief Macros for using npc executable with boost program options
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */
#ifndef IMPNPCTRANSPORT_GOOGLE_MAIN_H
#define IMPNPCTRANSPORT_GOOGLE_MAIN_H

#ifdef IMP_NPC_MAIN
#include "base/init_google.h"
#include "base/commandlineflags.h"

typedef int Int64Arg;
typedef unsigned int UInt64Arg;

#define IMP_NPC_PARAMETER(name, def, description)        \
  DEFINE_double(name, def, description);             \
  Pusher name##_pusher(#name, &FLAGS_##name)

#define IMP_NPC_PARAMETER_BOOL(name, def, description)        \
  DEFINE_bool(name, def, description)

#define IMP_NPC_PARAMETER_INT64(name, def, description)        \
  DEFINE_int64(name, def, description)

#define IMP_NPC_PARAMETER_UINT64(name, def, description)        \
  DEFINE_int64(name, def, description)

#define IMP_NPC_PARAMETER_STRING(name, def, description) \
  DEFINE_string(name, def, description)

#define IMP_NPC_PARSE_OPTIONS(argc, argv)       \
  InitGoogle(argv[0],                           \
             &argc, &argv, true);

#endif
#endif /* IMPNPCTRANSPORT_GOOGLE_MAIN_H */
