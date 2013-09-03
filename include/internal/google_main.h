/**
 *  \file google_map.h
 *  \brief Macros for using npc executable with boost program options
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */
#ifndef IMPNPCTRANSPORT_GOOGLE_MAIN_H
#define IMPNPCTRANSPORT_GOOGLE_MAIN_H

#include "base/init_google.h"
#include "base/commandlineflags.h"
#include <IMP/base/flags.h>
#include <boost/scoped_array.hpp>

inline void parse_argv_argc(int argc, char** argv) {
  IMP::Strings gargs =
      IMP::base::setup_from_argv_allowing_unknown(argc, argv, "something npc");
  gargs.insert(gargs.begin(), argv[0]);
  int gargc = gargs.size();
  boost::scoped_array<char*> gargv(new char* [gargs.size()]);
  for (unsigned int i = 0; i < gargs.size(); ++i) {
    gargv[i] = const_cast<char*>(gargs[i].c_str());
  }
  char** gargvp = gargv.get();
  InitGoogle("Something npc related", &gargc, &gargvp, true);
}

#define IMP_NPC_PARSE_OPTIONS(argc, argv) parse_argv_argc(argc, argv);

#endif /* IMPNPCTRANSPORT_GOOGLE_MAIN_H */
