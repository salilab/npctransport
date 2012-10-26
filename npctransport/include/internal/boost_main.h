/**
 *  \file boost_map.h
 *  \brief Macros for using npc executable with boost program options
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */
#ifndef IMPNPCTRANSPORT_BOOST_MAIN_H
#define IMPNPCTRANSPORT_BOOST_MAIN_H

#ifdef IMP_NPC_MAIN
#include <boost/program_options.hpp>

typedef int IntArg;
typedef unsigned int UIntArg;


#define IMP_NPC_PARAMETER(name, def, description)                       \
  double FLAGS_##name=def;                                              \
  Adder name##add(#name,                                                \
                                     &FLAGS_##name,                     \
                                     description);                      \
  Pusher name##_pusher(#name, &FLAGS_##name)


#define IMP_NPC_PARAMETER_BOOL(name, def, description)    \
  bool FLAGS_##name=def;                                  \
  Adder name##add(#name, &FLAGS_##name,                   \
                                         description)

#define IMP_NPC_PARAMETER_INT(name, def, description)                   \
  int FLAGS_##name=def;                                                 \
  Adder name##add(#name, &FLAGS_##name,                                 \
                  description);                                         \

#define IMP_NPC_PARAMETER_UINT(name, def, description)                   \
  unsigned int FLAGS_##name=def;                                                 \
  Adder name##add(#name, &FLAGS_##name,                                 \
                  description);                                         \

#define IMP_NPC_PARAMETER_STRING(name, def, description)                \
  std::string FLAGS_##name=def;                                         \
  Adder name##add(#name, &FLAGS_##name,                                 \
                  description);

#define IMP_NPC_PARSE_OPTIONS(argc, argv)                               \
  boost::program_options::variables_map vm;                             \
  boost::program_options::store(boost                                   \
                                ::program_options::parse_command_line(argc, \
                                                                      argv, \
                                                                      desc), \
                                vm);                                    \
  boost::program_options::notify(vm);                                   \
  if (vm.count("help")) {                                               \
    std::cout << desc << "\n";                                          \
    return IMP_NULLPTR;                                                 \
  }

#define IMP_NPC_PRINTHELP                       \
  if (FLAGS_help) {                             \
    std::cerr << desc;                          \
    return IMP_NULLPTR;                         \
  }

#endif
boost::program_options::options_description desc;
struct Adder {
  template <class V>
  Adder(std::string name, V*v, std::string description) {
    desc.add_options()(name.c_str(), boost::program_options::value<V>(v),
                       description.c_str());
  }
  Adder(std::string name, bool *b, std::string description) {
    desc.add_options()(name.c_str(),
                       boost::program_options::value<bool>(b)->zero_tokens(),
                       description.c_str());
  }
};

IMP_NPC_PARAMETER_BOOL(help, false, "Print help");




#endif /* IMPNPCTRANSPORT_BOOST_MAIN_H */
