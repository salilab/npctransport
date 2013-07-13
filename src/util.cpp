/**
# * \file fg_simulation.cpp
# * \brief Simulate an fg and a kap interacting
#
# * Copyright 2007-2012 IMP Inventors. All rights reserved.
# */

#include <IMP/npctransport/util.h>
#include <IMP/base/check_macros.h>
#include <IMP/base/exception.h>
#include <IMP/base_types.h>
#include <cmath>
#include <iostream>

#ifdef IMP_NPC_GOOGLE
IMP_GCC_PUSH_POP(diagnostic push)
IMP_GCC_PRAGMA(diagnostic ignored "-Wsign-compare")
#include "third_party/npc/npctransport/data/npctransport.pb.h"
IMP_GCC_PUSH_POP(diagnostic pop)
#include <IMP/npctransport/internal/google_main.h>
using namespace proto2;
#else
#include <IMP/npctransport/internal/boost_main.h>
#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <IMP/npctransport/internal/npctransport.pb.h>
using namespace ::google::protobuf;
#endif

IMPNPCTRANSPORT_BEGIN_NAMESPACE

// Converts protobuf configuration file from txt to pb format
void configuration_txt2pb
(std::string config_txt, std::string config_pb)
{
  npctransport_proto::Configuration config;
  std::ifstream ifs_txt(config_txt.c_str());
  IMP_ALWAYS_CHECK(ifs_txt, "File " << config_txt << " not found",
                   IMP::base::IOException);
  io::IstreamInputStream isis_txt(&ifs_txt);
  TextFormat::Parse(&isis_txt, &config);
  ifs_txt.close();
  std::ofstream ofs_pb(config_pb.c_str());
  config.SerializeToOstream(&ofs_pb);
  ofs_pb.close();
 }

IMPNPCTRANSPORT_END_NAMESPACE
