/**
# * \file fg_simulation.cpp
# * \brief Simulate an fg and a kap interacting
#
# * Copyright 2007-2012 IMP Inventors. All rights reserved.
# */

#include <IMP/npctransport/util.h>
#include <IMP/npctransport/typedefs.h>
#include <IMP/base/check_macros.h>
#include <IMP/base/exception.h>
#include <IMP/core/XYZR.h>
#include <IMP/core/Typed.h>
#include <IMP/base_types.h>
#include <cmath>
#include <iostream>

#ifdef IMP_NPC_GOOGLE
IMP_GCC_PUSH_POP(diagnostic push)
IMP_GCC_PRAGMA(diagnostic ignored "-Wsign-compare")
#include "third_party/npc/npctransport/data/npctransport.pb.h"
IMP_GCC_PUSH_POP(diagnostic pop)
#else
#include <IMP/npctransport/internal/npctransport.pb.h>
#endif

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

ParticlesTemp get_optimizable_particles
(ParticlesTemp const& particles)
{
  ParticlesTemp optimizables;
  for(unsigned int i = 0; i < particles.size(); i++)
    {
      Particle* p = particles[i];
      if(core::XYZ::get_is_setup(p)) {
        core::XYZ p_xyz(p);
        if(p_xyz.get_coordinates_are_optimized()){
          optimizables.push_back ( p );
        }
      }
    }
  return optimizables;
}


// finds the index of s->fgs() whose type equals pt.get_string()
// if it does not exist, add it to s
// @param s the statistics message for searching the fg
// @param pt the type of fg to look for
unsigned int find_or_add_fg_of_type(::npctransport_proto::Statistics* s,
                                    IMP::core::ParticleType pt)
{
  std::string type = pt.get_string();
  // find fg with same type
  unsigned int n = s->fgs_size();
  for(unsigned int i=0 ; i < n ; i++){
    if(s->fgs(i).type() == type){
      return i;
    }
  }
  // add new one if not found
  ::npctransport_proto::Statistics_FGStats* sfgs = s->add_fgs();
  sfgs->set_type( type );
  IMP_USAGE_CHECK(s->fgs(n).type() == type,
                  "bug - new floaters type should have been" << type);
  return n;
}

// finds the index of s->floaters() whose type equals pt.get_string()
// if it does not exist, add it to s
// @param s the statistics message for searching t
// @param pt the type to look for
unsigned int find_or_add_floater_of_type(::npctransport_proto::Statistics* s,
                                         IMP::core::ParticleType pt)
{
  std::string type = pt.get_string();
  // find fg with same type
  unsigned int n = s->floaters_size();
  for(unsigned int i=0 ; i < n ; i++){
    if(s->floaters(i).type() == type){
      return i;
    }
  }
  // add new one if not found
  ::npctransport_proto::Statistics_FloatStats* sfs = s->add_floaters();
  sfs->set_type( type );
  IMP_USAGE_CHECK(s->floaters(n).type() == type,
                  "bug - new type should have been" << type);
  return n;
}

// finds the index of s->interactions() whose type0 and type1 particle
// type equivalents are equale to it. If it does not exist, add it to s
// @param s the statistics message for searching t
// @param it the type to look for
unsigned int find_or_add_interaction_of_type
( ::npctransport_proto::Statistics* s, IMP::npctransport::InteractionType it)
{
  using namespace IMP;
  unsigned int n = s->interactions_size();
  for(unsigned int i=0 ; i < n ; i++){
    core::ParticleType pt0( s->interactions(i).type0() );
    core::ParticleType pt1( s->interactions(i).type1() );
    if(npctransport::InteractionType(pt0, pt1) == it) {
      return i;
    }
  }
  // add new one if not found
  ::npctransport_proto::Statistics_InteractionStats*
      sis = s->add_interactions();
  sis->set_type0( it.first.get_string() );
  sis->set_type1( it.second.get_string() );
  IMP_IF_CHECK(USAGE) {
    core::ParticleType pt0( s->interactions(n).type0() );
    core::ParticleType pt1( s->interactions(n).type1() );
    IMP_USAGE_CHECK( npctransport::InteractionType(pt0, pt1) == it,
                     "bug - new type should have been" << it );
  }
  return n;
}


IMPNPCTRANSPORT_END_NAMESPACE
