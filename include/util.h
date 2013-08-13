/**
# * \file util.h
# * \brief utility methods for npctransport

Simulate an fg and a kap interacting
#
# * Copyright 2007-2012 IMP Inventors. All rights reserved.
# */

#ifndef IMPNPCTRANSPORT_UTIL_H
#define IMPNPCTRANSPORT_UTIL_H

#include "npctransport_config.h"
#include <IMP/npctransport/npctransport_proto.fwd.h>
#include <string>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/**
   Converts protobuf configuration file config_txt (which is in pretty
   protobuf textual output format) to binary protobuf format
   (in file config_pb)

   @param config_txt the input textual protobuf config file
   @param config_pb the output binary protobuf file
 */
void configuration_txt2pb
(std::string config_txt, std::string config_pb);

/** returns particles with optimizable coordinates from particles */
ParticlesTemp get_optimizable_particles
(ParticlesTemp const& particles);


#ifndef SWIG

/** finds the index of s.fgs() whose type equals pt.get_string()
    if it does not exist, add it to s
    @param s the statistics message for searching the fg
    @param pt the type of fg to look for
*/
unsigned int find_or_add_fg_of_type(::npctransport_proto::Statistics* s,
                                    IMP::core::ParticleType pt);

/** finds the index of s->floaters() whose type equals pt.get_string()
    if it does not exist, add it to s
    @param s the statistics message for searching t
    @param pt the type to look for
*/
unsigned int find_or_add_floater_of_type(::npctransport_proto::Statistics* s,
                                         IMP::core::ParticleType pt);


/** finds the index of s.interactions() whose type0 and type1 particle
  type equivalents are equale to it. If it does not exist, add it to s
  @param s the statistics message for searching t
  @param it the type to look for
  @return the index of the interaction type in s.interaction()
*/
unsigned int find_or_add_interaction_of_type
( ::npctransport_proto::Statistics* s,
  IMP::npctransport::InteractionType it);

#endif  // SWIG

IMPNPCTRANSPORT_END_NAMESPACE


#endif /* IMPNPCTRANSPORT_UTIL_H */
