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

IMPNPCTRANSPORT_END_NAMESPACE


#endif /* IMPNPCTRANSPORT_UTIL_H */
