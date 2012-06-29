/**
 *  \file protobuf.h
 *  \brief A Score on the distance between a pair of particles.
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_PROTOBUF_H
#define IMPNPCTRANSPORT_PROTOBUF_H

#include "npctransport_config.h"

namespace npctransport {
class OutputConfiguration;
}

IMPNPCTRANSPORT_BEGIN_NAMESPACE

IMPNPCTRANSPORTEXPORT void
show_ranges(std::string fname);

IMPNPCTRANSPORTEXPORT int
assign_ranges(std::string fname, std::string output, unsigned int work_unit,
              bool show_steps);

IMPNPCTRANSPORTEXPORT int
get_number_of_work_units(std::string configuration_file);

IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_PROTOBUF_H */
