/**
 *  \file protobuf.h
 *  \brief A Score on the distance between a pair of particles.
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_PROTOBUF_H
#define IMPNPCTRANSPORT_PROTOBUF_H

#include "npctransport_config.h"
#include <boost/cstdint.hpp>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

class OutputConfiguration; // TODO: is this needed?

IMPNPCTRANSPORTEXPORT void
show_ranges(std::string fname);

/**
   reads the protobuf message in [fname], which may contain submessages with
   ranges of values (indicated by presence of .upper and .lower fields, with
   .steps possible steps for each such field). The output is a message with
   the [work_unit]th possible combination of these ranges, to the file [ofname].

   Note: the range values are enumerated as if they lie on a grid with
   log-evenly distributed axis-aligned grid points, using the .base field
   as the log base for each ranged field, such that e.g. iterating over the
   range [1..8] with 3 steps and base 2 will be enumerated as (1,4,8)

   @param fname input file name
   @param ofname output filename
   @param work_unit the index of combination of range values to be used. If the
                    total of possible combinations of all fields with ranges is k,
                    it is guaranteed that iterating over work_unit between 0..k-1
                    will enumerate over all possible combinations, and that
                    work_unit and (work_unit % k) will return the same output for
                    the same input.
   @param random_seed the random seed used to initialize the IMP random number
                      generator for this simulation

   @throw IMP::base::ValueException if any of the values in the configuration file
              are in conflict (e.g., simulation time and maximal number of frames)
*/
// Each range field also has .steps and .base field.
IMPNPCTRANSPORTEXPORT int
assign_ranges(std::string fname, std::string output, unsigned int work_unit,
              bool show_steps,
              boost::uint64_t random_seed // do not use boost::uint64_t cause of SWIG
              );

IMPNPCTRANSPORTEXPORT int
get_number_of_work_units(std::string configuration_file);

IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_PROTOBUF_H */
