/**
 *  \file protobuf.h
 *  \brief A Score on the distance between a pair of particles.
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_PROTOBUF_H
#define IMPNPCTRANSPORT_PROTOBUF_H

#include "npctransport_config.h"
#include <google/protobuf/repeated_field.h>
#include <boost/cstdint.hpp>
#include <set>

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

/** returns an STL set with all unique (no duplicates) values in
    ::google::protobuf::RepeatedField object

    @param[in] src repeated field from which values of type value_type
                   are copied
*/
template<class value_type>
std::set< value_type >
get_unique_set_from_repeated_field
( ::google::protobuf::RepeatedField< value_type > const& src )
{
  std::set< value_type > ret_set;
  for(int i = 0; i < src.size(); i++){
    ret_set.insert( src.Get(i) );
  }
  return ret_set;
}

/** returns ::google::protobuf::RepeatedField object
    based on an STL set

    @param[in] src set from which values of type value_type
                   are copied
    @param[out] trg repeated field to which values are copied
*/
/* template<class value_type> */
/* void copy_unique_set_to_repeated_field */
/* ( std::set< value_type > const& src, */
/*   ::google::protobuf::RepeatedField< value_type >& trg ) */
/* { */
/*   ::google::protobuf::RepeatedField< value_type > ret_val; */
/*   std::set< value_type >::const_iterator;// iter; */
/*   for(iter = src.begin(); iter != src.end(); iter++) { */
/*     ret_val.push_back( *iter ); */
/*   } */
/*   return ret_val; */
/* } */

IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_PROTOBUF_H */
