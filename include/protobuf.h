/**
 *  \file protobuf.h
 *  \brief A Score on the distance between a pair of particles.
 *
 *  Copyright 2007-2018 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_PROTOBUF_H
#define IMPNPCTRANSPORT_PROTOBUF_H

#include "npctransport_config.h"
#include <IMP/SingletonContainer.h>
#include <IMP/core/Typed.h>
#include <boost/cstdint.hpp>
#include <set>

#ifndef SWIG
// instead of including protobuf header, which is problematic due to
// minor issue with google headers namespaces
namespace npctransport_proto {
  class Conformation;
  class Output;
}
#endif

IMPNPCTRANSPORT_BEGIN_NAMESPACE


IMPNPCTRANSPORTEXPORT void show_ranges(std::string fname);

/**
   reads the protobuf message in [input_config_fname], which may contain submessages with
   ranges of values (indicated by presence of .upper and .lower fields, with
   .steps possible steps for each such field). The output is a message with
   the [work_unit]'th possible combination of these ranges, to the file
   [output_assignment_fname].

   Note: the range values are enumerated as if they lie on a grid with
   log-evenly distributed axis-aligned grid points, using the .base field
   as the log base for each ranged field, such that e.g. iterating over the
   range [1..8] with 3 steps and base 2 will be enumerated as (1,4,8)

   @param fname input_config_fname configuration file name
   @param output output_assignment_fname assignment file name
   @param work_unit the index of combination of range values to be used. If the
                    total of possible combinations of all fields with ranges is
   k,
                    it is guaranteed that iterating over work_unit between
   0..k-1
                    will enumerate over all possible combinations, and that
                    work_unit and (work_unit % k) will return the same output
   for
                    the same input.
   @param show_steps show the steps that occur
   @param random_seed the random seed used to initialize the IMP random number
                      generator for this simulation

   @throw IMP::ValueException if any of the values in the configuration
   file are in conflict (e.g., simulation time and maximal number of
   frames)
*/
// Each range field also has .steps and .base field.
IMPNPCTRANSPORTEXPORT int assign_ranges(
    std::string input_config_fname, std::string output_assignment_fname, unsigned int work_unit,
    bool show_steps,
    boost::uint64_t random_seed  // do not use boost::uint64_t cause of SWIG
    );

IMPNPCTRANSPORTEXPORT int get_number_of_work_units(
    std::string configuration_file);

#ifndef SWIG
/**
   Loads a protobuf conformation into the diffusers and sites

   @param conformation the saved conformation protobuf message
   @param beads corresponding diffusers to be updated
   @param sites a map of sites for each diffuser particle type
                to be updated

   @note the beads and sites must have the same structure
                       as the ones used when saving (e.g. their
                       non-changing variables are expected to
                       be identical, and they differ only in the
                       dynamic ones)
   \see save_pb_conformation
*/
void load_pb_conformation
( const ::npctransport_proto::Conformation &conformation,
  IMP::SingletonContainerAdaptor beads,
  boost::unordered_map<core::ParticleType, algebra::Sphere3Ds> &sites);

/**
   Saves a protobuf conformation from the diffusers and sites

   @param beads beads to save
   @param sites a map of sites for each diffuser particle type
                to be saved
   @param conformation the conformation protobuf message to be save

   \see load_pb_conformation
 */
void save_pb_conformation
( IMP::SingletonContainerAdaptor beads,
  const boost::unordered_map<core::ParticleType, algebra::Sphere3Ds> &sites,
  ::npctransport_proto::Conformation *conformation );

//! load file output_fname into protobuf output object output
//! return true if succesful
bool load_output_protobuf(std::string output_fname,
                          ::npctransport_proto::Output& output);
#endif

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_PROTOBUF_H */
