/**
 *  \file io.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/protobuf.h>
#include <IMP/npctransport/automatic_parameters.h>
#ifdef IMP_NPC_GOOGLE
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#include "third_party/npc/npctransport/data/npctransport.pb.h"
#pragma GCC diagnostic pop
#else
#include <IMP/npctransport/internal/npctransport.pb.h>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/message.h>
#endif
#include <IMP/utility.h>
#include <IMP/algebra/GridD.h>
#include <IMP/atom/estimates.h>
#include <IMP/algebra/grid_storages.h>
#include <IMP/base/CreateLogContext.h>
#include <IMP/base/SetLogState.h>
#include <boost/scoped_ptr.hpp>
#include <fstream>
#include <iostream>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
#ifdef IMP_NPC_GOOGLE
using namespace proto2;
#else
using namespace ::google::protobuf;
#endif
namespace {
void show_ranges(std::string name, const Message* message) {
  const Reflection* r(message->GetReflection());
  /*std::cout << "inspecting " << name << " " << message->GetTypeName()
    << std::endl;*/
  const Descriptor* d(message->GetDescriptor());
  const FieldDescriptor* lfd(d->FindFieldByName("lower"));
  if (lfd) {
    const FieldDescriptor* fd(d->FindFieldByName("upper"));
    if (r->HasField(*message, fd)) {
      std::cout << "found " << name << std::endl;
    }
  } else {
    for (int i = 0; i < d->field_count(); ++i) {
      const FieldDescriptor* fd(d->field(i));
      if (fd->type() == FieldDescriptor::TYPE_MESSAGE) {
        if (fd->is_repeated()) {
          int sz = r->FieldSize(*message, fd);
          for (int i = 0; i < sz; ++i) {
            std::ostringstream oss;
            oss << fd->name() << i;
            show_ranges(oss.str(), &r->GetRepeatedMessage(*message, fd, i));
          }
        } else {
          show_ranges(fd->name(), &r->GetMessage(*message, fd));
        }
      }
    }
  }
}
}
void show_ranges(std::string fname) {
  std::fstream in(fname.c_str(), std::ios::in | std::ios::binary);
  ::npctransport_proto::Configuration input;
  input.ParseFromIstream(&in);
  show_ranges("root", &input);
}

// a range of double values [lb..ub] with <steps> discrete steps,
// assicuated with message m
namespace {
struct Range {
  std::string name;
  double lb, ub;
  int steps;
  double base;
  Message* m;
  IMP_SHOWABLE_INLINE(Range, {
    out << name << " " << lb << "-" << ub << " in " << steps;
  });
};
typedef base::Vector<Range> Ranges;

double get_value(const Reflection* r, const Message* m,
                 const FieldDescriptor* fd) {
  if (fd->type() == FieldDescriptor::TYPE_DOUBLE) {
    return r->GetDouble(*m, fd);
  } else {
    return r->GetInt32(*m, fd);
  }
}
void set_value(const Reflection* r, Message* m, const FieldDescriptor* fd,
               double v) {
  if (fd->type() == FieldDescriptor::TYPE_DOUBLE) {
    return r->SetDouble(m, fd, v);
  } else {
    return r->SetInt32(m, fd, v);
  }
}

// return a vector of ranges for any message that descends from the input
// protobuf message <message> and contains a range of values
// (those are indicated by having 'lower' and 'upper' fields).
// In addition, writes all constant messages to out_message, including
// messages with 'degenerate' ranges (only a lower bound)
//
// @param name the name to be given to the range (range.name),
//             only if it is a direct child of current message.
// @param message the input message
// @param out_message the message to which constand and 'degenerate' range
//                    messages are written (in the same hierarchy as in message)
//                    In addition, range.m (message field) of each range
//                    that is returned is associated with the descendent of
//                    out_message where it should be added to maintain the
//                    correct hierarchy
//        returned (range.m). In addition, constant and 'degenerate' range
//        messages are written to out_message.
//
// @return A vector of ranges. The range.name field of each Range is set to the
//         name of its direct parent message field (or the input parameter
//         <name>, if the range belong directly to the message).
//

Ranges get_ranges(std::string name, const Message* message,
                  Message* out_message) {
  std::string contextname = (std::string("get_ranges: ") + name);
  IMP::base::CreateLogContext gr(contextname.c_str());
  const Reflection* r(message->GetReflection());
  const Reflection* out_r(out_message->GetReflection());
  /*IMP_LOG(VERBOSE, "Inspecting " << name << " "
    << message->GetTypeName() << std::endl);*/
  /*std::cout << "inspecting " << name << " " << message->GetTypeName()
    << std::endl;*/
  Ranges ret;
  const Descriptor* d(message->GetDescriptor());
  const Descriptor* out_d(out_message->GetDescriptor());
  const FieldDescriptor* lfd(d->FindFieldByName("lower"));
  const FieldDescriptor* ufd(d->FindFieldByName("upper"));
  if (lfd && ufd) {  // current message directly contains range
    if (r->HasField(*message, ufd)) {
      IMP_LOG(VERBOSE, "Found range " << name << std::endl);
      Range cur;
      const FieldDescriptor* sfd(d->FindFieldByName("steps"));
      const FieldDescriptor* bfd(d->FindFieldByName("base"));
      cur.name = name;
      cur.lb = get_value(r, message, lfd);
      cur.ub = get_value(r, message, ufd);
      if (lfd->type() == FieldDescriptor::TYPE_INT32) {
        cur.base = 1;  // evenly spaced
      } else {
        cur.base = r->GetDouble(*message, bfd);
      }
      cur.steps = r->GetInt32(*message, sfd);
      cur.m = out_message;
      ret.push_back(cur);
      // IMP_LOG(VERBOSE, "Ret is " << IMP::Showable(ret) << std::endl);
    } else {
      IMP_LOG(VERBOSE, "Found value " << name << std::endl);
      const FieldDescriptor* out_lfd(out_d->FindFieldByName("value"));
      set_value(out_r, out_message, out_lfd, get_value(r, message, lfd));
    }
  } else {  // recursively look for descendent messages with ranges
    for (int i = 0; i < d->field_count(); ++i) {
      const FieldDescriptor* fd(d->field(i));
      const FieldDescriptor* out_fd(out_d->FindFieldByName(fd->name()));
      IMP_INTERNAL_CHECK(out_fd, "No field named " << fd->name());
      if (fd->type() == FieldDescriptor::TYPE_MESSAGE) {
        if (fd->is_repeated()) {
          int sz = r->FieldSize(*message, fd);
          for (int i = 0; i < sz; ++i) {
            ret +=
                get_ranges(fd->name(), &r->GetRepeatedMessage(*message, fd, i),
                           out_r->AddMessage(out_message, out_fd));
            // IMP_LOG(VERBOSE, "Got " << IMP::Showable(ret) << std::endl);
          }
        } else {
          ret += get_ranges(fd->name(), &r->GetMessage(*message, fd),
                            out_r->MutableMessage(out_message, out_fd));
          // IMP_LOG(VERBOSE, "Got " << IMP::Showable(ret) << std::endl);
        }
      } else {
        if (out_fd) {
          IMP_LOG(VERBOSE, "Found constant " << fd->name() << std::endl);
          if (out_fd->type() == FieldDescriptor::TYPE_STRING) {
            out_r->SetString(out_message, out_fd, r->GetString(*message, fd));
          } else {
            set_value(out_r, out_message, out_fd, get_value(r, message, fd));
          }
        }
      }
    }
  }
  // IMP_LOG(VERBOSE, "Returning " << IMP::Showable(ret) << std::endl);
  return ret;
}
template <class R, class Out>
void copy(R r, Out out) {
  std::copy(r.begin(), r.end(), out);
}

// Compute the [work_unit]'th combination of values from the set of ranges r
// using r[i].steps in each dimension, with points evenly distributed for each
// entry in r[i] on a log scale, using r[i].base log base.
//
// It is guaranteed that if r induces k possible combinations,
// then iterating over work_unit from 0 to (k-1) will enumerate every possible
// combination of these values.
// The enumeration is cyclic, that is, work_unit % k and work_unit will return
// the same result (over the same range)
//
// @param[in] r a vector with the set of value ranges
// @param[in] work_unit the index of the combination of ranges that will be
//                      stored to values and indexesd. This is a cyclic value.
// @param[out] values The emumerated combination of values within the ranges
//                    defined by r
// @param[out] indexes The emumerated combination of indexes for work_unit
// @param[in]  show_steps whether to display the grid
int assign_internal(const Ranges& r, int work_unit, Floats& values,
                    Ints& indexes, bool show_steps) {
  typedef IMP::algebra::DenseGridStorageD<-1, int> Storage;
  typedef IMP::algebra::LogEmbeddingD<-1> Embedding;
  typedef IMP::algebra::GridD<-1, Storage, int, Embedding> Grid;
  IMP::Ints steps(r.size());
  // read bounding box for grid from r
  IMP::algebra::VectorKD lb = IMP::algebra::get_zero_vector_kd(r.size());
  IMP::algebra::VectorKD ub = IMP::algebra::get_zero_vector_kd(r.size());
  IMP::algebra::VectorKD factors(IMP::Floats(lb.get_dimension(), 2.0));
  for (unsigned int i = 0; i < r.size(); ++i) {
    lb[i] = r[i].lb;
    ub[i] = r[i].ub;
    factors[i] = r[i].base;
    steps[i] = r[i].steps;
  }
  Grid g(Storage(steps, 0),
         Embedding(IMP::algebra::BoundingBoxKD(IMP::algebra::VectorKD(lb),
                                               IMP::algebra::VectorKD(ub)),
                   factors, steps, true));
  if (show_steps) {
    Grid::ExtendedIndex ei(Ints(g.get_dimension(), 0));
    for (unsigned int i = 0; i < g.get_dimension(); ++i) {
      std::cout << r[i].name << ": ";
      for (unsigned int j = 0; j < g.get_number_of_voxels(i); ++j) {
        ei[i] = j;
        std::cout << g.get_center(ei)[i] << " ";
      }
      std::cout << std::endl;
      ei[i] = 0;
    }
  }
  Grid::AllIndexIterator it = g.all_indexes_begin();
  unsigned int nv = std::distance(g.all_indexes_begin(), g.all_indexes_end());
  std::advance(it, work_unit % nv);
  IMP::algebra::VectorKD center = g.get_center(*it);
  values = Floats(center.coordinates_begin(), center.coordinates_end());
  indexes.resize(values.size());
  // it-> returns a copy
  copy(*it, indexes.begin());
  return g.get_number_of_voxels();
}
}

// see documentation in .h file
int assign_ranges(std::string fname, std::string ofname, unsigned int work_unit,
                  bool show_steps, boost::uint64_t random_seed) {
  IMP_FUNCTION_LOG;
  std::fstream in(fname.c_str(), std::ios::in | std::ios::binary);
  if (!in) {
    IMP_THROW("Could not open file " << fname, IOException);
  }
  ::npctransport_proto::Configuration input;
  bool success = input.ParseFromIstream(&in);
  if (!success) {
    IMP_THROW("Unable to read from protobuf " << fname, IOException);
  }
  npctransport_proto::Output output_output;
  npctransport_proto::Assignment& output = *output_output.mutable_assignment();
  base::SetLogState sls(base::VERBOSE);
  Ranges ranges = get_ranges("all", &input, &output);
  /*for (unsigned int i=0; i< ranges.size(); ++i) {
    std::cout << ranges[i].lb << " " << ranges[i].ub << " " << ranges[i].steps
              << " " << ranges[i].base << std::endl;
              }*/
  int ret = 0;
  if (ranges.empty()) {
    IMP_WARN("No message with value ranges detected for file " << fname);
    // return 0;
  } else {
    Floats values;
    Ints indexes;
    // assign the [work_unit]'th combination of ranges into the
    // message indicated for each range entry (ranges[i].m)
    ret = assign_internal(ranges, work_unit, values, indexes, show_steps);
    for (unsigned int i = 0; i < ranges.size(); ++i) {
      const Reflection* r(ranges[i].m->GetReflection());
      const Descriptor* d(ranges[i].m->GetDescriptor());
      IMP_LOG(VERBOSE, "Assigning range " << ranges[i].name << std::endl);
      const FieldDescriptor* vfd(d->FindFieldByName("value"));
      set_value(r, ranges[i].m, vfd, values[i]);
      const FieldDescriptor* ifd(d->FindFieldByName("index"));
      IMP_INTERNAL_CHECK(ifd, "No index found?");
      r->SetInt32(ranges[i].m, ifd, indexes[i]);
    }
  }
  // output work units and automatic parameters
  double time_step = get_time_step(output);
  std::cout << "DEBUG protobug auto" << std::endl;
  output.set_time_step(time_step);
  output.set_work_unit(work_unit);
  output.set_number_of_frames(get_number_of_frames(output, time_step));
  output.set_dump_interval_frames(
      get_dump_interval_in_frames(output, time_step));
  output.set_statistics_interval_frames(
      get_statistics_interval_in_frames(output, time_step));
  output.set_range(get_close_pairs_range(output));
  output.set_random_seed(random_seed);
  std::fstream out(ofname.c_str(), std::ios::out | std::ios::binary);
  if (!out) {
    IMP_THROW("Could not open file " << ofname, IOException);
  }

  {
    for (int i = 0; i < output.fgs_size(); ++i) {
      if (output.fgs(i).number().value() != 0) {
        output_output.mutable_statistics()->add_fgs();
        IMP_USAGE_CHECK(
            output_output.mutable_statistics()->fgs_size() ==
                static_cast<int>(i + 1),
            "Wrong size: " << output_output.mutable_statistics()->fgs_size());
      }
    }
    for (int i = 0; i < output.floaters_size(); ++i) {
      if (output.floaters(i).number().value() != 0) {
        output_output.mutable_statistics()->add_floaters();
      }
    }
    for (int i = 0; i < output.interactions_size(); i++) {
      if (output.interactions(i).is_on().value()) {
        ::npctransport_proto::Statistics_InteractionStats* pOutStats_i =
            output_output.mutable_statistics()->add_interactions();
        pOutStats_i->set_type0(output.interactions(i).type0());
        pOutStats_i->set_type1(output.interactions(i).type1());
      }
    }
  }

  bool written = output_output.SerializeToOstream(&out);
  if (!written) {
    IMP_THROW("Unable to write to " << ofname, IOException);
  }

  return ret;
}

int get_number_of_work_units(std::string assignment_file) {
  ::npctransport_proto::Configuration input;
  std::fstream in(assignment_file.c_str(), std::ios::in | std::ios::binary);
  if (!in) {
    IMP_THROW("Could not open file " << assignment_file, IOException);
  }
  input.ParseFromIstream(&in);
  npctransport_proto::Assignment output;
  base::SetLogState sls(base::VERBOSE);
  Ranges ranges = get_ranges("all", &input, &output);
  Floats values;
  Ints indexes;
  int ret = assign_internal(ranges, 0, values, indexes, false);
  return ret;
}

IMPNPCTRANSPORT_END_NAMESPACE
