/**
 *  \file io.cpp
 *  \brief description.
 *
 *  Copyright 2007-2020 IMP Inventors. All rights reserved.
 *
 */

// IMP headers:
#include <IMP/npctransport/protobuf.h>
#include <IMP/npctransport/automatic_parameters.h>
#include <IMP/npctransport/internal/npctransport.pb.h>
#include <IMP/SingletonContainer.h>
#include <IMP/utility.h>
#include <IMP/algebra/GridD.h>
#include <IMP/atom/estimates.h>
#include <IMP/algebra/grid_storages.h>
#include <IMP/CreateLogContext.h>
#include <IMP/SetLogState.h>
// Protobuf headers:
#include <google/protobuf/descriptor.h>
#include <google/protobuf/message.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/coded_stream.h>
// C++ and boost headers:
#include <boost/scoped_ptr.hpp>
#include <fstream>
#include <iostream>
#include <fcntl.h>
#if defined(_MSC_VER)
#include <io.h>
#else
#include <sys/types.h>
#include <sys/stat.h>
#endif


IMPNPCTRANSPORT_BEGIN_NAMESPACE
using namespace ::google::protobuf;
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
        IMP_LOG(PROGRESS, "Found " << name << std::endl);
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
  ::npctransport_proto::Configuration config;
  config.ParseFromIstream(&in);
  show_ranges("root", &config);
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
  typedef Vector<Range> Ranges;

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
  double get_repeated_value
  (const Reflection* r, const Message* m,
   const FieldDescriptor* fd,
   int i)
  {
    if (fd->type() == FieldDescriptor::TYPE_DOUBLE) {
      return r->GetRepeatedDouble(*m, fd, i);
    } else {
      return r->GetRepeatedInt32(*m, fd, i);
    }
  }
  void add_repeated_value
  (const Reflection* r, Message* m, const FieldDescriptor* fd,
   double v)
  {
    if (fd->type() == FieldDescriptor::TYPE_DOUBLE) {
      return r->AddDouble(m, fd, v);
    } else {
      return r->AddInt32(m, fd, v);
    }
  }


  // return a vector of ranges for any message that descends from the config
  // protobuf message <message> and contains a range of values
  // (those are indicated by having 'lower' and 'upper' fields).
  // In addition, writes all constant messages to out_message, including
  // messages with 'degenerate' ranges (only a lower bound).
  // Copying is by the field name (not its index!)
  //
  // @param name the name to be given to the range (range.name),
  //             only if it is a direct child of current message.
  // @param in_message the config message
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
  //         name of its direct parent message field (or the config parameter
  //         <name>, if the range belong directly to the message).
  //

  Ranges get_ranges(std::string name, const Message* in_message,
                    Message* out_message) {
    std::string contextname = (std::string("get_ranges: ") + name);
    IMP::CreateLogContext gr(contextname.c_str());
    const Reflection* in_r(in_message->GetReflection());
    const Reflection* out_r(out_message->GetReflection());
    IMP_LOG(VERBOSE,"Inspecting in message " << name << " "
      << in_message->GetTypeName() << std::endl);
    Ranges ret;
    const Descriptor* in_d(in_message->GetDescriptor());
    const Descriptor* out_d(out_message->GetDescriptor());
    const FieldDescriptor* lfd(in_d->FindFieldByName("lower"));
    const FieldDescriptor* ufd(in_d->FindFieldByName("upper"));
    if (lfd && ufd) {
      // Handle message that directly contains range (actual or degenerate):
      if (in_r->HasField(*in_message, ufd)) {
        // Actual range from lower to upper:
        IMP_LOG(VERBOSE, "Found range " << name << std::endl);
        Range cur;
        const FieldDescriptor* sfd(in_d->FindFieldByName("steps"));
        const FieldDescriptor* bfd(in_d->FindFieldByName("base"));
        cur.name = name;
        cur.lb = get_value(in_r, in_message, lfd);
        cur.ub = get_value(in_r, in_message, ufd);
        //        if (lfd->type() == FieldDescriptor::TYPE_INT32) {
        //   cur.base = 1;  // evenly spaced
        //} else {
        cur.base = in_r->GetDouble(*in_message, bfd);
        //}
        cur.steps = in_r->GetInt32(*in_message, sfd);
        cur.m = out_message;
        ret.push_back(cur);
        // IMP_LOG(VERBOSE, "Ret is " << IMP::Showable(ret) << std::endl);
      } else {
        // Degenerate range with a single value (=lower):
        IMP_LOG(VERBOSE, "Found value " << name << std::endl);
        const FieldDescriptor* out_lfd(out_d->FindFieldByName("value"));
        set_value(out_r, out_message, out_lfd, get_value(in_r, in_message, lfd));
      }
    } else { // message does not directly contain a range
      // Recursively look for descendent messages with ranges:
      for (int i = 0; i < in_d->field_count(); ++i) {
        const FieldDescriptor* in_fd(in_d->field(i));
        if(!in_fd->is_repeated()){
          if(!in_r->HasField(*in_message, in_fd)){
            continue; // skip unset field
          }
        }
        const FieldDescriptor* out_fd(out_d->FindFieldByName(in_fd->name()));
        IMP_INTERNAL_CHECK(out_fd, "No field named " << in_fd->name()
                           << " in assignment message");
        if (in_fd->type() == FieldDescriptor::TYPE_MESSAGE) {
          if (in_fd->is_repeated()) {
            int sz = in_r->FieldSize(*in_message, in_fd);
            for (int i = 0; i < sz; ++i) {
              ret +=
                get_ranges(in_fd->name(), &in_r->GetRepeatedMessage(*in_message, in_fd, i),
                           out_r->AddMessage(out_message, out_fd));
              // IMP_LOG(VERBOSE, "Got " << IMP::Showable(ret) << std::endl);
            }
          } else { // not repeated:
            ret += get_ranges(in_fd->name(), &in_r->GetMessage(*in_message, in_fd),
                              out_r->MutableMessage(out_message, out_fd));
            // IMP_LOG(VERBOSE, "Got " << IMP::Showable(ret) << std::endl);
          }
        } else { // not a message:
          if (out_fd) {
            IMP_LOG(VERBOSE, "Found constant " << in_fd->name() << std::endl);
            if (out_fd->type() == FieldDescriptor::TYPE_STRING) {
              if (in_fd->is_repeated()) {
                int sz = in_r->FieldSize(*in_message, in_fd);
                for(int i = 0; i < sz; ++i) {
                  std::string str= in_r->GetRepeatedString(*in_message, in_fd, i);
                  out_r->AddString(out_message, out_fd, str);
                }
              } else { // not a repeated string:
                out_r->SetString(out_message, out_fd,
                                 in_r->GetString(*in_message, in_fd));
              }
            } else { // not a string:
              if (in_fd->is_repeated()) {
                int sz= in_r->FieldSize(*in_message, in_fd);
                for(int i= 0; i<sz; ++i) {
                  add_repeated_value(out_r, out_message, out_fd,
                                     get_repeated_value(in_r, in_message, in_fd, i));
                }
              } else { // not a repeated field
                set_value(out_r, out_message, out_fd,
                          get_value(in_r, in_message, in_fd));
              }
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
      Ints dzeros(g.get_dimension(), 0);
      Grid::ExtendedIndex ei(dzeros.begin(), dzeros.end());
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
    values = Floats(center.begin(), center.end());
    indexes.resize(values.size());
    // it-> returns a copy
    copy(*it, indexes.begin());
    return g.get_number_of_voxels();
  }
} // anonymous namespace


namespace {

  //! set default value of protobuf message 'prefix' to passed dafault value
  //! the value must be positive (otherwise default value is used)
#define SET_DEFAULT_POSITIVE_FLOAT(prefix, default_value)                              \
  {                                                                     \
    bool has_value = ( *prefix ).has_value();                           \
    bool invalid_value = has_value ? ( *prefix ).value() <= 0.0 : true; \
    if(invalid_value){                                                  \
      ( *prefix ).set_value(default_value);                             \
    }                                                                   \
  }

  // set default options for temperature accelerated MD version
  void set_default_tamd_options(::npctransport_proto::Assignment& assign) {
    unsigned int n = assign.fgs_size();
    for(unsigned int i = 0 ; i < n ; i++) {
      ::npctransport_proto::Assignment_FGAssignment* mfg =
        assign.mutable_fgs(i);
      if(mfg->is_tamd()){
        SET_DEFAULT_POSITIVE_FLOAT(mfg->mutable_tamd_t_factor_coeff(), 1.0);
        SET_DEFAULT_POSITIVE_FLOAT(mfg->mutable_tamd_t_factor_base(), 1.0);
        SET_DEFAULT_POSITIVE_FLOAT(mfg->mutable_tamd_f_factor_coeff(), 1.0);
        SET_DEFAULT_POSITIVE_FLOAT(mfg->mutable_tamd_f_factor_base(), 1.0);
        SET_DEFAULT_POSITIVE_FLOAT(mfg->mutable_tamd_k(), 1.0);
      } // if mfg->is_tamd
    } // for
  }
#undef SET_DEFAULT_POSITIVE_FLOAT
}; // anonymous namespace


// see documentation in .h file
int assign_ranges(std::string ifname, std::string ofname, unsigned int work_unit,
                  bool show_steps, boost::uint64_t random_seed) {
  IMP_FUNCTION_LOG;
  std::fstream in(ifname.c_str(), std::ios::in | std::ios::binary);
  if (!in) {\
    IMP_THROW("Could not open file " << ifname, IOException);
  }
  ::npctransport_proto::Configuration config;
  bool success = config.ParseFromIstream(&in);
  if (!success) {
    IMP_THROW("Unable to read from protobuf " << ifname, IOException);
  }
  npctransport_proto::Output output;
  npctransport_proto::Assignment& assignment = *output.mutable_assignment();
  output.mutable_statistics(); // create if not there
  //  SetLogState sls(WARNING);
  Ranges ranges = get_ranges("all", &config, &assignment);
  /*for (unsigned int i=0; i< ranges.size(); ++i) {
    std::cout << ranges[i].lb << " " << ranges[i].ub << " " << ranges[i].steps
              << " " << ranges[i].base << std::endl;
              }*/
  int ret = 0;
  if (ranges.empty()) {
    IMP_WARN("No message with value ranges detected for file '" << ifname
             << "' - this is probably fine" << std::endl);
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
  // assignment work units and automatic parameters
  double max_trans_relative_to_radius = 0.05; // TODO: param?
  double time_step = get_time_step(assignment, max_trans_relative_to_radius);
  assignment.set_time_step(time_step);
  assignment.set_work_unit(work_unit);
  assignment.set_number_of_frames(get_number_of_frames(assignment, time_step));
  assignment.set_dump_interval_frames
    ( get_dump_interval_in_frames(assignment, time_step) );
  assignment.set_statistics_interval_frames
    ( get_statistics_interval_in_frames(assignment, time_step) );
  assignment.set_output_statistics_interval_frames
    ( get_output_statistics_interval_in_frames(assignment, time_step) );
  assignment.set_range(get_close_pairs_range(assignment));
  IMP_LOG(VERBOSE, "Close pair sphere-sphere distance range in A: "
            << assignment.range() << std::endl);
  assignment.set_random_seed(random_seed);

  set_default_tamd_options(assignment);
  std::fstream out(ofname.c_str(), std::ios::out | std::ios::binary);
  if (!out) {
    IMP_THROW("Could not open file " << ofname, IOException);
  }

  // Fill in types for fgs and floaters if needed
  // TODO: this is just for backward support -
  //       should be removed once deprecated
  {
    for (int i = 0; i < assignment.fgs_size(); ++i)
      {
        // store default type if one does not exist
        bool has_type =  assignment.fgs(i).has_type();
        if(has_type) {
          has_type = ( assignment.fgs(i).type() != "" );
        }
        IMP_ALWAYS_CHECK(has_type, "fg " << i << " lacking type",
                         ValueException);
      } // for i
    for (int i = 0; i < assignment.floaters_size(); ++i)
      {
        bool has_type = assignment.floaters(i).has_type();
        if(has_type) {
          has_type = (assignment.floaters(i).type() != "" );
        }
        IMP_ALWAYS_CHECK(has_type, "floater " << i << " lacking type",
                         ValueException);
      } // for i
  }

  bool written = output.SerializeToOstream(&out);
  if (!written) {
    IMP_THROW("Unable to write to " << ofname, IOException);
  }

  return ret;
}

int get_number_of_work_units(std::string assignment_file) {
  ::npctransport_proto::Configuration config;
  std::fstream in(assignment_file.c_str(), std::ios::in | std::ios::binary);
  if (!in) {
    IMP_THROW("Could not open file " << assignment_file, IOException);
  }
  config.ParseFromIstream(&in);
  npctransport_proto::Assignment assignment;
  SetLogState sls(VERBOSE);
  Ranges ranges = get_ranges("all", &config, &assignment);
  Floats values;
  Ints indexes;
  int ret = assign_internal(ranges, 0, values, indexes, false);
  return ret;
}

void load_pb_conformation
( const ::npctransport_proto::Conformation &conformation,
  IMP::SingletonContainerAdaptor beads,
  boost::unordered_map<core::ParticleType, algebra::Sphere3Ds> &sites)
{
  IMP_CONTAINER_FOREACH(IMP::SingletonContainer, beads, {
      const ::npctransport_proto::Conformation_Particle &pcur =
        conformation.particle(_2);
      core::RigidBody rb(beads->get_model(), _1);
      algebra::Vector3D translation(pcur.x(), pcur.y(), pcur.z());
      algebra::Rotation3D rotation
        ( algebra::Vector4D(pcur.r(), pcur.i(), pcur.j(), pcur.k()) );

      algebra::Transformation3D tr(rotation, translation);
      rb.set_reference_frame(algebra::ReferenceFrame3D(tr));
    });
  for (int i = 0; i < conformation.sites_size(); ++i) {
    const ::npctransport_proto::Conformation::Sites &cur =
      conformation.sites(i);
    core::ParticleType pt(cur.name());
    for (int j = 0; j < cur.coordinates_size(); ++j) {
      const ::npctransport_proto::Conformation::Coordinates &coords =
        cur.coordinates(j);
      sites[pt][j]._set_center // might not work?
        (algebra::Vector3D(coords.x(), coords.y(), coords.z()));
      // radius?
    }
  }
}

void save_pb_conformation
( SingletonContainerAdaptor beads,
  const boost::unordered_map<core::ParticleType, algebra::Sphere3Ds> &sites,
  ::npctransport_proto::Conformation *conformation )
{
  conformation->clear_sites();
  conformation->clear_particle();
  IMP_CONTAINER_FOREACH
    (IMP::SingletonContainer, beads,
     {
       ::npctransport_proto::Conformation_Particle *pcur =
         conformation->add_particle();
       core::RigidBody rb(beads->get_model(), _1);
       algebra::Transformation3D tr =
         rb.get_reference_frame().get_transformation_to();
       pcur->set_x(tr.get_translation()[0]);
       pcur->set_y(tr.get_translation()[1]);
       pcur->set_z(tr.get_translation()[2]);
       pcur->set_r(tr.get_rotation().get_quaternion()[0]);
       pcur->set_i(tr.get_rotation().get_quaternion()[1]);
       pcur->set_j(tr.get_rotation().get_quaternion()[2]);
       pcur->set_k(tr.get_rotation().get_quaternion()[3]);
     });
  typedef boost::unordered_map<core::ParticleType, algebra::Sphere3Ds> M;
  for (M::const_iterator it = sites.begin(); it != sites.end(); it++)
    {
      ::npctransport_proto::Conformation::Sites *cur
        = conformation->add_sites();
      cur->set_name(it->first.get_string());
      algebra::Sphere3Ds spheres = it->second;
      for (algebra::Sphere3Ds::const_iterator sphere = spheres.begin();
           sphere != spheres.end(); sphere++)
        {
          ::npctransport_proto::Conformation::Coordinates *out_coords =
            cur->add_coordinates();
          algebra::Vector3D coord = sphere->get_center();
          out_coords->set_x(coord[0]);
          out_coords->set_y(coord[1]);
          out_coords->set_z(coord[2]);
        }
    }
}

//! load file output_fname into protobuf output object output
bool load_output_protobuf
(std::string output_fname,
 ::npctransport_proto::Output& output)
{
  bool is_ok(false);
  int fd=IMP_C_OPEN(output_fname.c_str(),
                    IMP_C_OPEN_FLAG(O_RDONLY) | IMP_C_OPEN_BINARY);
  if(fd!=-1) {
    google::protobuf::io::FileInputStream fis(fd);
    google::protobuf::io::CodedInputStream cis(&fis);
    cis.SetTotalBytesLimit(500000000,200000000);
    is_ok=output.ParseFromCodedStream(&cis);
    IMP_C_CLOSE(fd);
  }
  return is_ok;
}



IMPNPCTRANSPORT_END_NAMESPACE
