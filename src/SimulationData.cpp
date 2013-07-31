/**
 *  \file SimulationData.cpp
 *  \brief description.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/SimulationData.h>
#include <IMP/npctransport/SitesGeometry.h>
#include <IMP/npctransport/SlabGeometry.h>
#include <IMP/npctransport/SlabSingletonScore.h>
#include <IMP/npctransport/particle_types.h>
#include <IMP/npctransport/protobuf.h>
#ifdef IMP_NPC_GOOGLE
IMP_GCC_PUSH_POP(diagnostic push)
IMP_GCC_PRAGMA(diagnostic ignored "-Wsign-compare")
#include "third_party/npc/npctransport/data/npctransport.pb.h"
IMP_GCC_PUSH_POP(diagnostic pop)
#else
#include <IMP/npctransport/internal/npctransport.pb.h>
#endif
#include <IMP/npctransport/creating_particles.h>
#include <IMP/npctransport/io.h>
#include <IMP/npctransport/typedefs.h>
#include <IMP/algebra/vector_generators.h>
#include <IMP/atom/estimates.h>
#include <IMP/atom/distance.h>
#include <IMP/atom/Diffusion.h>
#include <IMP/atom/Selection.h>
#include <IMP/base/log.h>
#include <IMP/core/HarmonicUpperBound.h>
#include <IMP/core/pair_predicates.h>
#include <IMP/core/RestraintsScoringFunction.h>
#include <IMP/core/XYZR.h>
#include <IMP/core/generic.h>
#include <IMP/display/LogOptimizerState.h>
#include <IMP/display/PymolWriter.h>
#include <IMP/display/primitive_geometries.h>
#include <IMP/display/restraint_geometry.h>
//#include <IMP/example/optimizing.h>
//#include <IMP/example/randomizing.h>
#include <IMP/npctransport/rmf_links.h>
#include <RMF/FileHandle.h>
#include <RMF/FileConstHandle.h>
#include <IMP/rmf/atom_io.h>
#include <IMP/rmf/frames.h>
#include <numeric>

#include <set>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
#define GET_ASSIGNMENT(name) name##_ = pb_assignment.name().value()
#define GET_VALUE(name) name##_ = pb_assignment.name()

SimulationData::SimulationData(std::string output_file, bool quick,
                               std::string rmf_file_name)
: Object("SimulationData%1%"),
  m_(nullptr),
  bd_(nullptr),
  scoring_(nullptr),
  diffusers_changed_(false),
  obstacles_changed_(false),
  optimizable_diffusers_(nullptr),
  diffusers_(nullptr),
  rmf_file_name_(rmf_file_name),
  is_stats_reset_(false)
{
  initialize(output_file, quick);
}

void SimulationData::initialize(std::string output_file, bool quick) {
  ::npctransport_proto::Output pb_data_;
  output_file_name_ = output_file;
  std::ifstream file(output_file_name_.c_str(), std::ios::binary);
  bool read = pb_data_.ParseFromIstream(&file);
  if (!read) {
    IMP_THROW("Unable to read data from protobuf" << output_file_name_,
              base::IOException);
  }
  pb_data_.mutable_statistics(); // create it if not there
  const ::npctransport_proto::Assignment &
      pb_assignment= pb_data_.assignment();
  GET_ASSIGNMENT(box_side);
  GET_ASSIGNMENT(tunnel_radius);
  GET_ASSIGNMENT(slab_thickness);
  GET_ASSIGNMENT(slab_is_on);
  GET_ASSIGNMENT(box_is_on);
  GET_VALUE(number_of_trials);
  GET_VALUE(number_of_frames);
  GET_VALUE(dump_interval_frames);
  GET_ASSIGNMENT(angular_d_factor);
  GET_VALUE(range);
  GET_VALUE(statistics_interval_frames);
  GET_VALUE(time_step);
  GET_ASSIGNMENT(statistics_fraction);
  GET_VALUE(maximum_number_of_minutes);
  if (quick) {
    number_of_frames_ = 2;
    number_of_trials_ = 1;
  }
  // initialize scoring
  scoring_ = new Scoring(this, pb_assignment);

  // create particles hierarchy
  root_ = new Particle(get_model());
  root_->add_attribute(get_simulation_data_key(), this);
  atom::Hierarchy hr = atom::Hierarchy::setup_particle(root_);
  root_->set_name("root");
  for (int i = 0; i < pb_assignment.fgs_size(); ++i) {
    if (!pb_assignment.fgs(i).has_type()) {
      std::string type_i = type_of_fg[i].get_string();
      pb_data_.mutable_assignment()->mutable_fgs(i)
        ->set_type(type_i);
    }
    create_fgs(pb_assignment.fgs(i));
  }
  for (int i = 0; i < pb_assignment.floaters_size(); ++i) {
    create_floaters(pb_assignment.floaters(i), type_of_float[i],
                    display::get_display_color(i));
  }
  for (int i = 0; i < pb_assignment.obstacles_size(); ++i) {
    create_obstacles(pb_assignment.obstacles(i));
  }

  IMP_LOG(TERSE, "   SimulationData before adding interactions" << std::endl);
  for (int i = 0; i < pb_assignment.interactions_size(); ++i) {
    const ::npctransport_proto::Assignment_InteractionAssignment &
        interaction_i = pb_assignment.interactions(i);
    add_interaction(interaction_i);
  }

  // bounding box / slab constraints on diffusers

  if (pb_data_.has_rmf_conformation()) {
    RMF::FileConstHandle fh =
        RMF::open_rmf_buffer_read_only(pb_data_.rmf_conformation());
    initialize_positions_from_rmf(fh, 0);
    // load from output file
  } else if (pb_data_.has_conformation()) {
    std::cout << "Loading from output file " << std::endl ;
    load_pb_conformation(pb_data_.conformation(), get_diffusers(), sites_);
  }
  if (pb_data_.has_statistics()) {
    if (pb_data_.statistics().has_bd_simulation_time_ns()) {
      const double fs_in_ns = 1.0E+6;
      get_bd()->set_current_time
        (pb_data_.statistics().bd_simulation_time_ns()
         * fs_in_ns);
    }
  }
}


/**
   Adds the FG Nup chains to the model hierarchy,
   based on the settings in data
*/
void SimulationData::create_fgs
( const ::npctransport_proto::Assignment_FGAssignment &fg_data)
{
  // set type
  core:: ParticleType type(fg_data.type());
  if (fg_data.number().value() > 0) {
    fg_types_.insert(type);
    fgs_stats_.push_back(base::Vector<BodyStatisticsOptimizerStates>());
    chain_stats_.push_back(ChainStatisticsOptimizerStates());
    ParticlesTemp fg_particles;
    IMP::Particle* p_fg_root = new IMP::Particle( get_model() );
    atom::Hierarchy fg_root = atom::Hierarchy::setup_particle(p_fg_root);
    fg_root->set_name(type.get_string());
    atom::Hierarchy(get_root()).add_child(fg_root);
    for (int j = 0; j < fg_data.number().value(); ++j) {
      IMP::Particle* pc=
        create_fg_chain(this, fg_data, display::Color(.3, .3, .3));
      atom::Hierarchy hc(pc);
      fg_root.add_child(hc);
      fg_particles.push_back(pc);
      ParticlesTemp chain_beads = hc.get_children();
      // set chain anchors if specified
      if (fg_data.anchor_coordinates_size() > j) {
        ::npctransport_proto::Assignment_XYZ xyz =
          fg_data.anchor_coordinates(j);
        core::XYZ d(chain_beads[0]);
        d.set_coordinates(algebra::Vector3D(xyz.x(), xyz.y(), xyz.z()));
        d.set_coordinates_are_optimized(false);
      }
      // stats
      chain_stats_.back().push_back
        ( new ChainStatisticsOptimizerState(chain_beads) );
      chain_stats_.back().back()->set_period(statistics_interval_frames_);
      fgs_stats_.back().push_back(BodyStatisticsOptimizerStates());
      for (unsigned int k = 0; k < chain_beads.size(); ++k) {
        fgs_stats_.back().back()
            .push_back( new BodyStatisticsOptimizerState(chain_beads[k]) );
        fgs_stats_.back().back().back()
            ->set_period(statistics_interval_frames_);
      }
    } // for j
    particles_[type] = fg_particles;
    // add sites for this type
    if (fg_data.interactions().value() > 0) {
      int nsites = fg_data.interactions().value();
      set_sites(type, nsites, fg_data.radius().value());
    }
    // add general scoring scale factors
    get_scoring()->set_interaction_range_factor
    ( type, fg_data.interaction_range_factor().value());
    get_scoring()->set_interaction_k_factor
      ( type, fg_data.interaction_k_factor().value());
    set_diffusers_changed(true);
  } // if
}

atom::Hierarchies SimulationData::get_fg_chains(atom::Hierarchy root) const
{
  IMP_INTERNAL_CHECK(root, "root for SimulationData::get_fg_chains() is null");
  // TODO: maybe FGs should just be marked by a decorator, and not identified by
  //       type alone
  atom::Hierarchies ret;
  if (root.get_number_of_children() == 0) {
    return ret;
  }
  // I. return root itself if the type of its first direct child
  // is contained in [type_of_fg]
  atom::Hierarchy c = root.get_child(0);
  if (core::Typed::get_is_setup(c)) {
    core::ParticleType t = core::Typed(c).get_type();
    if (fg_types_.find( t ) != fg_types_.end() ) {
      return atom::Hierarchies(1, root);
    }
  }
  // II. If not returned yet, recurse on all of root's children
  for (unsigned int i = 0; i < root.get_number_of_children(); ++i) {
    ret += get_fg_chains(root.get_child(i));
  }
  return ret;
}


/**
   Adds the 'floaters' (free diffusing particles) to the model hierarchy,
   based on the settings in data
*/
void SimulationData::create_floaters(
    const ::npctransport_proto::Assignment_FloaterAssignment &f_data,
    core::ParticleType default_type, display::Color color) {
  core:: ParticleType type(default_type);
  if (f_data.has_type()){
    type = core::ParticleType(f_data.type());
  }
  if (f_data.number().value() > 0) {
    // prepare statistics for this type of floaters:
    float_stats_.push_back(BodyStatisticsOptimizerStates());
    if (slab_is_on_) {  // if has tunnel, create a list of particle stats
      float_transport_stats_.push_back(
          ParticleTransportStatisticsOptimizerStates());
    }
    // create a sub hierarchy with this type of floaters:
    atom::Hierarchy cur_root =
        atom::Hierarchy::setup_particle(new Particle(get_model()));
    IMP_LOG(WARNING,  "   type " << type.get_string() << std::endl);
    cur_root->set_name(type.get_string());
    atom::Hierarchy(get_root()).add_child(cur_root);
    // populate hierarchy with particles:
    ParticlesTemp cur_particles;
    for (int j = 0; j < f_data.number().value(); ++j) {
      double dc = f_data.d_factor().value();
      Particle *cur_p =
          create_particle(this, f_data.radius().value(),
                          angular_d_factor_, dc,
                          color, type);
      cur_particles.push_back(cur_p);
      cur_root.add_child(atom::Hierarchy::setup_particle(cur_p));
      // add particle statistics:
      IMP_NEW(BodyStatisticsOptimizerState, bsos, (cur_p));
      bsos->set_period(statistics_interval_frames_);
      float_stats_.back().push_back(bsos);
      if (slab_is_on_) {  // only if has tunnel
        IMP_NEW(ParticleTransportStatisticsOptimizerState, ptsos,
                (cur_p, -0.5 * slab_thickness_,  // tunnel bottom
                 0.5 * slab_thickness_)          // tunnel top
                );
        ptsos->set_period(statistics_interval_frames_);
        float_transport_stats_.back().push_back(ptsos);
      }
    }
    particles_[type] = cur_particles;
    // add interaction sites to particles of this type:
    if (f_data.interactions().value() > 0) {
      int nsites = f_data.interactions().value();
      IMP_LOG(WARNING, nsites << " sites added " << std::endl);
      set_sites(type, nsites, f_data.radius().value());
    }
    // add type-specific scoring scale factors
    get_scoring()->set_interaction_range_factor
    ( type, f_data.interaction_range_factor().value());
    get_scoring()->set_interaction_k_factor
      ( type, f_data.interaction_k_factor().value());

    set_diffusers_changed(true);
  }
}

/**
   Adds the 'obstacles' (possibly static e.g. nups that make the pore)
   to the model hierarchy, based on the settings in data
*/
void SimulationData::create_obstacles
( const ::npctransport_proto::Assignment_ObstacleAssignment &o_data)
{
  core::ParticleType type = core::ParticleType(o_data.type());
  // create a sub hierarchy with this type of floaters:
  atom::Hierarchy cur_root =
    atom::Hierarchy::setup_particle(new Particle(get_model()));
  IMP_LOG(WARNING,  "   obstacle type " << type.get_string() << std::endl);
  cur_root->set_name(type.get_string());
  atom::Hierarchy(get_root()).add_child(cur_root);
  // populate hierarchy with particles:
  ParticlesTemp cur_particles;
  for (int j = 0; j < o_data.xyzs_size(); ++j) {
      double dc = o_data.d_factor().value();
      Particle *cur_p =
        create_particle(this, o_data.radius().value(),
                        angular_d_factor_, dc,
                        display::Color(.3, .6, .6), type);
      cur_particles.push_back( cur_p );
      cur_root.add_child(atom::Hierarchy::setup_particle(cur_p));
      ::npctransport_proto::Assignment_XYZ xyz = o_data.xyzs(j);
      core::XYZ p_xyz(cur_p);
      p_xyz.set_coordinates(algebra::Vector3D(xyz.x(), xyz.y(), xyz.z()));
      p_xyz.set_coordinates_are_optimized( !o_data.is_static() );
   }
  particles_[type] = cur_particles;
  // add interaction sites to particles of this type:
  if (o_data.interactions().value() > 0) {
    int nsites = o_data.interactions().value();
    IMP_LOG(WARNING, nsites << " sites added " << std::endl);
    set_sites(type, nsites, o_data.radius().value());
  }
  // add type-specific scoring scale factors
  get_scoring()->set_interaction_range_factor
    ( type, o_data.interaction_range_factor().value());
  get_scoring()->set_interaction_k_factor
    ( type, o_data.interaction_k_factor().value());
  set_obstacles_changed(true);
}

void SimulationData::add_interaction
( const ::npctransport_proto::Assignment_InteractionAssignment & idata)
{
  // TODO: if diffusers added later, need to update interaction statistics!
  core::ParticleType type0(idata.type0());
  core::ParticleType type1(idata.type1());
  if (idata.is_on().value()) {
    IMP_LOG(TERSE, "   Adding interacton between " << type0
            << " and " << type1  << std::endl);
  } else {
    return;
  }
  get_scoring()->add_interaction(idata);
  // add statistics about this interaction to interactions_stats_
  // between all diffusing particles
  Particles set0, set1; // TODO: turn to a real set?!
  IMP_CONTAINER_FOREACH // _1 is the particle index
    ( SingletonContainer,
      get_diffusers(),
      {
        if (IMP::core::Typed(get_model(), _1).get_type() == type0)
          {
            set0.push_back( get_model()->get_particle(_1) );
          }
        if (IMP::core::Typed(get_model(), _1).get_type() == type1)
          {
            set1.push_back( get_model()->get_particle(_1) );
          }
      }
      );
  double stats_contact_range = 1.5;  // TODO: make a param
  double stats_slack = 30; // TODO: make a param - this is only efficiency of
                           //       ClosePairContainer
  IMP_LOG(PROGRESS,
          "Interaction " << type0.get_string() << ", " << type1.get_string()
          << "  sizes: " << set0.size() << ", " << set1.size()
          << " statistics range: " << stats_contact_range
          << std::endl);
  if (set0.size() > 0 && set1.size() > 0) {
    InteractionType interaction_type = std::make_pair(type0, type1);
    IMP_NEW(BipartitePairsStatisticsOptimizerState, bpsos,
            (get_model(), interaction_type, set0, set1,
             stats_contact_range, stats_slack));
    bpsos->set_period(statistics_interval_frames_);
    interactions_stats_.push_back(bpsos);
  }

}

Model *SimulationData::get_model() {
  set_was_used(true);
  if (!m_) {
    m_ = new Model("NPC model %1%");
  }
  return m_;
}

// Note and beware: this method assumes that the hierarchy in the RMF file
// was constructed in the same way as the hierarchy within this SimulationData
// object. Use with great caution, otherwise unexpected results may arise
void SimulationData::initialize_positions_from_rmf(RMF::FileConstHandle f,
                                                   int frame) {
  f.set_current_frame(f.get_number_of_frames() - 1);
  // RMF::show_hierarchy_with_values(f.get_root_node());
  link_hierarchies_with_sites(f, get_root().get_children());
  if (frame == -1) {
    std::cout << "Loading from last frame of RMF file with "
              << f.get_number_of_frames() << " frames" << std::endl;
    IMP::rmf::load_frame(f, f.get_number_of_frames() - 1);
  } else {
    IMP::rmf::load_frame(f, frame);
  }
}

void SimulationData::link_rmf_file_handle(RMF::FileHandle fh) {
  // TODO: this should be updated if particles are updated / restraints are
  //       probably restraints linkage should be handled internally in
  //       Scoring.cpp (for better encapsulation
  IMP_LOG(TERSE, "Setting up dump" << std::endl);
  Scoring* s=get_scoring();
  add_hierarchies_with_sites(fh, atom::Hierarchy(get_root()).get_children());
  IMP::rmf::add_restraints(fh, RestraintsTemp(1, s->get_predr()));
  IMP::rmf::add_restraints(fh, s->get_chain_restraints());
  if (s->get_has_bounding_box()) {
    IMP::rmf::add_restraints(fh, RestraintsTemp(1, s->get_box_restraint()));
    IMP_NEW(display::BoundingBoxGeometry, bbg, (get_box()));
    IMP::rmf::add_static_geometries(fh, display::Geometries(1, bbg));
  }
  if (s->get_has_slab()) {
    IMP::rmf::add_restraints(fh, RestraintsTemp(1, s->get_slab_restraint()));
    IMP_NEW(SlabWireGeometry, slab_geometry,
            (slab_thickness_, tunnel_radius_, box_side_));
    IMP::rmf::add_static_geometries(fh, slab_geometry->get_components());
  }
}

rmf::SaveOptimizerState *SimulationData::get_rmf_sos_writer() {
  if (!rmf_sos_writer_) {
    if (get_rmf_file_name().empty()) return nullptr;
    RMF::FileHandle fh = RMF::create_rmf_file(get_rmf_file_name());
    link_rmf_file_handle(fh);
    IMP_NEW(rmf::SaveOptimizerState, los, (fh));
    std::cout << "Dump interval for RMF SaveOptimizerState set to "
              << dump_interval_frames_ << std::endl;
    los->set_period(dump_interval_frames_);
    rmf_sos_writer_ = los;
  }
  return rmf_sos_writer_;
}

void SimulationData::set_rmf_file_name(const std::string &new_name) {
  if (get_rmf_file_name() == new_name) return;  // nothing to do
  if (rmf_sos_writer_ && bd_) {
    bd_->remove_optimizer_state(rmf_sos_writer_);
    rmf_sos_writer_ = nullptr;
  }
  rmf_file_name_ = new_name;
  if (bd_) {
    bd_->add_optimizer_state(get_rmf_sos_writer());
  }
}

void SimulationData::dump_geometry() {
  IMP_OBJECT_LOG;
  base::Pointer<display::Writer> w = display::create_writer("dump.pym");
  IMP_NEW(TypedSitesGeometry, g, (get_diffusers()));
  for (base::map<core::ParticleType, algebra::Vector3Ds>::const_iterator it =
           sites_.begin();
       it != sites_.end(); ++it) {
    g->set_sites(it->first, it->second);
  }
  w->add_geometry(g);
  if (box_is_on_) {
    IMP_NEW(display::BoundingBoxGeometry, bbg, (get_box()));
    bbg->set_was_used(true);
    w->add_geometry(bbg);
  }
  if (slab_is_on_) {
    IMP_NEW(display::CylinderGeometry, cyl_geom, (get_cylinder()));
    w->add_geometry(cyl_geom);
    // IMP_NEW(SlabWireGeometry, slab_geometry,
    //         (slab_height_ /* h */ ,
    //          slab_radius_ /* r */,
    //          slab_width_ /* w */) );
    // w->add_geometry( slab_geometry );
  }
}

void SimulationData::reset_rmf() {
  if (get_rmf_sos_writer()) {
    get_rmf_sos_writer()->reset();
  }
}

Scoring * SimulationData::get_scoring()
{
  IMP_USAGE_CHECK(scoring_ != nullptr, "Null scoring");
  return scoring_;
}


atom::BrownianDynamics *SimulationData::get_bd() {
  set_was_used(true);
  if (!bd_) {
    bd_ = new atom::BrownianDynamics(m_);
    bd_->set_maximum_time_step(time_step_);
    bd_->set_maximum_move(range_ / 4);
    bd_->set_current_time(0.0);
    //#ifdef _OPENMP
    if (dump_interval_frames_ > 0 && !get_rmf_file_name().empty()) {
      bd_->add_optimizer_state(get_rmf_sos_writer());
    }
    //#endif
    bd_->set_scoring_function
      ( get_scoring()->get_scoring_function() );
    // add all kind of observers to the optimization:
    for (unsigned int i = 0; i < fgs_stats_.size(); ++i) {
      for (unsigned int j = 0; j < fgs_stats_[i].size(); ++j) {
        bd_->add_optimizer_states(fgs_stats_[i][j]);
      }
    }
    for (unsigned int i = 0; i < float_stats_.size(); ++i) {
      bd_->add_optimizer_states(float_stats_[i]);
    }
    if (slab_is_on_) {
      for (unsigned int i = 0; i < float_transport_stats_.size(); ++i) {
        bd_->add_optimizer_states(float_transport_stats_[i]);
        // associate each with this bd_, so it can update transport times
        for (unsigned int j = 0; j < float_transport_stats_[i].size(); j++) {
          float_transport_stats_[i][j]->set_owner(bd_);
        }
      }
    }
    for (unsigned int i = 0; i < chain_stats_.size(); ++i) {
      bd_->add_optimizer_states(chain_stats_[i]);
    }
    bd_->add_optimizer_states(interactions_stats_);
  }
  return bd_;
}

container::ListSingletonContainer *
SimulationData::get_diffusers(){
  // TODO: obstacles ?
  bool new_diffusers = !diffusers_ || diffusers_changed_ || obstacles_changed_;
  if (!diffusers_) {
    diffusers_ = new container::ListSingletonContainer(m_);
  }
  if (new_diffusers) {
    // Add all leaves of the hierarchy as diffusers
    ParticlesTemp leaves = get_as<ParticlesTemp>(atom::get_leaves(get_root()));
    IMP_LOG(TERSE, "Leaves are " << leaves << std::endl);
    diffusers_->set_particles(leaves);
    IMP_USAGE_CHECK(leaves.size() == diffusers_->get_indexes().size(),
                    "Set and get particles don't match");
    // mark the update
    set_diffusers_changed(false);
    set_obstacles_changed(false);
  }
  return diffusers_;
}

container::ListSingletonContainer *
SimulationData::get_optimizable_diffusers()
{
  if (!optimizable_diffusers_) {
    optimizable_diffusers_ = new container::ListSingletonContainer(m_);
  }
  IMP::ParticlesTemp optimizable_diffusers;
  IMP_CONTAINER_FOREACH
    (container::ListSingletonContainer,
     get_diffusers(),
     {
       Particle* p = m_->get_particle(_1);
       if(core::XYZ::get_is_setup(p)) {
         core::XYZ p_xyz(p);
         if(p_xyz.get_coordinates_are_optimized()){
           optimizable_diffusers.push_back ( p );
         }
       }
     }
     );
  optimizable_diffusers_->set_particles(optimizable_diffusers);
  return optimizable_diffusers_;
}


void SimulationData::set_sites(core::ParticleType t0, unsigned int n,
                               double r) {
  algebra::Sphere3D s(algebra::get_zero_vector_d<3>(), r);
  algebra::Vector3Ds sites = algebra::get_uniform_surface_cover(s, n);
  sites_[t0] = algebra::Vector3Ds(sites.begin(), sites.begin() + n);
}


void SimulationData::write_geometry(std::string out) {
  IMP_OBJECT_LOG;
  base::Pointer<display::Writer> w = display::create_writer(out);
  {
    IMP_NEW(TypedSitesGeometry, g, (get_diffusers()));
    for (base::map<core::ParticleType, algebra::Vector3Ds>::const_iterator it =
             sites_.begin();
         it != sites_.end(); ++it) {
      g->set_sites(it->first, it->second);
    }
    w->add_geometry(g);
  }
  Scoring * s = get_scoring();
  for (unsigned int i = 0; i < s->get_chain_restraints().size(); ++i) {
    IMP_NEW(display::RestraintGeometry, rsg, (s->get_chain_restraints()[i]));
    w->add_geometry(rsg);
  }
  if (s->get_has_bounding_box()) {
    IMP_NEW(display::RestraintGeometry, rsg, (s->get_box_restraint()));
    w->add_geometry(rsg);
  }
  if (s->get_has_slab()) {
    IMP_NEW(display::RestraintGeometry, slab_rsg, (s->get_slab_restraint()));
    w->add_geometry(slab_rsg);
  }
  {
    IMP_NEW(display::RestraintGeometry, prsg, (s->get_predr()));
    w->add_geometry(prsg);
  }
  if (box_is_on_) {
    IMP_NEW(display::BoundingBoxGeometry, bbg, (get_box()));
    w->add_geometry(bbg);
  }
  if (slab_is_on_) {
    // IMP_NEW(SlabWireGeometry, slab_geometry,
    //         (1000 /* h */ , 100 /* r */, 300 /* w */) );
    // w->add_geometry(slab_geometry);
    IMP_NEW(display::CylinderGeometry, cyl_geom, (get_cylinder()));
    w->add_geometry(cyl_geom);
  }
}

// TODO: turn into a template inline in unamed space?
/**
   updates (message).field() with a weighted average of its current
   value and new_value, giving weight n_old_frames, n_new_frames to each,
   respectively.
*/
#define UPDATE_AVG(n_frames, n_new_frames, message, field, new_value)     \
  (message).set_##field(static_cast<double>(n_frames *(message).field() + \
                                            n_new_frames *new_value) /    \
                        (n_frames + n_new_frames));

// number of site-site interactions between a and b
int SimulationData::get_number_of_interactions(Particle *a, Particle *b) const {
  if (core::get_distance(core::XYZR(a), core::XYZR(b)) > range_) return 0;
  const algebra::Vector3Ds &sa = sites_.find(core::Typed(a).get_type())->second;
  const algebra::Vector3Ds &sb = sites_.find(core::Typed(b).get_type())->second;
  int ct = 0;
  for (unsigned int i = 0; i < sa.size(); ++i) {
    for (unsigned int j = 0; j < sb.size(); ++j) {
      if (algebra::get_distance(sa[i], sb[j]) < range_) {
        ++ct;
      }
    }
  }
  return ct;
}

// see doc in .h file
boost::tuple<double, double, double, double>
SimulationData::get_interactions_and_interacting(
    const ParticlesTemp &kaps, const base::Vector<ParticlesTemps> &fgs) const {
  double interactions = 0, interacting = 0, bead_partners = 0,
         chain_partners = 0;
  for (unsigned int i = 0; i < kaps.size(); ++i) {
    bool kap_found = false;
    for (unsigned int j = 0; j < fgs.size(); ++j) {
      for (unsigned int k = 0; k < fgs[j].size(); ++k) {
        bool chain_found = false;
        for (unsigned int l = 0; l < fgs[j][k].size(); ++l) {
          int num = get_number_of_interactions(kaps[i], fgs[j][k][l]);
          if (num > 0) {
            interactions += num;
            ++bead_partners;
            if (!kap_found) ++interacting;
            kap_found = true;
            if (!chain_found) ++chain_partners;
            chain_found = true;
          }
        }
      }
    }
  }
  return boost::make_tuple(interactions, interacting, bead_partners,
                           chain_partners);
}

void SimulationData::reset_statistics_optimizer_states() {
  is_stats_reset_ = true;  // indicate to update_statistics()
  get_bd()->set_current_time(0.0);
  for (unsigned int i = 0; i < fgs_stats_.size(); ++i) {
    for (unsigned int j = 0; j < fgs_stats_[i].size(); ++j) {
      for (unsigned int k = 0; k < fgs_stats_[i][j].size(); ++k) {
        fgs_stats_[i][j][k]->reset();
      }
    }
  }
  for (unsigned int i = 0; i < float_stats_.size(); ++i) {
    for (unsigned int j = 0; j < float_stats_[i].size(); ++j) {
      float_stats_[i][j]->reset();
    }
  }
  if (slab_is_on_) {
    for (unsigned int i = 0; i < float_transport_stats_.size(); ++i) {
      for (unsigned int j = 0; j < float_transport_stats_[i].size(); ++j) {
        float_transport_stats_[i][j]->reset();
      }
    }
  }
  for (unsigned int i = 0; i < chain_stats_.size(); ++i) {
    for (unsigned int j = 0; j < chain_stats_[i].size(); ++j) {
      chain_stats_[i][j]->reset();
    }
  }
  for (unsigned int i = 0; i < interactions_stats_.size(); ++i) {
    interactions_stats_[i]->reset();
  }
}


// @param nf_new number of new frames accounted for in this statistics update
void SimulationData::update_statistics(const boost::timer &timer,
                                       unsigned int nf_new) {
  IMP_OBJECT_LOG;
  ::npctransport_proto::Output output;

  std::ifstream inf(output_file_name_.c_str(), std::ios::binary);
  output.ParseFromIstream(&inf);
  inf.close();
  ::npctransport_proto::Statistics &stats = *output.mutable_statistics();

  int nf = stats.number_of_frames();
  if (is_stats_reset_) {  // restart statistics from scratch
    // TODO: what's if multiple trials?
    nf = 0;
    is_stats_reset_ = false;
    for (int i = 0; i < stats.floaters().size(); i++) {
      (*stats.mutable_floaters(i)).clear_transport_time_points_ns();
    }
  }
  std::cout << "Updating statistics file " << output_file_name_
            << " that currently has " << nf << " frames, with " << nf_new
            << " additional frames" << std::endl;
  ParticlesTemp all;
  ParticlesTemps floaters;
  base::Vector<ParticlesTemps> fgs;
  for (unsigned int i = 0; i < type_of_fg.size(); ++i) {
    if (particles_.find(type_of_fg[i]) != particles_.end()) {
      fgs.push_back(ParticlesTemps());
      ParticlesTemp fgi_particles = particles_.find(type_of_fg[i])->second;
      for (unsigned int j = 0; j < fgi_particles.size(); ++j) {
        atom::Hierarchy h(fgi_particles[j]);
        ParticlesTemp chain = get_as<ParticlesTemp>(atom::get_leaves(h));
        all += chain;
        fgs.back().push_back(chain);
        IMP_ALWAYS_CHECK(stats.fgs_size() > static_cast<int>(i),
                         "Not enough fgs: " << stats.fgs_size() << " vs "
                                            << static_cast<int>(i),
                         ValueException);
#ifdef IMP_NPCTRANSPORT_USE_IMP_CGAL
        double volume = atom::get_volume(h);
        double radius_of_gyration = atom::get_radius_of_gyration(chain);
#else
        double volume = -1.;
        double radius_of_gyration = -1.;
#endif
        UPDATE_AVG(nf, nf_new, *stats.mutable_fgs(i), volume, volume);
        double length =
            core::get_distance(core::XYZ(chain[0]), core::XYZ(chain.back()));
        UPDATE_AVG(nf, nf_new, *stats.mutable_fgs(i), length, length);
        UPDATE_AVG(nf, nf_new, *stats.mutable_fgs(i), radius_of_gyration,
                   radius_of_gyration);
      }
    }
  }
  // correlation times and diffusion coefficient
  for (unsigned int i = 0; i < float_stats_.size(); ++i) {
    for (unsigned int j = 0; j < float_stats_[i].size(); ++j) {
      int cnf = (nf) * float_stats_[i].size() + j;
      IMP_ALWAYS_CHECK(stats.floaters_size() > static_cast<int>(i),
                       "Not enough floaters", ValueException);
      UPDATE_AVG(cnf, nf_new,  // TODO: is nf_new correct? I think so
                 *stats.mutable_floaters(i), diffusion_coefficient,
                 float_stats_[i][j]->get_diffusion_coefficient());
      UPDATE_AVG(cnf, nf_new,  // TODO: is nf_new correct? I think so
                 *stats.mutable_floaters(i), correlation_time,
                 float_stats_[i][j]->get_correlation_time());
      float_stats_[i][j]->reset();
    }
  }
  // update avg number of transports per particle for each type of floats:
  if (slab_is_on_) {
    for (unsigned int type_i = 0; type_i < float_transport_stats_.size();
         ++type_i) {
      unsigned int n_particles = float_transport_stats_[type_i].size();
      // collect individual transport times in an ordered set,
      // and add them to the statistics file:
      std::set<double> times_i(
          stats.floaters(type_i).transport_time_points_ns().begin(),
          stats.floaters(type_i).transport_time_points_ns().end());
      for (unsigned int j = 0; j < n_particles; ++j) {  // add new
        Floats const &new_times_ij = float_transport_stats_[type_i][j]
            ->get_transport_time_points_in_ns();
        for (unsigned int k = 0; k < new_times_ij.size(); k++) {
          times_i.insert(new_times_ij[k]);
        }
      }
      (*stats.mutable_floaters(type_i)).clear_transport_time_points_ns();
      for (std::set<double>::const_iterator it = times_i.begin();
           it != times_i.end(); it++) {
        (*stats.mutable_floaters(type_i)).add_transport_time_points_ns(*it);
      }
      // update avg too
      double avg_n_transports_type_i = times_i.size() * 1.0 / n_particles;
      (*stats.mutable_floaters(type_i))
          .set_avg_n_transports(avg_n_transports_type_i);
    }
  }
  for (unsigned int i = 0; i < fgs_stats_.size(); ++i) {
    for (unsigned int j = 0; j < fgs_stats_[i].size(); ++j) {
      for (unsigned int k = 0; k < fgs_stats_[i][j].size(); ++k) {
        unsigned int n = fgs_stats_[i].size() * fgs_stats_[i][j].size();
        unsigned int cnf = (nf) * n + j * fgs_stats_[i][j].size() + k;
        ;
        IMP_ALWAYS_CHECK(stats.fgs_size() > static_cast<int>(i),
                         "Not enough fgs", ValueException);
        UPDATE_AVG(cnf, nf_new, *stats.mutable_fgs(i),
                   particle_correlation_time,
                   fgs_stats_[i][j][k]->get_correlation_time());
        UPDATE_AVG(cnf, nf_new, *stats.mutable_fgs(i),
                   particle_diffusion_coefficient,
                   fgs_stats_[i][j][k]->get_diffusion_coefficient());
        fgs_stats_[i][j][k]->reset();
      }
    }
  }
  for (unsigned int i = 0; i < chain_stats_.size(); ++i) {
    for (unsigned int j = 0; j < chain_stats_[i].size(); ++j) {
      unsigned int n = chain_stats_[i].size();
      unsigned int cnf = (nf) * n + j;
      IMP_ALWAYS_CHECK(stats.fgs_size() > static_cast<int>(i), "Not enough fgs",
                       ValueException);
      /*UPDATE_AVG(cnf, nf_new,
             *stats.mutable_fgs(i), chain_correlation_time,
             chain_stats_[i][j]->get_correlation_time());*/
      UPDATE_AVG(cnf, nf_new, *stats.mutable_fgs(i),
                 chain_diffusion_coefficient,
                 chain_stats_[i][j]->get_diffusion_coefficient());
      Floats df = chain_stats_[i][j]->get_diffusion_coefficients();
      UPDATE_AVG(cnf, nf_new, *stats.mutable_fgs(i),
                 local_diffusion_coefficient,
                 std::accumulate(df.begin(), df.end(), 0.0) / df.size());
      chain_stats_[i][j]->reset();
    }
  }

  for (unsigned int i = 0; i < type_of_float.size(); ++i) {
    if (particles_.find(type_of_float[i]) != particles_.end()) {
      floaters.push_back(particles_.find(type_of_float[i])->second);
      all += floaters.back();
      double interactions, interacting, bead_partners, chain_partners;
      boost::tie(interactions, interacting, bead_partners, chain_partners) =
          get_interactions_and_interacting(floaters.back(), fgs);
      IMP_ALWAYS_CHECK(stats.floaters_size() > static_cast<int>(i),
                       "Not enough fgs", ValueException);
      UPDATE_AVG(nf, nf_new, *stats.mutable_floaters(i), interactions,
                 interactions / floaters.back().size());
      UPDATE_AVG(nf, nf_new, *stats.mutable_floaters(i), interacting,
                 interacting / floaters.back().size());
      UPDATE_AVG(nf, nf_new, *stats.mutable_floaters(i),
                 interaction_partner_chains, chain_partners / interacting);
      UPDATE_AVG(nf, nf_new, *stats.mutable_floaters(i),
                 interaction_partner_beads, bead_partners / interacting);
    }
  }

  // update statistics gathered on interaction rates
  for (unsigned int i = 0; i < interactions_stats_.size(); i++) {
    IMP_LOG(PROGRESS, "adding interaction statistics # " << i << std::endl);
    IMP_ALWAYS_CHECK(stats.interactions_size() > static_cast<int>(i),
                     "Not enough fgs", ValueException);
    ::npctransport_proto::Statistics_InteractionStats *pOutStats_i =
        stats.mutable_interactions(i);
    base::Pointer<BipartitePairsStatisticsOptimizerState> pInStats_i =
        interactions_stats_[i];
    // verify correct interaction type is stored
    InteractionType itype = pInStats_i->get_interaction_type();
    std::string s_type0 = itype.first.get_string();
    std::string s_type1 = itype.second.get_string();
    if (std::string(pOutStats_i->type0()) != s_type0 ||
        std::string(pOutStats_i->type1()) != s_type1) {
      IMP_THROW("Incompatible interaction types in pOutStats_i ["
                    << pOutStats_i->type0() << ", " << pOutStats_i->type1()
                    << "] and pInStats_i " << s_type0 << ", " << s_type1 << "]"
                    << std::endl,
                IMP::base::ValueException);
    }
    // save the rest of the interactions info
    Int n0 = pInStats_i->get_number_of_particles_1();
    Int n1 = pInStats_i->get_number_of_particles_2();
    Float avg_contacts_num = pInStats_i->get_average_number_of_contacts();

    UPDATE_AVG(nf, nf_new, *pOutStats_i, avg_contacts_per_particle0,
               avg_contacts_num / n0);
    UPDATE_AVG(nf, nf_new, *pOutStats_i, avg_contacts_per_particle1,
               avg_contacts_num / n1);
    UPDATE_AVG(nf, nf_new, *pOutStats_i, avg_pct_bound_particles0,
               pInStats_i->get_average_percentage_bound_particles_1());
    UPDATE_AVG(nf, nf_new, *pOutStats_i, avg_pct_bound_particles1,
               pInStats_i->get_average_percentage_bound_particles_2());
    interactions_stats_[i]->reset();
  }

  {
    double total_energy  = get_bd()->get_scoring_function()->evaluate(false);
    UPDATE_AVG(nf, nf_new, stats, energy_per_particle,  // TODO: reset?
               total_energy / all.size());
  }

  // Todo: define better what we want of timer
  // UPDATE_AVG(nf, nf_new, stats, seconds_per_iteration, timer.elapsed());
  stats.set_seconds_per_iteration(timer.elapsed());

  stats.set_number_of_frames(nf + nf_new);
  const double fs_in_ns = 1.0E+6;
  stats.set_bd_simulation_time_ns(const_cast<SimulationData *>(this)
                                      ->get_bd()->get_current_time() /
                                  fs_in_ns);

  ::npctransport_proto::Conformation *conformation =
      output.mutable_conformation();
  save_pb_conformation(get_diffusers(), sites_, conformation);

  // save rmf too
  {
    std::string buf;
    RMF::FileHandle fh = RMF::create_rmf_buffer(buf);
    const_cast<SimulationData *>(this)->link_rmf_file_handle(fh);
    IMP::rmf::save_frame(fh, 0);
    output.set_rmf_conformation(buf);
  }

  // dump to file
  std::ofstream outf(output_file_name_.c_str(), std::ios::binary);
  output.SerializeToOstream(&outf);
}
void SimulationData::set_interrupted(bool tf) {
  ::npctransport_proto::Output output;
  std::ifstream inf(output_file_name_.c_str(), std::ios::binary);
  output.ParseFromIstream(&inf);
  inf.close();
  ::npctransport_proto::Statistics &stats = *output.mutable_statistics();
  stats.set_interrupted(tf ? 1 : 0);
  std::ofstream outf(output_file_name_.c_str(), std::ios::binary);
  stats.SerializeToOstream(&outf);
}

display::Geometry *SimulationData::get_static_geometry() {
  if (!static_geom_) {
    IMP_NEW(display::BoundingBoxGeometry, bbg, (this->get_box()));
    static_geom_ = bbg;
  }
  return static_geom_;
}

algebra::BoundingBox3D SimulationData::get_box() const {
  return algebra::get_cube_d<3>(.5 * box_side_);
}

algebra::Cylinder3D SimulationData::get_cylinder() const {
  algebra::Vector3D pt(0, 0, slab_thickness_ / 2.0);
  algebra::Segment3D seg(pt, -pt);
  return algebra::Cylinder3D(seg, tunnel_radius_);
}
#undef GET_ASSIGNMENT
#undef GET_VALUE

IMPNPCTRANSPORT_END_NAMESPACE
