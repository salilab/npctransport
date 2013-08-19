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
#include <IMP/npctransport/util.h>
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
#define GET_VALUE_DEF(name, default_value)          \
  {                                             \
  if(pb_assignment.has_##name() )               \
    { name##_ = pb_assignment.name(); }         \
  else                                          \
    { name##_ = default_value; }                \
  }

SimulationData::SimulationData(std::string output_file, bool quick,
                               std::string rmf_file_name)
: Object("SimulationData%1%"),
  m_(nullptr),
  bd_(nullptr),
  scoring_(nullptr),
  statistics_(nullptr),
  diffusers_changed_(false),
  obstacles_changed_(false),
  optimizable_diffusers_(nullptr),
  diffusers_(nullptr),
  rmf_file_name_(rmf_file_name)
{
  initialize(output_file, quick);
}

void SimulationData::initialize(std::string output_file, bool quick) {
  m_ = new Model("NPC model %1%");
  ::npctransport_proto::Output pb_data_;
  std::ifstream file(output_file.c_str(), std::ios::binary);
  bool read = pb_data_.ParseFromIstream(&file);
  if (!read) {
    IMP_THROW("Unable to read data from protobuf" << output_file,
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
  GET_ASSIGNMENT(statistics_fraction);
  GET_VALUE(time_step);
  GET_VALUE(maximum_number_of_minutes);
  GET_VALUE_DEF(fg_anchor_inflate_factor, 1.0);
  GET_VALUE_DEF(are_floaters_on_one_slab_side, false);
  if (quick) {
    number_of_frames_ = 2;
    number_of_trials_ = 1;
  }
  // initialize scoring
  scoring_ = new Scoring(this, pb_assignment);
  statistics_ = new Statistics(this, statistics_interval_frames_, output_file);

  // create particles hierarchy
  root_ = new Particle(get_model());
  root_->add_attribute(get_simulation_data_key(), this);
  atom::Hierarchy hr = atom::Hierarchy::setup_particle(root_);
  root_->set_name("root");
  for (int i = 0; i < pb_assignment.fgs_size(); ++i) {
    IMP_USAGE_CHECK( pb_assignment.fgs(i).has_type() &&
                     pb_assignment.fgs(i).type() != "",
                     "FG should've been assigned a type thru protobuf.h");
    create_fgs(pb_assignment.fgs(i));
  }
  for (int i = 0; i < pb_assignment.floaters_size(); ++i) {
    IMP_USAGE_CHECK( pb_assignment.floaters(i).has_type() &&
                     pb_assignment.floaters(i).type() != "",
                     "Floater should've been assigned a type thru protobuf.h");
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
  IMP_USAGE_CHECK(fg_data.has_type(), "It is assumed that fg data has"
                  "a type by now");
  IMP_LOG(PROGRESS, "creating FG of type '" << fg_data.type() << "'");
  core:: ParticleType type(fg_data.type());
  if (fg_data.number().value() > 0) {
    IMP_ALWAYS_CHECK( fg_types_.find(type) == fg_types_.end(),
                      "Currently support only single insertion of each type,"
                      " can be easily fixed in the future", base::ValueException);
    fg_types_.insert(type);
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
        core::XYZR d(chain_beads[0]);
        d.set_coordinates(algebra::Vector3D(xyz.x(), xyz.y(), xyz.z()));
        d.set_radius(d.get_radius() * fg_anchor_inflate_factor_); // inflate
        d.set_coordinates_are_optimized(false);
      }
      get_statistics()->add_fg_chain_stats( chain_beads ); // stats
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
  if (f_data.number().value() == 0) {
    return;
  }
  IMP_ALWAYS_CHECK( floater_types_.find(type) == floater_types_.end(),
                    "Currently support only single insertion of each type,"
                    " can be fixed in the future", base::ValueException);
  floater_types_.insert(type);
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
    get_statistics()->add_floater_stats(cur_p); // stats
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


/**
   Adds the 'obstacles' (possibly static e.g. nups that make the pore)
   to the model hierarchy, based on the settings in data
*/
void SimulationData::create_obstacles
( const ::npctransport_proto::Assignment_ObstacleAssignment &o_data)
{
  core::ParticleType type = core::ParticleType(o_data.type());
  if(o_data.xyzs_size() == 0){
    return;
  }
  obstacle_types_.insert(type);
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
  get_statistics()->add_interaction_stats(type0, type1);
}

Model *SimulationData::get_model() {
  set_was_used(true);
  IMP_USAGE_CHECK(m_, "model not initialized in get_model()");
  return m_;
}

/** returns true if particle type is of fg type
    (that is, particle was added within ::create_fgs()
*/
bool SimulationData::get_is_fg(ParticleIndex pi) const
{
  Model* m = get_model();
  if (core::Typed::get_is_setup(m,pi)) {
    core::ParticleType t = core::Typed(m,pi).get_type();
    if (fg_types_.find( t ) != fg_types_.end() ) {
      return true;
    }
  }
  return false;
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
  // I. return root itself if the type of its first direct child is fg
  ParticleIndex pi_child0 = root.get_child(0).get_particle_index();
  if (get_is_fg(pi_child0)) {
      return atom::Hierarchies(1, root);
  }
  // II. If not returned yet, recurse on all of root's children
  for (unsigned int i = 0; i < root.get_number_of_children(); ++i) {
    ret += get_fg_chains(root.get_child(i));
  }
  return ret;
}

/** return all the obstacle particles */
ParticlesTemp SimulationData::get_obstacle_particles() const
{
  ParticlesTemp ret;
  typedef base::set<core::ParticleType> ParticleTypeSet;
  ParticleTypeSet const& o_types = get_obstacle_types();
  for(ParticleTypeSet::const_iterator
        ti = o_types.begin(); ti != o_types.end(); ti++)
    {
      atom::Hierarchies hs = get_particles_of_type(*ti);
      ret += ParticlesTemp( atom::get_leaves( hs ) );
    }
  return ret;
}




// Note and beware: this method assumes that the hierarchy in the RMF file
// was constructed in the same way as the hierarchy within this SimulationData
// object. Use with great caution, otherwise unexpected results may arise
void SimulationData::initialize_positions_from_rmf(RMF::FileConstHandle f,
                                                   int frame) {
  f.set_current_frame( RMF::FrameID(f.get_number_of_frames() - 1) );
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
  IMP::rmf::add_restraints( fh, RestraintsTemp
                            (1, s->get_predicates_pair_restraint() )
                           );
  IMP::rmf::add_restraints(fh, s->get_all_chain_restraints() );
  if (s->get_has_bounding_box()) {
    IMP::rmf::add_restraints
      ( fh, RestraintsTemp(1, s->get_bounding_box_restraint()) );
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
    IMP_NEW(rmf::SaveOptimizerState, los, (get_model(), fh));
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
    //IMP_NEW(display::CylinderGeometry, cyl_geom, (get_cylinder()));
    //w->add_geometry(cyl_geom);
     IMP_NEW(SlabWireGeometry, slab_geometry,
             (slab_thickness_ /* h */ ,
              tunnel_radius_ /* r */,
              box_side_ /* w */) );
     w->add_geometry( slab_geometry );
  }
}

void SimulationData::reset_rmf() {
  if (get_rmf_sos_writer()) {
    get_rmf_sos_writer()->reset();
  }
}

//     temporarily suspend output to RMF file
void SimulationData::switch_suspend_rmf(bool suspend)
{
  if(!rmf_sos_writer_)
    return;
  if(suspend){
    get_bd(false)->remove_optimizer_state(rmf_sos_writer_);
  } else {
    get_bd(false)->add_optimizer_state(rmf_sos_writer_);
  }
}


Scoring * SimulationData::get_scoring()
{
  IMP_USAGE_CHECK(scoring_ != nullptr, "Null scoring");
  return scoring_;
}

Statistics * SimulationData::get_statistics()
{
  IMP_USAGE_CHECK(statistics_ != nullptr, "Null stats");
  return statistics_;
}


atom::BrownianDynamics *SimulationData::get_bd(bool recreate) {
  if (!bd_ || recreate) {
    bd_ = new atom::BrownianDynamics(m_);
    bd_->set_maximum_time_step(time_step_);
    bd_->set_maximum_move(range_ / 4);
    bd_->set_current_time(0.0);
    bd_->set_scoring_function
      ( get_scoring()->get_scoring_function() );
    //#ifdef _OPENMP
    if (dump_interval_frames_ > 0 && !get_rmf_file_name().empty()) {
      bd_->add_optimizer_state(get_rmf_sos_writer());
    }
    //#endif
    get_statistics()->add_optimizer_states( bd_ );
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
  optimizable_diffusers_->set_particles(get_diffusers()->get_particles());
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
  for (unsigned int i = 0; i < s->get_all_chain_restraints().size(); ++i) {
    IMP_NEW(display::RestraintGeometry, rsg,
            (s->get_all_chain_restraints()[i]));
    w->add_geometry(rsg);
  }
  if (s->get_has_bounding_box()) {
    IMP_NEW(display::RestraintGeometry, rsg, (s->get_bounding_box_restraint()));
    w->add_geometry(rsg);
  }
  if (s->get_has_slab()) {
    IMP_NEW(display::RestraintGeometry, slab_rsg, (s->get_slab_restraint()));
    w->add_geometry(slab_rsg);
  }
  {
    IMP_NEW(display::RestraintGeometry, prsg,
            (s->get_predicates_pair_restraint()));
    w->add_geometry(prsg);
  }
  if (box_is_on_) {
    IMP_NEW(display::BoundingBoxGeometry, bbg, (get_box()));
    w->add_geometry(bbg);
  }
  if (slab_is_on_) {
    IMP_NEW(SlabWireGeometry, slab_geometry,
            (slab_thickness_, tunnel_radius_, box_side_));
    w->add_geometry(slab_geometry);
    //IMP_NEW(display::CylinderGeometry, cyl_geom, (get_cylinder()));
    //w->add_geometry(cyl_geom);
  }
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
#undef GET_VALUE_DEF

IMPNPCTRANSPORT_END_NAMESPACE