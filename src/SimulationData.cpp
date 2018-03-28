/**
 *  \file SimulationData.cpp
 *  \brief description.
 *
 *  Copyright 2007-2018 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/BrownianDynamicsTAMDWithSlabSupport.h>
#include <IMP/npctransport/FGChain.h>
#include <IMP/npctransport/ParticleFactory.h>
#include <IMP/npctransport/protobuf.h>
#include <IMP/npctransport/SimulationData.h>
#include <IMP/npctransport/SitesGeometry.h>
#include <IMP/npctransport/SlabWithCylindricalPoreGeometry.h>
#include <IMP/npctransport/SlabWithToroidalPoreGeometry.h>
#include <IMP/npctransport/SlabWithCylindricalPorePairScore.h>
#include <IMP/npctransport/SlabWithToroidalPorePairScore.h>
#include <IMP/npctransport/SlabWithCylindricalPore.h>
#include <IMP/npctransport/SlabWithToroidalPore.h>
#include <IMP/npctransport/ParticleFactory.h>
// #include <IMP/npctransport/particle_types.h>
#include <IMP/npctransport/protobuf.h>
#include <IMP/npctransport/internal/npctransport.pb.h>
#include <IMP/npctransport/enums.h>
#include <IMP/npctransport/io.h>
#include <IMP/npctransport/typedefs.h>
#include <IMP/npctransport/util.h>
#include <IMP/algebra/vector_generators.h>
#include <IMP/atom/distance.h>
#include <IMP/atom/Diffusion.h>
#include <IMP/atom/Mass.h>
#include <IMP/atom/estimates.h>
#include <IMP/atom/Selection.h>
#include <IMP/compiler_macros.h>
#include <IMP/log.h>
#include <IMP/internal/units.h>
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
#include <algorithm>
#include <iterator>
#include <numeric>
#include <set>

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/coded_stream.h>
#include <fcntl.h>
#if defined(_MSC_VER)
#include <io.h>
#else
#include <sys/types.h>
#include <sys/stat.h>
#endif

IMPNPCTRANSPORT_BEGIN_NAMESPACE
#define GET_ASSIGNMENT(name) name##_ = pb_assignment.name().value()
#define GET_ASSIGNMENT_DEF(name, default_value)                         \
  {                                                                     \
     if(!pb_assignment.has_##name() )                                    \
       {                                                                 \
         pb_mutable_assignment->mutable_##name()                         \
           ->set_value(default_value);                                   \
       }                                                                 \
     name##_ = pb_assignment.name().value();                             \
   }
#define GET_VALUE(name) name##_ = pb_assignment.name()
#define GET_VALUE_DEF(name, default_value)               \
  {                                                      \
    if(!pb_assignment.has_##name() )                     \
      { pb_mutable_assignment                            \
          ->set_##name(default_value); }                 \
    name##_ = pb_assignment.name();                      \
  }

//  #define GET_ASSIGNMENT_DEF(name, default_value) \
//   {                                             \
//     if(pb_assignment.has_##name() )              \
//       { name##_ = pb_assignment.name().value(); }       \
//     else                                                \
//       { name##_ = default_value; }                      \
//   }
// #define GET_VALUE_DEF(name, default_value)      \
//   {                                             \
//     if(pb_assignment.has_##name() )             \
//       { name##_ = pb_assignment.name(); }       \
//     else                                        \
//       { name##_ = default_value; }              \
//   }

#define SLAB_TYPE_STRING "slab"
namespace {
  bool is_slab_particle(Particle* p)
  {
    return SlabWithPore::get_is_setup(p);
  }
};

SimulationData::SimulationData(std::string prev_output_file, bool quick,
                               std::string rmf_file_name,
                               std::string new_output_file)
: Object("SimulationData%1%"),
  m_(nullptr),
  bd_(nullptr),
  scoring_(nullptr),
  statistics_(nullptr),
  root_(nullptr),
  slab_particle_(nullptr),
  rmf_file_name_(rmf_file_name),
  is_save_restraints_to_rmf_(true)
{
  if(new_output_file=="") {
    new_output_file = prev_output_file;
  }
  initialize(prev_output_file, new_output_file, quick);
}

void SimulationData::initialize(std::string prev_output_file,
                                std::string new_output_file,
                                bool quick) {
  m_ = new Model("NPC model %1%");
  ::npctransport_proto::Output pb_data;
  //  std::ifstream file(prev_output_file.c_str(), std::ios::binary);
  //  bool read = pb_data.ParseFromIstream(&file);
  bool read(false);
  int fd= IMP_C_OPEN(prev_output_file.c_str(),
                     IMP_C_OPEN_FLAG(O_RDONLY) | IMP_C_OPEN_BINARY);
  if(fd!=-1){
    google::protobuf::io::FileInputStream fis(fd);
    google::protobuf::io::CodedInputStream cis(&fis);
    cis.SetTotalBytesLimit(500000000,200000000);
    read= pb_data.ParseFromCodedStream(&cis);
    IMP_C_CLOSE(fd);
  }
  IMP_ALWAYS_CHECK(read,
                   "Unable to read data from protobuf" << prev_output_file,
                   IMP::IOException);
  pb_data.mutable_statistics(); // create it if not there
  const ::npctransport_proto::Assignment &
      pb_assignment= pb_data.assignment();
  ::npctransport_proto::Assignment*
      pb_mutable_assignment= pb_data.mutable_assignment();
  GET_VALUE_DEF(output_npctransport_version, 1.0); // old configuration files that don't have output version should be treated as version 1.0
  GET_ASSIGNMENT(box_side);
  GET_ASSIGNMENT(tunnel_radius);
  GET_ASSIGNMENT_DEF(tunnel_radius_k, -1.0);
  GET_ASSIGNMENT_DEF(pore_anchored_beads_k, -1.0);
  GET_ASSIGNMENT(slab_thickness);
  GET_ASSIGNMENT(slab_is_on);
  GET_ASSIGNMENT(box_is_on);
  GET_VALUE(number_of_trials);
  GET_VALUE(number_of_frames);
  GET_VALUE(dump_interval_frames);
  GET_ASSIGNMENT(angular_d_factor);
  GET_VALUE(range);
  GET_VALUE(statistics_interval_frames);
  GET_VALUE_DEF(output_statistics_interval_frames,10000);
  GET_ASSIGNMENT(statistics_fraction);
  GET_VALUE(time_step);
  GET_ASSIGNMENT_DEF(time_step_wave_factor, 0.0);
  GET_VALUE(maximum_number_of_minutes);
  GET_VALUE_DEF(fg_anchor_inflate_factor, 1.0);
  GET_VALUE_DEF(are_floaters_on_one_slab_side, false);
  GET_VALUE_DEF(is_exclude_floaters_from_slab_initially, false);
  GET_ASSIGNMENT_DEF(temperature_k, strip_units(IMP::internal::DEFAULT_TEMPERATURE));
  GET_VALUE_DEF(is_xyz_hist_stats, false)
  GET_VALUE_DEF(is_backbone_harmonic, false);
  GET_ASSIGNMENT_DEF(backbone_tau_ns, 1.0);
  initial_simulation_time_ns_ = 0.0; // default
  if (pb_data.has_statistics()) {
    if (pb_data.statistics().has_bd_simulation_time_ns()) {
      initial_simulation_time_ns_ =
        (pb_data.statistics().bd_simulation_time_ns() * FS_IN_NS);
    }
  }
  if (quick) {
    number_of_frames_ = 2;
    number_of_trials_ = 1;
  }
  // initialize scoring
  scoring_ = new Scoring(this, pb_assignment);
  statistics_ = new Statistics(this,
                               statistics_interval_frames_,
                               new_output_file);

  // create particles hierarchy
  root_ = new Particle(get_model());
  root_->add_attribute(get_simulation_data_key(), this);
  atom::Hierarchy hr = atom::Hierarchy::setup_particle(root_);
  root_->set_name("root");
  if(get_has_slab()){
    // TODO: at some point, obstacles should be the children of slab (cause in some sense, they refine it)
    atom::Hierarchy h_slab_particle=
      atom::Hierarchy::setup_particle( get_slab_particle() );
    // TODO: XYZR is fake - just to make it compatible with atom hierarchy (note pore radius is a different attribute than 'radius')
    core::XYZ::setup_particle( get_slab_particle(),
                               algebra::Vector3D(0,0,0) );
    core::XYZR d=
      core::XYZR::setup_particle( get_slab_particle(),
                                  0.0 );
    d.set_coordinates_are_optimized(false);
    atom::Mass::setup_particle( get_slab_particle(),
                                1.0 );
    core::Typed::setup_particle( get_slab_particle(),
                                 core::ParticleType(SLAB_TYPE_STRING));
    //    get_root().add_child( h_slab_particle ); // TODO: include it as child of root? This may have unforeseen consequences
  }
  for (int i = 0; i < pb_assignment.fgs_size(); ++i)
    {
      // verify type first
      IMP_ALWAYS_CHECK(pb_assignment.fgs(i).has_type(),
                       "old or corrupt assignment file"
                       << ", which lacks fg types" << std::endl,
                       ValueException);
      IMP_ALWAYS_CHECK(pb_assignment.fgs(i).type() != "",
                       "FG should've been assigned a valued type,"
                       " possiblu thru protobuf.h", ValueException);
      pb_mutable_assignment->mutable_fgs(i)->site_relative_distance(); // initiate it with default value for backward compatability
      create_fgs(pb_assignment.fgs(i));
    }
  for (int i = 0; i < pb_assignment.floaters_size(); ++i)
    {
      // verify type first
      IMP_ALWAYS_CHECK(pb_assignment.floaters(i).has_type(),
                       "old or corrupt assignment file"
                       << ", which lacks floater types"
                       << std::endl, ValueException);
      IMP_USAGE_CHECK( pb_assignment.floaters(i).type() != "",
                       "Floater should've been assigned a type"
                       << " thru protobuf.h");
      pb_mutable_assignment->mutable_floaters(i)->site_relative_distance(); // initiate it with default value for backward compatability
      create_floaters(pb_assignment.floaters(i),
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

  // Load from RMF conformation (or conformation) if they exist in protobuf
  if (pb_data.has_rmf_conformation())
    {
      IMP_LOG(VERBOSE, "Restarting from output file internal RMF conformation"
                << std::endl);
      RMF::BufferConstHandle buffer(pb_data.rmf_conformation());
      RMF::FileConstHandle fh =
        RMF::open_rmf_buffer_read_only(buffer);
      initialize_positions_from_rmf(fh, 0);
    } else if (pb_data.has_conformation())
    {
      IMP_LOG(VERBOSE, "Restarting from output file conformation" << std::endl);
      load_pb_conformation(pb_data.conformation(), get_beads(), sites_);
    }

  get_bd()->set_current_time( initial_simulation_time_ns_ );
  pb_mutable_assignment->add_imp_module_version
    ( IMP::get_module_version() );
  pb_mutable_assignment->add_npc_module_version
    ( IMP::npctransport::get_module_version() );
  std::ofstream outf(new_output_file.c_str(), std::ios::binary);
  pb_data.SerializeToOstream(&outf);
}


void
SimulationData::create_slab_particle()
{
  slab_particle_=
    new Particle(get_model(), SLAB_TYPE_STRING);
  if(get_is_slab_with_cylindrical_pore()){
    SlabWithCylindricalPore::setup_particle(slab_particle_,
                                            slab_thickness_,
                                            tunnel_radius_);
  } else{
    SlabWithToroidalPore::setup_particle(slab_particle_,
                                         slab_thickness_,
                                         tunnel_radius_,
                                         1.0); // TODO: add support for skewed toroids
  }
  SlabWithPore(slab_particle_)
    .set_pore_radius_is_optimized(get_is_pore_radius_dynamic());
}

/**
   Adds the FG Nup chains to the model hierarchy,
   based on the settings in data
*/
void SimulationData::create_fgs
( const ::npctransport_proto::Assignment_FGAssignment &fg_data)
{
  // Save type:
  IMP_ALWAYS_CHECK(fg_data.has_type(), "fg type missing in fg_data",
                   ValueException);
  core::ParticleType fg_chain_type(fg_data.type());
  IMP_ALWAYS_CHECK( fg_chain_types_.find(fg_chain_type) == fg_chain_types_.end(),
                    "Currently support only single insertion of each type,"
                    " can be fixed if needed in the future",
                    ValueException);
  fg_chain_types_.insert(fg_chain_type);


  // Make main root or all chains of this type:
  atom::Hierarchy chains_root =
    atom::Hierarchy::setup_particle(new IMP::Particle( get_model() ) );
  core::Typed::setup_particle(chains_root, fg_chain_type);
  chains_root->set_name( fg_chain_type.get_string() ); //+ " root" );
  get_root().add_child(chains_root);

  // Add n chains:
  Pointer<FGChain> representative_chain= nullptr;
  for (int j = 0; j < fg_data.number().value(); ++j) {
    Pointer<FGChain> chain= create_fg_chain
      (this, chains_root, fg_data, display::Color(.3, .3, .3));
    if(j==0){
      representative_chain= chain;
    }
    // Beads book keeping:
    beads_ += chain->get_beads();
    // Set chain j anchors by fg_data if specified:
    if (fg_data.anchor_coordinates_size() > j) {
      ::npctransport_proto::Assignment_XYZ xyz =
        fg_data.anchor_coordinates(j);
      core::XYZR d(chain->get_bead(0));
      d.set_coordinates(algebra::Vector3D(xyz.x(), xyz.y(), xyz.z()));
      d.set_radius(d.get_radius() * fg_anchor_inflate_factor_); // inflate
      if (get_is_pore_radius_dynamic()
          && pore_anchored_beads_k_>0) {
        get_scoring()->add_restrained_anchor_bead(d); // TODO: think if problems might arise when restarting - if current radius of pore is updated in output file
      }else{
        d.set_coordinates_are_optimized(false);
      }
    }
    // add stats on chain
    get_statistics()->add_fg_chain_stats( chain );
  } // for j

  // Register bead types from representative chain:
  Particles representative_beads= representative_chain->get_beads();
  ParticleTypeSet new_fg_bead_types;
  ParticleTypeSet existing_fg_bead_types;
  for(unsigned int k= 0; k<representative_beads.size(); k++){
    new_fg_bead_types.insert(IMP::core::Typed(representative_beads[k]).get_type());
  }
  std::set_intersection(new_fg_bead_types.begin(), new_fg_bead_types.end(),
                        fg_bead_types_.begin(), fg_bead_types_.end(),
                        std::inserter(existing_fg_bead_types,
                                      existing_fg_bead_types.end())
                        );;
  IMP_ALWAYS_CHECK(existing_fg_bead_types.size()==0,
                   "Cannot add bead types that are already associated with"
                   " a different chain type (check if chain type + suffix are equal,"
                   "e.g. Nu+p=N+up)",
                   ValueException);
  fg_bead_types_.insert(new_fg_bead_types.begin(),
                        new_fg_bead_types.end());


  // Add sites and general scoring scale factors information to scoring
  // for each bead type in this chain type:
  // NOTE: although they're all the same, use per bead type info so
  //       can be differentiated in the future if need be
  for( ParticleTypeSet::const_iterator it_type= new_fg_bead_types.begin();
       it_type != new_fg_bead_types.end(); it_type++){
    // Add sites for this type:
    if (fg_data.interactions().value() != 0) {
      int nsites = fg_data.interactions().value();
      set_sites(*it_type, nsites,
                fg_data.radius().value() * fg_data.site_relative_distance(),
                fg_data.site_radius());
    }
    get_scoring()->set_interaction_range_factor
      ( *it_type,
        fg_data.interaction_range_factor().value());
    get_scoring()->set_interaction_k_factor
      ( *it_type,
        fg_data.interaction_k_factor().value());
  }

}


/**
   Adds the 'floaters' (free diffusing particles) to the model hierarchy,
   based on the settings in data
*/
void SimulationData::create_floaters
( const ::npctransport_proto::Assignment_FloaterAssignment &f_data,
  display::Color color)
{
  core::ParticleType type(f_data.type());
  if (f_data.number().value() == 0) {
    return;
  }
  IMP_ALWAYS_CHECK( floater_types_.find(type) == floater_types_.end(),
                    "Currently support only single insertion of each type,"
                    " can be fixed in the future", ValueException);
  floater_types_.insert(type);

  // Main root for type:
  atom::Hierarchy cur_root =
    atom::Hierarchy::setup_particle(new Particle(get_model()));
  core::Typed::setup_particle(cur_root, type);
  cur_root->set_name(type.get_string());
  get_root().add_child(cur_root);
  // Populate floaters under root:
  IMP_NEW(ParticleFactory, pf,
          (this, f_data.radius().value(),
           f_data.d_factor().value(),
           angular_d_factor_,
           color, type) );
  if(f_data.beads_tail_n()>0){ // TODO: make beads a hierarchical structure - for now it's simplest this way
    for (int i=0; i < f_data.beads_tail_n(); i++) {
      IMP_NEW(ParticleFactory, pf,
              (this, f_data.beads_tail_radius(),
               f_data.d_factor().value(),
               angular_d_factor_,
               color, type) );
    }
  }
  for (int j = 0; j < f_data.number().value(); ++j)
    {
      Particle *cur_p = pf->create();
      cur_root.add_child(atom::Hierarchy::setup_particle(cur_p));
      beads_.push_back(cur_p); // book keeping
      get_statistics()->add_floater_stats(cur_p); // stats
    }
  // add interaction sites to particles of this type:
  if (f_data.interactions().value() != 0)
    {
      int nsites = f_data.interactions().value();
      IMP_ALWAYS_CHECK(f_data.site_coordinates_size() == nsites ||
                       f_data.site_coordinates_size() == 0,
                       "The number of sites in sites_coordinates_size() must equal "
                       "the specified interactions number",
                       IMP::ValueException);
      if(f_data.site_coordinates_size() == 0) {
        set_sites(type,
                  nsites,
                  f_data.radius().value() * f_data.site_relative_distance(),
                  f_data.site_radius());
      }else{
        sites_[type]= algebra::Sphere3Ds();
        for(int ii= 0; ii<nsites; ii++) {
          ::npctransport_proto::Assignment_XYZ xyz =
            f_data.site_coordinates(ii);
          algebra::Vector3D site_ii_center(xyz.x(), xyz.y(), xyz.z());
          algebra::Sphere3D site_ii(site_ii_center, f_data.site_radius());
          sites_[type].push_back(site_ii);
        }
      }
      IMP_LOG(WARNING, nsites << " sites added " << std::endl);
    }

  // add type-specific scoring scale factors
  get_scoring()->set_interaction_range_factor
    ( type, f_data.interaction_range_factor().value());
  get_scoring()->set_interaction_k_factor
      ( type, f_data.interaction_k_factor().value());

  // add z-biasing potential for fraction of particles
  if(f_data.has_k_z_bias()) if (f_data.k_z_bias().value() != 0)
    {
      double k = f_data.k_z_bias().value();
      unsigned int n_bias = f_data.number().value();
      if(f_data.has_k_z_bias_fraction()) {
        double kzbf = f_data.k_z_bias_fraction().value();
        if(kzbf <= 0.0 || kzbf > 1.0) kzbf = 1.0;
        n_bias = (unsigned int)
          std::ceil( kzbf * n_bias );
      }
      ParticlesTemp ps;
      for (unsigned int i = 0; i < n_bias ; i++)
        {
          ps.push_back( cur_root.get_child(i) );
        }
      IMP_LOG(VERBOSE, "Biasing " << n_bias << " floaters"
                << " of type " << type
                << std::endl);
      get_scoring()->add_z_bias_restraint(ps, k);
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
  if(o_data.xyzs_size() == 0){
    return;
  }
  obstacle_types_.insert(type);
  // Main root for type:
  atom::Hierarchy cur_root =
    atom::Hierarchy::setup_particle(new Particle(get_model()));
  core::Typed::setup_particle(cur_root, type);
  IMP_LOG(WARNING,  "   obstacle type " << type.get_string() << std::endl);
  cur_root->set_name(type.get_string());
  get_root().add_child(cur_root);
  // Populate hierarchy with obstacles:
  IMP_NEW(ParticleFactory, pf,
          ( this, o_data.radius().value(),
            o_data.d_factor().value(),
            angular_d_factor_,
            display::Color(.3, .6, .6), type) );
  for (int j = 0; j < o_data.xyzs_size(); ++j) {
      Particle *cur_p = pf->create();
      cur_root.add_child(atom::Hierarchy::setup_particle(cur_p));
      beads_.push_back(cur_p); // book keeping
      ::npctransport_proto::Assignment_XYZ xyz = o_data.xyzs(j);
      core::XYZ p_xyz(cur_p);
      p_xyz.set_coordinates(algebra::Vector3D(xyz.x(), xyz.y(), xyz.z()));
      p_xyz.set_coordinates_are_optimized( !o_data.is_static() );
   }
  // add interaction sites to particles of this type:
  if (o_data.interactions().value() != 0) {
    int nsites = o_data.interactions().value();
    IMP_LOG(WARNING, nsites << " sites added " << std::endl);
    set_sites(type, nsites,
              o_data.radius().value(),
              0.0);
  }
  // add type-specific scoring scale factors
  get_scoring()->set_interaction_range_factor
    ( type, o_data.interaction_range_factor().value());
  get_scoring()->set_interaction_k_factor
    ( type, o_data.interaction_k_factor().value());
}

void SimulationData::add_interaction
( const ::npctransport_proto::Assignment_InteractionAssignment & idata)
{
  // TODO: if beads added later, need to update interaction statistics!
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

/** returns true if bead type is of fg type
    (that is, particle was added within ::create_fgs()
     including bead specific suffix)
*/
bool SimulationData::get_is_fg_bead(ParticleIndex pi) const
{
  Model* m = get_model();
  if (core::Typed::get_is_setup(m,pi)) {
    core::ParticleType t = core::Typed(m,pi).get_type();
    if (fg_bead_types_.find( t ) != fg_bead_types_.end() ) {
      return true;
    }
  }
  return false;
 }

/** returns true if chain type is of fg type
    (that is, an fg chain with same root type
    was added via ::create_fgs())
*/
bool SimulationData::get_is_fg_chain(ParticleIndex pi) const
{
  Model* m = get_model();
  if (core::Typed::get_is_setup(m,pi)) {
    core::ParticleType t = core::Typed(m,pi).get_type();
    if (fg_chain_types_.find( t ) != fg_chain_types_.end() ) {
      return true;
    }
  }
  return false;
 }


//!  return all the fg hierarchies in the simulation data Object
//!  that stand for individual FG chains
atom::Hierarchies SimulationData::get_fg_chain_roots() const
    {
      atom::Hierarchies ret;
      atom::Hierarchies C = get_root().get_children();
      for(unsigned int i =0 ; i < C.size(); i++) {
        if( get_is_fg_chain( C[i] ) ){
           ret += C[i].get_children();
          }
      }
      return ret;
    }


/** return all the obstacle particles */
ParticlesTemp SimulationData::get_obstacle_particles() const
{
  typedef boost::unordered_set<core::ParticleType> ParticleTypeSet;

  ParticlesTemp ret;
  ParticleTypeSet const& o_types = get_obstacle_types();
  for(unsigned int i = 0; i < get_root().get_number_of_children(); i++) {
    atom::Hierarchy root_of_type = get_root().get_child(i);
    core::ParticleType type = core::Typed(root_of_type).get_type();
    if(o_types.find(type) != o_types.end()) {
      Particles o_leaves = get_leaves(root_of_type) ;
      ret.insert( ret.begin(), o_leaves.begin(), o_leaves.end() );
    }
  }
  return ret;
}


// Behold and beware: this method assumes that the hierarchy in the RMF file
// was constructed in the same way as the hierarchy within this SimulationData
// object. Use with great caution, otherwise unexpected results may arise
void SimulationData::initialize_positions_from_rmf(RMF::FileConstHandle f,
                                                   int frame) {
  f.set_current_frame( RMF::FrameID(f.get_number_of_frames() - 1) );
  // RMF::show_hierarchy_with_values(f.get_root_node());
  ParticlesTemp linked_particles  = get_root().get_children();
  if(get_output_npctransport_version()<2.5){
    // slab particles not supported in earlier versions
    std::remove_if(linked_particles.begin(),
                   linked_particles.end(),
                   is_slab_particle);
  }
  link_hierarchies_with_sites(f, linked_particles);
  //  link_hierarchies_with_sites(f, get_root().get_children());
  if (frame == -1) {
    IMP_LOG(VERBOSE, "Loading from last frame of RMF file with "
              << f.get_number_of_frames() << " frames" << std::endl);
    frame = f.get_number_of_frames() - 1;
  }
  IMP::rmf::load_frame(f, RMF::FrameID(frame) );
}

void SimulationData::link_rmf_file_handle(RMF::FileHandle fh,
                                          bool with_restraints) {
  // TODO: this should be updated if particles are updated / restraints are
  //       probably restraints linkage should be handled internally in
  //       Scoring.cpp (for better encapsulation
  IMP_LOG(TERSE, "Setting up dump" << std::endl);
  Scoring* s=get_scoring();
  //  std::cout << "Linking " << get_root().get_number_of_children() << " particles and their children" << std::endl;
  add_hierarchies_with_sites(fh, get_root().get_children());
  if(with_restraints) {
    IMP::rmf::add_restraints
      (fh, s->get_scoring_function_restraints(true));
  }
  if(with_restraints && false) {
    IMP::rmf::add_restraints(fh, RestraintsTemp
                             (1, s->get_predicates_pair_restraint() )
                             );
    IMP::rmf::add_restraints(fh, s->get_all_chain_restraints() );
    IMP::rmf::add_restraints(fh, s->get_z_bias_restraints() );
    IMP::rmf::add_restraints(fh, s->get_custom_restraints() );
  }
  if (s->get_has_bounding_box()) {
    if(with_restraints && false) {
      IMP::rmf::add_restraints
        ( fh, RestraintsTemp(1, s->get_bounding_box_restraint()) );
    }
    IMP_NEW(display::BoundingBoxGeometry, bbg, (get_box()));
    IMP::rmf::add_static_geometries(fh, display::Geometries(1, bbg));
  }
  if (get_has_slab()) {
    if(with_restraints && false) {
      IMP::rmf::add_restraints(fh, RestraintsTemp(1, s->get_slab_restraint()));
    }
    IMP::Pointer<display::Geometry> slab_geometry;
    if(get_is_slab_with_cylindrical_pore()){
      slab_geometry = new SlabWithCylindricalPoreWireGeometry
        (get_slab_thickness(), get_pore_radius(), box_side_);
    } else {
      slab_geometry = new SlabWithToroidalPoreWireGeometry
	(get_slab_thickness(), get_pore_radius(), box_side_);
    }
    IMP::rmf::add_static_geometries(fh, slab_geometry->get_components());
  }
}

rmf::SaveOptimizerState *SimulationData::get_rmf_sos_writer() {
  if (!rmf_sos_writer_) {
    if (get_rmf_file_name().empty()) return nullptr;
    RMF::FileHandle fh = RMF::create_rmf_file(get_rmf_file_name());
    link_rmf_file_handle(fh, is_save_restraints_to_rmf_);
    IMP_NEW(rmf::SaveOptimizerState, los, (get_model(), fh));
    IMP_LOG(VERBOSE, "Dump interval for RMF SaveOptimizerState set to "
              << dump_interval_frames_ << std::endl);
    los->set_period(dump_interval_frames_);
    rmf_sos_writer_ = los;
  }
  return rmf_sos_writer_;
}

void SimulationData::set_rmf_file(const std::string &new_name,
                                  bool is_save_restraints_to_rmf)
{
  if (get_rmf_file_name() == new_name &&
      is_save_restraints_to_rmf_ == is_save_restraints_to_rmf)
    {  // nothing to do
      return ;
    }
  // Delete old one
  if (rmf_sos_writer_ && bd_) {
    bd_->remove_optimizer_state(rmf_sos_writer_);
    rmf_sos_writer_ = nullptr;
  }
  // Update vars
  rmf_file_name_ = new_name;
  is_save_restraints_to_rmf_ = is_save_restraints_to_rmf;
  // Construct new
  if (bd_) {
    bd_->add_optimizer_state(get_rmf_sos_writer());
  }
}

void SimulationData::dump_geometry() {
  IMP_OBJECT_LOG;
  Pointer<display::Writer> w = display::create_writer("dump.pym");
  IMP_NEW(TypedSitesGeometry, g, (get_beads()));
  for (boost::unordered_map<core::ParticleType, algebra::Sphere3Ds>::const_iterator it =
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
    IMP::Pointer<display::Geometry> slab_geometry;
    if(get_is_slab_with_cylindrical_pore()) {
      slab_geometry= new SlabWithCylindricalPoreWireGeometry
	(get_slab_thickness(), get_pore_radius(), box_side_);
    } else {
      slab_geometry= new SlabWithToroidalPoreWireGeometry
	(get_slab_thickness(), get_pore_radius(), box_side_);
    }
    w->add_geometry( slab_geometry );
  }
}


/**
   Returns the root of all chains of type 'type'
*/
atom::Hierarchy
  SimulationData::get_root_of_type(core::ParticleType type) const
{
  atom::Hierarchies H = get_root().get_children();
  for(unsigned int i = 0; i < H.size(); i++) {
    if(core::Typed(H[i]).get_type() == type) {
      return H[i];
    }
  }
  IMP_THROW("particle type " << type << " not found",
            ValueException);
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

Scoring const * SimulationData::get_scoring() const
{
  IMP_USAGE_CHECK(scoring_ != nullptr, "Null scoring");
  return scoring_;
}

Statistics * SimulationData::get_statistics()
{
  IMP_USAGE_CHECK(statistics_ != nullptr, "Null stats");
  return statistics_;
}

Statistics const * SimulationData::get_statistics() const
{
  IMP_USAGE_CHECK(statistics_ != nullptr, "Null stats");
  return statistics_;
}

atom::BrownianDynamics
*SimulationData::get_bd
(bool recreate)
{
  if (!bd_ || recreate) {
    bd_ = new BrownianDynamicsTAMDWithSlabSupport(m_, "BD_tamd_slab%1%", time_step_wave_factor_);
    bd_->set_maximum_time_step(time_step_);
    bd_->set_maximum_move(range_ / 4);
    bd_->set_current_time(0.0);
    bd_->set_scoring_function
      ( get_scoring()->get_scoring_function() );
    bd_->set_temperature(temperature_k_);
    //#ifdef _OPENMP
    if (dump_interval_frames_ > 0 && !get_rmf_file_name().empty()) {
      bd_->add_optimizer_state(get_rmf_sos_writer());
    }
    //#endif
  }
  return bd_;
}

void SimulationData::activate_statistics(){
  if(! get_statistics()->get_is_activated()) {
    get_statistics()->add_optimizer_states( get_bd() );
  }
}

ParticlesTemp
SimulationData::get_optimizable_beads()
{
  IMP::ParticlesTemp ret;
  for(unsigned int i = 0; i < beads_.size(); i++)
    {
      if(core::XYZ::get_is_setup(beads_[i]))
        {
          core::XYZ xyz(beads_[i]);
          if(xyz.get_coordinates_are_optimized())
            {  ret.push_back ( beads_[i] );  }
        }
    }
  return ret;
}

ParticlesTemp
SimulationData::get_non_optimizable_beads()
{
  IMP::ParticlesTemp ret;
  for(unsigned int i = 0; i < beads_.size(); i++)
    {
      if(core::XYZ::get_is_setup(beads_[i]))
        {
          core::XYZ xyz(beads_[i]);
          if(!xyz.get_coordinates_are_optimized())
            {  ret.push_back ( beads_[i] );  }
        }
    }
  return ret;
}


void SimulationData::set_sites(core::ParticleType t, int n,
                               double r, double sr) {
  sites_[t] = algebra::Sphere3Ds();
  if(n>0 && r > 0.0) {
    algebra::Sphere3D s(algebra::get_zero_vector_d<3>(), r);
    algebra::Vector3Ds sites_pos = algebra::get_uniform_surface_cover(s, n);
    for(int i = 0; i < n; i++) {
      sites_[t].push_back(algebra::Sphere3D(sites_pos[i], sr));
    }
  }
  if(n == -1 ||  r == 0.0) { // centered at bead
    IMP_ALWAYS_CHECK(n<2, "Cannot set more than one site at bead center",
                     ValueException);
    algebra::Sphere3D zero(algebra::Vector3D(0,0,0), sr);
    sites_[t] = algebra::Sphere3Ds(1, zero);
  }
}


void SimulationData::write_geometry(std::string out) {
  IMP_OBJECT_LOG;
  Pointer<display::Writer> w = display::create_writer(out);
  {
    IMP_NEW(TypedSitesGeometry, g, (get_beads()));
    for (boost::unordered_map<core::ParticleType, algebra::Sphere3Ds>
           ::const_iterator it = sites_.begin(); it != sites_.end(); ++it){
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
  if (get_has_slab()) {
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
    IMP::Pointer<display::Geometry> slab_geometry;
    if(get_is_slab_with_cylindrical_pore()) {
      slab_geometry= new SlabWithCylindricalPoreWireGeometry
	(get_slab_thickness(), get_pore_radius(), box_side_);
    } else {
      slab_geometry= new SlabWithToroidalPoreWireGeometry
	(get_slab_thickness(), get_pore_radius(), box_side_);
    }
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
  IMP_USAGE_CHECK(get_has_bounding_box(), "no bounding box defined");
  return algebra::get_cube_d<3>(.5 * box_side_);
}

double SimulationData::get_box_size() const{
  IMP_USAGE_CHECK(get_has_bounding_box(), "no bounding box defined");
  return box_side_;
}

void SimulationData::set_box_size(double box_size) {
  box_side_= box_size;
}

algebra::Cylinder3D SimulationData::get_cylinder() const {
  IMP_USAGE_CHECK(get_is_slab_with_cylindrical_pore(),
                  "no slab with cylidrical pore defined");
  algebra::Vector3D pt(0, 0, get_slab_thickness() / 2.0);
  algebra::Segment3D seg(pt, -pt);
  return algebra::Cylinder3D(seg, get_pore_radius());
}

double
SimulationData::get_slab_thickness() const
{
  IMP_USAGE_CHECK(get_slab_particle() != nullptr && get_has_slab(),
                  "invalid slab - can't get thickness");
  SlabWithPore swp(get_slab_particle());
  return swp.get_thickness();
}

//! returns te current tunnel radius
double
SimulationData::get_tunnel_radius() const
{
  IMP_USAGE_CHECK(get_slab_particle() != nullptr && get_has_slab(),
                  "invalid slab - can't get pore radius");
  SlabWithPore swp(get_slab_particle());
  return swp.get_pore_radius();
}

#undef GET_ASSIGNMENT
#undef GET_ASSIGNMENT_DEF
#undef GET_VALUE
#undef GET_VALUE_DEF

IMPNPCTRANSPORT_END_NAMESPACE
