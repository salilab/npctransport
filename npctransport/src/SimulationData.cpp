/**
 *  \file simulation_data.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/SimulationData.h>
#include <IMP/npctransport/SitesPairScore.h>
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
#include <IMP/core/BoundingBox3DSingletonScore.h>
#include <IMP/atom/Selection.h>
#include <IMP/container/ConsecutivePairContainer.h>
#include <IMP/core/DistancePairScore.h>
#include <IMP/core/SphereDistancePairScore.h>
#include <IMP/core/HarmonicUpperBound.h>
#include <IMP/core/RestraintsScoringFunction.h>
#include <IMP/core/XYZR.h>
#include <IMP/core/generic.h>
#include <IMP/display/LogOptimizerState.h>
#include <IMP/display/PymolWriter.h>
#include <IMP/display/primitive_geometries.h>
#include <IMP/display/restraint_geometry.h>
#include <IMP/example/optimizing.h>
#include <IMP/example/randomizing.h>
#include <IMP/npctransport/rmf_links.h>
#include <RMF/FileHandle.h>
#include <RMF/FileConstHandle.h>
#include <IMP/rmf/atom_io.h>
#include <IMP/rmf/frames.h>
#include <numeric>

#include <set>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
#define GET_ASSIGNMENT(name)                    \
  name##_= data.assignment().name().value()
#define GET_VALUE(name)                         \
  name##_= data.assignment().name()

SimulationData::SimulationData(std::string output_file,
                               bool quick,
                               std::string rmf_file_name):
  Object("SimulationData%1%"),
  rmf_file_name_( rmf_file_name ),
  is_stats_reset_( false )
{
  initialize(output_file, quick);
}

namespace {

 void load_conformation(const ::npctransport_proto::Conformation &conformation,
                         SingletonContainer *diffusers,
                         compatibility::map<core::ParticleType,
                                            algebra::Vector3Ds> &sites) {
    IMP_CONTAINER_FOREACH(SingletonContainer, diffusers,
                          {
                            const ::npctransport_proto::Conformation_Particle& pcur
                              = conformation.particle(_2);
                            core::RigidBody rb(diffusers->get_model(),
                                               _1);
                            algebra::Vector3D translation(pcur.x(),
                                                          pcur.y(),
                                                          pcur.z());
                            algebra::Rotation3D rotation(algebra::Vector4D(pcur.r(),
                                                                           pcur.i(),
                                                                           pcur.j(),
                                                                           pcur.k()));

                            algebra::Transformation3D tr(rotation, translation);
                            rb.set_reference_frame(algebra::ReferenceFrame3D(tr));
                          });
    for ( int i=0; i< conformation.sites_size(); ++i) {
      const ::npctransport_proto::Conformation::Sites &cur
        = conformation.sites(i);
      core::ParticleType pt(cur.name());
      sites[pt].clear();
      for (int j=0; j< cur.coordinates_size(); ++j) {
        const ::npctransport_proto::Conformation::Coordinates
          &coords= cur.coordinates(j);
        sites[pt].push_back(algebra::Vector3D(coords.x(),
                                              coords.y(),
                                              coords.z()));
      }
    }
  }

}


void SimulationData::initialize(std::string output_file,
                               bool quick) {
  output_file_name_= output_file;
  ::npctransport_proto::Output data;
  data.mutable_statistics();
  std::ifstream file(output_file_name_.c_str(), std::ios::binary);
  bool read=data.ParseFromIstream(&file);
  if (!read) {
    IMP_THROW("Unable to read from protobuf", base::IOException);
  }
  GET_ASSIGNMENT(interaction_k);
  GET_ASSIGNMENT(interaction_range);
  GET_ASSIGNMENT(backbone_k);
  GET_ASSIGNMENT(box_side);
  GET_ASSIGNMENT(tunnel_radius);
  GET_ASSIGNMENT(slab_thickness);
  GET_ASSIGNMENT(slab_is_on);
  GET_ASSIGNMENT(box_is_on);
  GET_ASSIGNMENT(slack);
  GET_VALUE(number_of_trials);
  GET_VALUE(number_of_frames);
  GET_VALUE(dump_interval_frames);
  GET_ASSIGNMENT(nonspecific_k);
  GET_ASSIGNMENT(nonspecific_range);
  GET_ASSIGNMENT(angular_d_factor);
  GET_VALUE(statistics_interval_frames);
  GET_ASSIGNMENT(excluded_volume_k);
  GET_VALUE(range);
  GET_VALUE(time_step);
  GET_ASSIGNMENT(statistics_fraction);
  GET_VALUE(maximum_number_of_minutes);
  if(quick){
    number_of_frames_ = 2;
    number_of_trials_ = 1;
  }

  // create particles hierarchy
  root_= new Particle(get_m());
  root_->add_attribute(get_simulation_data_key(), this);
  atom::Hierarchy hr=atom::Hierarchy::setup_particle(root_);
  root_->set_name("root");
  for (int i=0; i< data.assignment().fgs_size(); ++i) {
    create_fgs(data.assignment().fgs(i), type_of_fg[i]);
  }
  for (int i=0; i< data.assignment().floaters_size(); ++i) {
    create_floaters(data.assignment().floaters(i), type_of_float[i],
                    display::get_display_color(i));
  }
  IMP_INTERNAL_CHECK(get_diffusers()->get_indexes().empty(),
                     "Particles should not be in diffusers yet");

  // Add all leaves of th hierarchy as the set of diffusers returned by
  // get_diffusers()
  ParticlesTemp leaves= get_as<ParticlesTemp>(atom::get_leaves(get_root()));
  IMP_LOG(TERSE, "Leaves are " << leaves << std::endl);
  get_diffusers()->set_particles(leaves);
  IMP_USAGE_CHECK(leaves.size() == get_diffusers()->get_indexes().size(),
                  "Set and get don't match");

  IMP_LOG(TERSE, "   SimulationData before adding interactions" <<std::endl);
  for (int i=0; i< data.assignment().interactions_size(); ++i) {
    const ::npctransport_proto::Assignment_InteractionAssignment&
        interaction_i = data.assignment().interactions(i);
    if (interaction_i.is_on().value()) {
      IMP_LOG(TERSE, "   Adding interacton " << i <<std::endl);
      add_interaction( interaction_i );
    }
  }

  // bounding box / slab constraints on diffusers
  if (box_is_on_) {
    create_bounding_box_restraint_on_diffusers();
  }
  if (slab_is_on_) {
    create_slab_restraint_on_diffusers();
  }

  if (data.has_rmf_conformation()) {
    RMF::FileConstHandle fh
      = RMF::open_rmf_buffer_read_only(data.rmf_conformation());
    initialize_positions_from_rmf(fh, 0);
  // load from output file
  } else if (data.has_conformation()) {
    std::cout << "Loading from output file " << std::endl;
     load_conformation(data.conformation(),
                       get_diffusers(),
                       sites_);
   }
  if (data.has_statistics()) {
    if (data.statistics().has_bd_simulation_time_ns()) {
      const double fs_in_ns = 1.0E+6;
      get_bd()->set_current_time
        ( data.statistics().bd_simulation_time_ns() * fs_in_ns );
    }
  }
}

/**
   Adds the 'floaters' (free diffusing particles) to the model hierarchy,
   based on the settings in data
*/
void SimulationData::
create_floaters(const  ::npctransport_proto::Assignment_FloaterAssignment&data,
                core::ParticleType type, display::Color color) {
  if (data.number().value() > 0) {
    // prepare statistics for this type of floaters:
    float_stats_.push_back
      (BodyStatisticsOptimizerStates());
    if( slab_is_on_ ){ // if has tunnel, create a list of particle stats
      float_transport_stats_.push_back
        (ParticleTransportStatisticsOptimizerStates());
    }
    // create a sub hierarchy with this type of floaters:
    atom::Hierarchy cur_root
      = atom::Hierarchy::setup_particle(new Particle(get_m()));
    std::cout << "   type " << type.get_string() << std::endl;
    cur_root->set_name(type.get_string());
    atom::Hierarchy(get_root()).add_child(cur_root);
    // populate hierarchy with particles:
    ParticlesTemp cur_particles;
    for (int j=0; j< data.number().value(); ++j) {
      double dc=data.d_factor().value();
      Particle* cur_p =
        create_particle(this, data.radius().value(),
                        angular_d_factor_, dc,
                        color,
                        type, type.get_string());
      cur_particles.push_back( cur_p );
      cur_root.add_child( atom::Hierarchy::setup_particle( cur_p ) );
      // add particle statistics:
      IMP_NEW( BodyStatisticsOptimizerState, bsos, ( cur_p ) );
      bsos->set_period( statistics_interval_frames_ );
      float_stats_.back().push_back( bsos );
      if( slab_is_on_ ) { // only if has tunnel
        IMP_NEW(ParticleTransportStatisticsOptimizerState, ptsos,
                ( cur_p,
                  -0.5 * slab_thickness_, // tunnel bottom
                  0.5 * slab_thickness_) // tunnel top
                );
        ptsos->set_period( statistics_interval_frames_ );
        float_transport_stats_.back().push_back( ptsos );
      }
      // add interaction sites to particles of this type:
      if (data.interactions().value() >0) {
        int nsites= data.interactions().value();
        std::cout << nsites << " sites added " << std::endl;
        set_sites(type, nsites, data.radius().value());
      }
    }
    particles_[type]=cur_particles;
    interaction_range_factors_[type]= data.interaction_range_factor().value();
    interaction_k_factors_[type]= data.interaction_k_factor().value();
  }
}

/**
   Adds the FG Nup chains to the model hierarchy,
   based on the settings in data
*/
void SimulationData::
create_fgs(const ::npctransport_proto::Assignment_FGAssignment&data,
           core::ParticleType type) {
  if (data.number().value() > 0) {
    fgs_stats_.push_back(base::Vector<BodyStatisticsOptimizerStates>());
    chain_stats_.push_back(ChainStatisticsOptimizerStates());
    ParticlesTemp cur_particles;
    atom::Hierarchy hi= atom::Hierarchy::setup_particle(new Particle(get_m()));
    hi->set_name(type.get_string());
    atom::Hierarchy(get_root()).add_child(hi);
    double rlf=data.rest_length_factor().value();
    backbone_scores_
      .push_back(new LinearWellPairScore(rlf*data.radius().value()*2.0,
                                         backbone_k_));
    for (int j=0; j< data.number().value(); ++j) {
      double dc=data.d_factor().value();
      atom::Hierarchy hc(create_chain(this,
                                      data.number_of_beads().value(),
                                      data.radius().value(),
                                      angular_d_factor_, dc,
                                      backbone_scores_.back(),
                                      display::Color(.3,.3,.3), type, "fg"));
      cur_particles.push_back(hc);
      ParticlesTemp chain=hc.get_children();
      chain_stats_.back()
        .push_back(new ChainStatisticsOptimizerState(chain));
      chain_stats_.back().back()->set_period(statistics_interval_frames_);
      fgs_stats_.back().push_back(BodyStatisticsOptimizerStates());
      for (unsigned int k=0; k < chain.size(); ++k) {
        fgs_stats_.back().back()
          .push_back(new BodyStatisticsOptimizerState(chain[k]));
        fgs_stats_.back().back()
          .back()->set_period(statistics_interval_frames_);
      }
      hi.add_child(atom::Hierarchy(cur_particles.back()));
      if (data.interactions().value() > 0) {
        int nsites= data.interactions().value();
        set_sites(type, nsites, data.radius().value());
      }
    }
    particles_[type]=cur_particles;
    interaction_range_factors_[type]= data.interaction_range_factor().value();
    interaction_k_factors_[type]= data.interaction_k_factor().value();
  }
}


/**
   Creates bounding volume restraints such as box restraint and slab restraints,
   based on the box_size_, slab_height_, slab_radius_, etc. class variables
*/
void SimulationData::create_bounding_box_restraint_on_diffusers()
{
  // Add bounding box restraint
  // TODO: what does backbone_spring_k_ has to do
  //       with bounding box constraint?
  IMP_NEW(core::HarmonicUpperBound, hub, (0, excluded_volume_k_));
  IMP_NEW(core::GenericBoundingBox3DSingletonScore<core::HarmonicUpperBound>,
          bbss,
          (hub.get(), get_box()));
  box_restraint_=container::create_restraint(bbss.get(),
                                             get_diffusers(),
                                             "bounding box");
}

void SimulationData::create_slab_restraint_on_diffusers() {
    // Add cylinder restraint
  IMP_NEW(SlabSingletonScore,
          slab_score,
          (slab_thickness_ /* h */, tunnel_radius_ /*r*/, excluded_volume_k_) );
  slab_restraint_=container::create_restraint(slab_score.get(),
                                              get_diffusers(),
                                              "bounding slab");
}

Model *SimulationData::get_m() {
  set_was_used(true);
  if (!m_) {
    m_= new Model("NPC model %1%");
  }
  return m_;
}

// Note and beware: this method assumes that the hierarchy in the RMF file
// was constructed in the same way as the hierarchy within this SimulationData
// object. Use with great caution, otherwise unexpected results may come
void
SimulationData::initialize_positions_from_rmf(RMF::FileConstHandle f, int frame)
{
  f.set_current_frame(f.get_number_of_frames() - 1);
  RMF::show_hierarchy_with_values(f.get_root_node());
  link_hierarchies_with_sites( f, get_root().get_children() );
  if (frame==-1) {
    std::cout << "Loading from last frame of RMF file with "
              << f.get_number_of_frames() << " frames" << std::endl;
    IMP::rmf::load_frame( f, f.get_number_of_frames() - 1 );
  } else {
    IMP::rmf::load_frame( f, frame );
  }
}

void
SimulationData::link_rmf_file_handle(RMF::FileHandle fh)
{
  IMP_LOG(TERSE, "Setting up dump" << std::endl);
  add_hierarchies_with_sites(fh, atom::Hierarchy(get_root()).get_children());
  IMP::rmf::add_restraints(fh, RestraintsTemp(1, get_predr()));
  IMP::rmf::add_restraints(fh, chain_restraints_);
  if(get_has_bounding_box()){
    IMP::rmf::add_restraints(fh, RestraintsTemp(1, box_restraint_));
    IMP_NEW(display::BoundingBoxGeometry, bbg, (get_box()));
    IMP::rmf::add_static_geometries(fh, display::Geometries(1, bbg));
  }
  if( get_has_slab() ){
    IMP::rmf::add_restraints(fh, RestraintsTemp(1, slab_restraint_));
    IMP_NEW(SlabWireGeometry, slab_geometry,
            ( slab_thickness_ , tunnel_radius_ , box_side_ ) );
    IMP::rmf::add_static_geometries
      (fh, slab_geometry->get_components() );
  }
}


rmf::SaveOptimizerState *
SimulationData::get_rmf_sos_writer()
{
  if (!rmf_sos_writer_) {
    IMP_ALWAYS_CHECK( !get_rmf_file_name().empty(),
                      "RMF file name was not set", IMP::base::ValueException );
    RMF::FileHandle fh=RMF::create_rmf_file(get_rmf_file_name());
    link_rmf_file_handle(fh);
    IMP_NEW(rmf::SaveOptimizerState, los, (fh));
    std::cout << "Dump interval for RMF SaveOptimizerState set to "
              << dump_interval_frames_ << std::endl;
    los->set_period( dump_interval_frames_ );
    rmf_sos_writer_ = los;
  }
  return rmf_sos_writer_;
}

void SimulationData::set_rmf_file_name(const std::string& new_name)
{
  if( get_rmf_file_name() == new_name )
    return;// nothing to do
  if( rmf_sos_writer_ && bd_ ) {
    bd_->remove_optimizer_state( rmf_sos_writer_ );
    rmf_sos_writer_ = nullptr;
  }
  rmf_file_name_ = new_name;
  if( bd_ ) {
    bd_->add_optimizer_state( get_rmf_sos_writer() );
  }
}

void SimulationData::dump_geometry()
{
  IMP_OBJECT_LOG;
  Pointer<display::Writer> w= display::create_writer("dump.pym");
  IMP_NEW(TypedSitesGeometry, g, (get_diffusers()));
  for (compatibility::map<core::ParticleType,
         algebra::Vector3Ds>::const_iterator it= sites_.begin();
       it != sites_.end(); ++it) {
    g->set_sites(it->first, it->second);
  }
  w->add_geometry(g);
  if(get_has_bounding_box()){
    IMP_NEW(display::BoundingBoxGeometry, bbg, (get_box()));
    bbg->set_was_used(true);
    w->add_geometry(bbg);
  }
  if( get_has_slab() ){
    IMP_NEW(display::CylinderGeometry, cyl_geom, (get_cylinder()) );
    w->add_geometry( cyl_geom);
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

atom::BrownianDynamics *SimulationData::get_bd() {
  set_was_used(true);
  if (!bd_) {
    bd_=new atom::BrownianDynamics(m_);
    bd_->set_maximum_time_step(time_step_);
    bd_->set_maximum_move(range_/4);
    bd_->set_current_time( 0.0 );
    //#ifdef _OPENMP
    if (dump_interval_frames_ > 0 && !get_rmf_file_name().empty()) {
      bd_->add_optimizer_state(get_rmf_sos_writer());
    }
    //#endif
    // set up the restraints for the BD simulation:
    RestraintsTemp rs= chain_restraints_;
    if(get_has_bounding_box()) rs.push_back(box_restraint_);
    if( get_has_slab() ) rs.push_back(slab_restraint_);
    rs.push_back( this->get_predr() );
    IMP_NEW(core::RestraintsScoringFunction, rsf, (rs));
    bd_->set_scoring_function(rsf);
    // add all kind of observers to the optimization:
    for (unsigned int i=0; i< fgs_stats_.size(); ++i) {
      for (unsigned int j=0; j< fgs_stats_[i].size(); ++j) {
        bd_->add_optimizer_states(fgs_stats_[i][j]);
      }
    }
    for (unsigned int i=0; i< float_stats_.size(); ++i) {
      bd_->add_optimizer_states(float_stats_[i]);
    }
    if( slab_is_on_ ){
      for (unsigned int i=0; i< float_transport_stats_.size(); ++i) {
        bd_->add_optimizer_states(float_transport_stats_[i]);
        // associate each with this bd_, so it can update transport times
        for(unsigned int j=0 ; j < float_transport_stats_[i].size() ; j++) {
          float_transport_stats_[i][j]->set_owner(bd_);
        }
      }
    }
    for (unsigned int i=0; i< chain_stats_.size(); ++i) {
      bd_->add_optimizer_states(chain_stats_[i]);
    }
    bd_->add_optimizer_states(interactions_stats_);
  }
  return bd_;
}



container::ListSingletonContainer *
SimulationData::get_diffusers() {
  if (!diffusers_) {
    diffusers_ = new container::ListSingletonContainer(m_);
  }
  return diffusers_;
}

// a close pair container for all diffusers
container::ClosePairContainer* SimulationData::get_cpc() {
  if (!cpc_) {
    cpc_= new container::ClosePairContainer(get_diffusers(),
                                            range_,
                                            slack_);
    cpc_->add_pair_filter(new container::ExclusiveConsecutivePairFilter());
  }
  return cpc_;
}

container::PredicatePairsRestraint* SimulationData::get_predr() {
  if (!predr_) {
    // set linear repulsion upon penetration between all close pairs
    // returned by get_cpc(), with different scores for interactions
    // between particles of different (ordered) types
    IMP_NEW(core::OrderedTypePairPredicate, otpp, ());
    otpp_=otpp;
    IMP_NEW(container::PredicatePairsRestraint, ppr, (otpp, get_cpc()));
    predr_=ppr;
    IMP_NEW(LinearSoftSpherePairScore, ssps, (excluded_volume_k_));
    ppr->set_unknown_score(ssps.get() );
  }
  return predr_;
}

void SimulationData::set_sites(core::ParticleType t0,
                               unsigned int n, double r) {
  algebra::Sphere3D s(algebra::get_zero_vector_d<3>(), r);
  algebra::Vector3Ds sites
      = algebra::get_uniform_surface_cover(s, n);
  sites_[t0]=algebra::Vector3Ds(sites.begin(), sites.begin()+n);
}

/**
   add a SitesPairScore restraint that applies to particles of
   types t0 and t1 to the PredicatePairsRestraint object returned by
   get_predr().

   A SitesPairScore interaction means site-specific
   attractive forces between bidning sites on each particle,
   and non-specific attraction and repulsion (upon penetraion)
   between the particles themselves.
*/
void
SimulationData::add_interaction
(const ::npctransport_proto::Assignment_InteractionAssignment& idata )
{
  // extract interaction params
  core::ParticleType type0(idata.type0());
  core::ParticleType type1(idata.type1());
  double base_k=interaction_k_;
  if (idata.has_interaction_k()) {
    base_k= idata.interaction_k().value();
  }
  // no particles so drop it
  if (interaction_k_factors_.find(type0)== interaction_k_factors_.end()
      || interaction_k_factors_.find(type1)== interaction_k_factors_.end()) {
    return;
  }
  double interaction_k= base_k
    * interaction_k_factors_.find(type0)->second // TODO: validate type exists
    * interaction_k_factors_.find(type1)->second;
  double base_range=interaction_range_;
  if (idata.has_interaction_range()) {
    base_range= idata.interaction_range().value();
  }
  double interaction_range= base_range
    * interaction_range_factors_.find(type0)->second
    * interaction_range_factors_.find(type1)->second;

  std::cout << "creating interaction "
            << idata.type0() << "," << idata.type1()
            << " effective_k = " << interaction_k
            << ", effective_range = " << interaction_range
            << ", nonspecific k = " << nonspecific_k_
            << ", nonspecific range = " << nonspecific_range_
            << ", excluded volume k = " << excluded_volume_k_
            << std::endl;

  // create interaction
  container::PredicatePairsRestraint *ppr= get_predr();
  // add the interaction restraint both for (t0,t1) and (t1,t0)
  {
    // TODO: repulsion was also added in get_predr - do we double count here?
    core::ParticleTypes ts;
    ts.push_back(type0);
    ts.push_back(type1);
    int interaction_id= otpp_->get_value(ts);
    set_sites_score(interaction_range, // site-specific
                    interaction_k,
                    nonspecific_range_, // non-specific = entire particle
                    nonspecific_k_,
                    excluded_volume_k_,
                    sites_[type0], sites_[type1],
                    interaction_id , ppr);
  }
  {
    core::ParticleTypes ts;
    ts.push_back(type1);
    ts.push_back(type0);
    int interaction_id= otpp_->get_value(ts);
    set_sites_score(interaction_range, // site-specific
                    interaction_k,
                    nonspecific_range_, // non-specific = entire particle
                    nonspecific_k_,
                    excluded_volume_k_,
                    sites_[type0], sites_[type1],
                    interaction_id , ppr);
  }
  // add statistics about this interaction to interactions_stats_
  // between all diffusing particles
  ParticlesTemp set0, set1;
  for(unsigned int i = 0; i < diffusers_->get_particles().size(); i++) {
    if(IMP::core::Typed(diffusers_->get_particles()[i]).get_type() == type0) {
      set0.push_back( diffusers_->get_particles()[i] );
    }
    if(IMP::core::Typed(diffusers_->get_particles()[i]).get_type() == type1) {
      set1.push_back( diffusers_->get_particles()[i] );
    }
  }
  double stats_contact_range = 1.5; // TODO: make a param
  IMP_LOG( PROGRESS,
           "Interaction "
           << type0.get_string() << ", "
           << type1.get_string()
           << "  sizes: " << set0.size() << ", " << set1.size()
           << " statistics range: " << stats_contact_range << std::endl );
  if(set0.size() > 0 && set1.size() > 0) {
    InteractionType interaction_type = std::make_pair(type0,type1);
    IMP_NEW( BipartitePairsStatisticsOptimizerState,
             bpsos ,
             ( get_m(), interaction_type,
               set0, set1,
               stats_contact_range ) );
    bpsos->set_period(statistics_interval_frames_);
    interactions_stats_.push_back (bpsos);
  }
}


void SimulationData::write_geometry(std::string out) {
  IMP_OBJECT_LOG;
  Pointer<display::Writer> w= display::create_writer(out);
  {
    IMP_NEW(TypedSitesGeometry, g, (get_diffusers()));
    for (compatibility::map<core::ParticleType,
           algebra::Vector3Ds>::const_iterator it= sites_.begin();
         it != sites_.end(); ++it) {
      g->set_sites(it->first, it->second);
    }
    w->add_geometry(g);
  }
  for (unsigned int i=0; i< chain_restraints_.size(); ++i) {
    IMP_NEW(display::RestraintGeometry, rsg,
            (chain_restraints_[i]));
    w->add_geometry(rsg);
  }
  if(get_has_bounding_box()){
    IMP_NEW(display::RestraintGeometry, rsg,
            (box_restraint_));
    w->add_geometry(rsg);
  }
  if( get_has_slab() ){
    IMP_NEW(display::RestraintGeometry, slab_rsg,
            (slab_restraint_));
    w->add_geometry(slab_rsg);
  }
  {
    IMP_NEW(display::RestraintGeometry, prsg,
            (predr_));
    w->add_geometry(prsg);
  }
  if(get_has_bounding_box()){
    IMP_NEW(display::BoundingBoxGeometry, bbg, (get_box()));
    w->add_geometry(bbg);
  }
  if( get_has_slab() ){
    // IMP_NEW(SlabWireGeometry, slab_geometry,
    //         (1000 /* h */ , 100 /* r */, 300 /* w */) );
    // w->add_geometry(slab_geometry);
    IMP_NEW(display::CylinderGeometry, cyl_geom, (get_cylinder()) );
    w->add_geometry( cyl_geom);
  }

}

// TODO: turn into a template inline in unamed space?
/**
   updates (message).field() with a weighted average of its current
   value and new_value, giving weight n_old_frames, n_new_frames to each,
   respectively.
*/
#define UPDATE_AVG(n_frames, n_new_frames, message, field, new_value)   \
  (message).set_##field                                                 \
  ( static_cast<double>                                                 \
  (n_frames*(message).field() + n_new_frames * new_value) /             \
    (n_frames + n_new_frames) );


int SimulationData::get_number_of_interactions(Particle *a, Particle *b) const{
  if (core::get_distance(core::XYZR(a), core::XYZR(b)) > range_) return 0;
  const algebra::Vector3Ds &sa= sites_.find(core::Typed(a).get_type())->second;
  const algebra::Vector3Ds &sb= sites_.find(core::Typed(b).get_type())->second;
  int ct=0;
  for (unsigned int i=0; i< sa.size(); ++i) {
    for (unsigned int j=0; j< sb.size(); ++j) {
      if (algebra::get_distance(sa[i], sb[j]) < range_) {
        ++ct;
      }
    }
  }
  return ct;
}

boost::tuple<double,double,double, double>
SimulationData
::get_interactions_and_interacting(const ParticlesTemp &kaps,
                                   const base::Vector<ParticlesTemps> &fgs)
  const {
  double interactions=0, interacting=0, bead_partners=0, chain_partners=0;
  for (unsigned int i=0; i < kaps.size(); ++i) {
    bool found=false;
    for (unsigned int j=0; j< fgs.size(); ++j) {
      for (unsigned int k=0; k < fgs[j].size(); ++k) {
        bool chain_found=false;
        for (unsigned int l=0; l< fgs[j][k].size(); ++l) {
          int num= get_number_of_interactions(kaps[i], fgs[j][k][l]);
          if (num>0) {
            interactions+=num;
            ++bead_partners;
            if (!found) ++interacting;
            found=true;
            if (!chain_found) ++chain_partners;
            chain_found=true;
          }
        }
      }
    }
  }
  return boost::make_tuple(interactions, interacting, bead_partners,
                           chain_partners);
}


void SimulationData::reset_statistics_optimizer_states() {
  is_stats_reset_ = true; // indicate to update_statistics()
  get_bd()->set_current_time( 0.0 );
  for (unsigned int i=0; i< fgs_stats_.size(); ++i) {
    for (unsigned int j=0; j < fgs_stats_[i].size(); ++j) {
      for (unsigned int k=0; k < fgs_stats_[i][j].size(); ++k) {
        fgs_stats_[i][j][k]->reset();
      }
    }
  }
  for (unsigned int i=0; i< float_stats_.size(); ++i) {
    for (unsigned int j=0; j < float_stats_[i].size(); ++j) {
      float_stats_[i][j]->reset();
    }
  }
  if( slab_is_on_ ) {
    for (unsigned int i=0; i< float_transport_stats_.size(); ++i) {
      for (unsigned int j=0; j < float_transport_stats_[i].size(); ++j) {
        float_transport_stats_[i][j]->reset();
      }
    }
  }
  for (unsigned int i=0; i< chain_stats_.size(); ++i) {
    for (unsigned int j=0; j < chain_stats_[i].size(); ++j) {
      chain_stats_[i][j]->reset();
    }
  }
  for (unsigned int i=0; i< interactions_stats_.size(); ++i) {
    interactions_stats_[i]->reset();
  }
}

namespace {
  void save_conformation(SingletonContainer *diffusers,
                         const compatibility::map<core::ParticleType,
                                                  algebra::Vector3Ds> &sites,
                         ::npctransport_proto::Conformation *conformation) {
    conformation->clear_sites();
    conformation->clear_particle();
    IMP_CONTAINER_FOREACH(SingletonContainer, diffusers,
                          {
                            ::npctransport_proto::Conformation_Particle* pcur
                              = conformation->add_particle();
                            core::RigidBody rb(diffusers->get_model(),
                                               _1);
                            algebra::Transformation3D tr
                              = rb.get_reference_frame().get_transformation_to();
                            pcur->set_x(tr.get_translation()[0]);
                            pcur->set_y(tr.get_translation()[1]);
                            pcur->set_z(tr.get_translation()[2]);
                            pcur->set_r(tr.get_rotation().get_quaternion()[0]);
                            pcur->set_i(tr.get_rotation().get_quaternion()[1]);
                            pcur->set_j(tr.get_rotation().get_quaternion()[2]);
                            pcur->set_k(tr.get_rotation().get_quaternion()[3]);
                          });
    typedef compatibility::map<core::ParticleType,
                               algebra::Vector3Ds>  M;
    for ( M::const_iterator it = sites.begin(); it != sites.end() ;it++) {
      ::npctransport_proto::Conformation::Sites *cur
        = conformation->add_sites();
      cur->set_name(it->first.get_string());
      algebra::Vector3Ds coords = it->second;
      for (algebra::Vector3Ds::const_iterator coord = coords.begin();
             coord != coords.end(); coord++) {
        ::npctransport_proto::Conformation::Coordinates
          *out_coords= cur->add_coordinates();
        out_coords->set_x((*coord)[0]);
        out_coords->set_y((*coord)[1]);
        out_coords->set_z((*coord)[2]);
      }
    }
  }
}


// @param nf_new number of new frames accounted for in this statistics update
void
SimulationData::update_statistics
(const boost::timer &timer, unsigned int nf_new) const
{
  IMP_OBJECT_LOG;
  ::npctransport_proto::Output output;

  std::ifstream inf(output_file_name_.c_str(), std::ios::binary);
  output.ParseFromIstream(&inf);
  inf.close();
  ::npctransport_proto::Statistics &stats=*output.mutable_statistics();

  int nf= stats.number_of_frames();
  if(is_stats_reset_){ // restart statistics from scratch
    // TODO: what's if multiple trials?
    nf = 0;
    is_stats_reset_ = false;
    for (int i =0; i < stats.floaters().size(); i++){
      (*stats.mutable_floaters(i)).clear_transport_time_points_ns();
    }
  }
  std::cout << "Updating statistics file " << output_file_name_
            << " that currently has " << nf << " frames, with "
            << nf_new << " additional frames" << std::endl;
  ParticlesTemp all;
  ParticlesTemps floaters;
  base::Vector<ParticlesTemps> fgs;
  for (unsigned int i=0; i<type_of_fg.size(); ++i) {
    if (particles_.find(type_of_fg[i]) != particles_.end()) {
      fgs.push_back(ParticlesTemps());
      ParticlesTemp fgi_particles= particles_.find(type_of_fg[i])->second;
      for (unsigned int j=0; j< fgi_particles.size(); ++j) {
        atom::Hierarchy h(fgi_particles[j]);
        ParticlesTemp chain
            = get_as<ParticlesTemp>(atom::get_leaves(h));
        all+=chain;
        fgs.back().push_back(chain);
        IMP_ALWAYS_CHECK(stats.fgs_size() > static_cast<int>(i),
                         "Not enough fgs: " << stats.fgs_size()
                         << " vs " << static_cast<int>(i),
                         ValueException);
#ifdef IMP_NPCTRANSPORT_USE_IMP_CGAL
        double volume= atom::get_volume(h);
        double radius_of_gyration= atom::get_radius_of_gyration(chain);
#else
        double volume= -1.;
        double radius_of_gyration= -1.;
#endif
        UPDATE_AVG(nf, nf_new, *stats.mutable_fgs(i), volume, volume);
        double length= core::get_distance(core::XYZ(chain[0]),
                                          core::XYZ(chain.back()));
        UPDATE_AVG(nf, nf_new, *stats.mutable_fgs(i), length, length);
        UPDATE_AVG(nf, nf_new, *stats.mutable_fgs(i), radius_of_gyration,
               radius_of_gyration);
      }
    }
  }
  // correlation times and diffusion coefficient
  for (unsigned int i=0; i<float_stats_.size(); ++i) {
    for (unsigned int j=0; j<float_stats_[i].size(); ++j) {
      int cnf=(nf)*float_stats_[i].size()+j;
      IMP_ALWAYS_CHECK(stats.floaters_size() >  static_cast<int>(i),
                       "Not enough floaters",
                       ValueException);
      UPDATE_AVG(cnf, nf_new, // TODO: is nf_new correct? I think so
             *stats.mutable_floaters(i), diffusion_coefficient,
             float_stats_[i][j]->get_diffusion_coefficient());
      UPDATE_AVG(cnf, nf_new, // TODO: is nf_new correct? I think so
             *stats.mutable_floaters(i), correlation_time,
             float_stats_[i][j]->get_correlation_time());
      float_stats_[i][j]->reset();
    }
  }
  // update avg number of transports per particle for each type of floats:
  if( slab_is_on_ ) {
    for (unsigned int type_i=0; type_i<float_transport_stats_.size(); ++type_i) {
      unsigned int n_particles = float_transport_stats_[type_i].size();
      // collect individual transport times in an ordered set,
      // and add them to the statistics file:
      std::set<double>
          times_i( stats.floaters(type_i).transport_time_points_ns().begin(),
                   stats.floaters(type_i).transport_time_points_ns().end());
      for(unsigned int j = 0; j < n_particles ; ++j) { // add new
        Floats const& new_times_ij =
          float_transport_stats_[type_i][j]->get_transport_time_points_in_ns();
        for(unsigned int k = 0; k < new_times_ij.size(); k++) {
          times_i.insert( new_times_ij[k] );
        }
      }
      (*stats.mutable_floaters(type_i)).clear_transport_time_points_ns();
      for(std::set<double>::const_iterator it = times_i.begin() ;
          it != times_i.end() ; it++) {
        (*stats.mutable_floaters(type_i))
          .add_transport_time_points_ns( *it );
      }
      // update avg too
      double avg_n_transports_type_i =
        times_i.size() * 1.0 / n_particles;
      (*stats.mutable_floaters(type_i)).set_avg_n_transports
        ( avg_n_transports_type_i );
    }
  }
  for (unsigned int i=0; i<fgs_stats_.size(); ++i) {
    for (unsigned int j=0; j<fgs_stats_[i].size(); ++j) {
      for (unsigned int k=0; k< fgs_stats_[i][j].size(); ++k) {
        unsigned int n=fgs_stats_[i].size()*fgs_stats_[i][j].size();
        unsigned int cnf=(nf)*n+j*fgs_stats_[i][j].size()+k;;
        IMP_ALWAYS_CHECK(stats.fgs_size() >  static_cast<int>(i), "Not enough fgs",
                         ValueException);
        UPDATE_AVG(cnf, nf_new,
                   *stats.mutable_fgs(i),particle_correlation_time,
                   fgs_stats_[i][j][k]->get_correlation_time());
        UPDATE_AVG(cnf, nf_new,
                   *stats.mutable_fgs(i),particle_diffusion_coefficient,
                   fgs_stats_[i][j][k]->get_diffusion_coefficient());
        fgs_stats_[i][j][k]->reset();
      }
    }
  }
  for (unsigned int i=0; i<chain_stats_.size(); ++i) {
    for (unsigned int j=0; j<chain_stats_[i].size(); ++j) {
      unsigned int n=chain_stats_[i].size();
      unsigned int cnf=(nf)*n+j;
      IMP_ALWAYS_CHECK(stats.fgs_size() > static_cast<int>(i), "Not enough fgs",
                       ValueException);
      /*UPDATE_AVG(cnf, nf_new,
             *stats.mutable_fgs(i), chain_correlation_time,
             chain_stats_[i][j]->get_correlation_time());*/
      UPDATE_AVG(cnf, nf_new,
             *stats.mutable_fgs(i), chain_diffusion_coefficient,
             chain_stats_[i][j]->get_diffusion_coefficient());
      Floats df=chain_stats_[i][j]->get_diffusion_coefficients();
      UPDATE_AVG(cnf, nf_new,
             *stats.mutable_fgs(i), local_diffusion_coefficient,
             std::accumulate(df.begin(), df.end(), 0.0)/df.size());
      chain_stats_[i][j]->reset();
    }
  }

  for (unsigned int i=0; i<type_of_float.size(); ++i) {
    if (particles_.find(type_of_float[i]) != particles_.end()) {
      floaters.push_back(particles_.find(type_of_float[i])->second);
      all+=floaters.back();
      double interactions, interacting, bead_partners, chain_partners;
      boost::tie(interactions, interacting, bead_partners, chain_partners)
          =get_interactions_and_interacting(floaters.back(), fgs);
      IMP_ALWAYS_CHECK(stats.floaters_size() >  static_cast<int>(i),
                       "Not enough fgs", ValueException);
      UPDATE_AVG(nf, nf_new, *stats.mutable_floaters(i), interactions,
             interactions/floaters.back().size());
      UPDATE_AVG(nf, nf_new, *stats.mutable_floaters(i), interacting,
             interacting/floaters.back().size());
      UPDATE_AVG(nf, nf_new, *stats.mutable_floaters(i),
                 interaction_partner_chains,
             chain_partners/interacting);
      UPDATE_AVG(nf, nf_new, *stats.mutable_floaters(i),
                 interaction_partner_beads,
             bead_partners/interacting);
    }
  }

  // update statistics gathered on interaction rates
  for(unsigned int i = 0; i < interactions_stats_.size(); i++){
    IMP_LOG( PROGRESS, "adding interaction statistics # " << i << std::endl );
    IMP_ALWAYS_CHECK(stats.interactions_size() > static_cast<int>(i),
                         "Not enough fgs", ValueException);
    ::npctransport_proto::Statistics_InteractionStats*
      pOutStats_i = stats.mutable_interactions(i);
    IMP::Pointer<BipartitePairsStatisticsOptimizerState>
      pInStats_i = interactions_stats_[i];
    // verify correct interaction type is stored
    InteractionType itype = pInStats_i->get_interaction_type();
    std::string s_type0 = itype.first.get_string();
    std::string s_type1 = itype.second.get_string();
    if( std::string(pOutStats_i->type0()) != s_type0
        || std::string(pOutStats_i->type1()) != s_type1) {
      IMP_THROW("Incompatible interaction types in pOutStats_i ["
                << pOutStats_i->type0() << ", " << pOutStats_i->type1()
                << "] and pInStats_i " << s_type0 << ", " << s_type1
                << "]" << std::endl ,
                IMP::base::ValueException );
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
           pInStats_i->get_average_percentage_bound_particles_1() );
    UPDATE_AVG(nf, nf_new, *pOutStats_i, avg_pct_bound_particles1,
           pInStats_i->get_average_percentage_bound_particles_2() );
    interactions_stats_[i]->reset();
  }

  UPDATE_AVG(nf, nf_new, stats, energy_per_particle, // TODO: reset?
             get_m()->evaluate(false)/all.size());

  // Todo: define better what we want of timer
  //UPDATE_AVG(nf, nf_new, stats, seconds_per_iteration, timer.elapsed());
  stats.set_seconds_per_iteration( timer.elapsed() );

  stats.set_number_of_frames( nf + nf_new );
  const double fs_in_ns = 1.0E+6;
  stats.set_bd_simulation_time_ns
    ( const_cast<SimulationData*>(this)->
      get_bd()->get_current_time() / fs_in_ns );

  ::npctransport_proto::Conformation *conformation
      = output.mutable_conformation();
  save_conformation(diffusers_,
                    sites_, conformation);

  // save rmf too
  {
    std::string buf;
    RMF::FileHandle fh= RMF::create_rmf_buffer(buf);
    const_cast<SimulationData*>(this)->
      link_rmf_file_handle(fh);
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
  ::npctransport_proto::Statistics &stats= *output.mutable_statistics();
  stats.set_interrupted(tf?1:0);
  std::ofstream outf(output_file_name_.c_str(), std::ios::binary);
  stats.SerializeToOstream(&outf);
}


display::Geometry* SimulationData::get_static_geometry() {
  if (!static_geom_) {
    IMP_NEW(display::BoundingBoxGeometry, bbg, (this->get_box()));
    static_geom_=bbg;
  }
  return static_geom_;
}


algebra::Cylinder3D SimulationData::get_cylinder() const {
  algebra::Vector3D pt(0,0, slab_thickness_/2.0);
  algebra::Segment3D seg(pt, -pt);
  return algebra::Cylinder3D(seg, tunnel_radius_);
}


IMPNPCTRANSPORT_END_NAMESPACE
