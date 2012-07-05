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
#ifdef IMP_NPC_GOOGLE
#include "third_party/npc/module/data/npctransport.pb.h"
#else
#include "npctransport.pb.h"
#endif
#include <IMP/npctransport/creating_particles.h>
#include <IMP/npctransport/io.h>
#include <IMP/npctransport/typedefs.h>
#include <IMP/algebra/vector_generators.h>
#include <IMP/atom/estimates.h>
#include <IMP/atom/distance.h>
#include <IMP/atom/Diffusion.h>
#include <IMP/rmf/atom_io.h>
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


IMPNPCTRANSPORT_BEGIN_NAMESPACE
#define GET_ASSIGNMENT(name)                    \
name##_= data.name().value()
#define GET_VALUE(name)                         \
  name##_= data.name()

SimulationData::SimulationData(std::string assignment_file,
                               std::string statistics_file,
                               bool quick,
                               std::string rmf_file_name):
  Object("SimulationData%1%"),
  rmf_file_name_( rmf_file_name )
{
  statistics_file_name_= statistics_file;
  ::npctransport::Assignment data;
  std::ifstream file(assignment_file.c_str(), std::ios::binary);
  data.ParseFromIstream(&file);
  GET_ASSIGNMENT(interaction_k);
  GET_ASSIGNMENT(interaction_range);
  GET_ASSIGNMENT(backbone_k);
  GET_ASSIGNMENT(box_side);
  GET_ASSIGNMENT(tunnel_radius);
  GET_ASSIGNMENT(slab_thickness);
  GET_ASSIGNMENT(slack);
  GET_VALUE(number_of_trials);
  GET_VALUE(number_of_frames);
  GET_VALUE(dump_interval);
  GET_ASSIGNMENT(nonspecific_k);
  GET_ASSIGNMENT(nonspecific_range);
  GET_ASSIGNMENT(angular_d_factor);
  GET_VALUE(statistics_interval);
  GET_ASSIGNMENT(excluded_volume_k);
  GET_VALUE(range);
  GET_VALUE(time_step);
  if(quick){
    number_of_frames_ = 2;
    number_of_trials_ = 1;
  }
  first_stats_=true;

  // create particles hierarchy
  root_= new Particle(get_m());
  root_->add_attribute(get_simulation_data_key(), this);
  atom::Hierarchy hr=atom::Hierarchy::setup_particle(root_);
  root_->set_name("root");
  bonds_= new container::PairContainerSet(get_m());
  for (int i=0; i< data.fgs_size(); ++i) {
    create_fgs(data.fgs(i), type_of_fg[i]);
  }
  for (int i=0; i< data.floaters_size(); ++i) {
    create_floaters(data.floaters(i), type_of_float[i],
                    display::get_display_color(i));
  }
  for (int i=0; i< data.interactions_size(); ++i) {
    const ::npctransport::Assignment_InteractionAssignment&
      interaction_i = data.interactions(i);
    if (interaction_i.is_on().value()) {
      add_interaction( interaction_i );
    }
  }
  ParticlesTemp leaves= get_as<ParticlesTemp>(atom::get_leaves(get_root()));
  get_diffusers()->set_particles(leaves);
  if (data.box_is_on().value()) {
    create_bounding_box_restraint();
  }
  if (data.slab_is_on().value()) {
    create_slab_restraint();
 }
}

/**
   Adds the 'floaters' (free diffusing particles) to the model hierarchy,
   based on the settings in data
*/
void SimulationData::
create_floaters(const  ::npctransport::Assignment_FloaterAssignment&data,
                core::ParticleType type, display::Color color) {
  if (data.number().value() > 0) {
    float_stats_.push_back(BodyStatisticsOptimizerStates());
    atom::Hierarchy cur_root
      = atom::Hierarchy::setup_particle(new Particle(get_m()));
    IMP_LOG(TERSE, "   type " << type.get_string() <<std::endl);
    cur_root->set_name(type.get_string());
    atom::Hierarchy(get_root()).add_child(cur_root);
    ParticlesTemp cur;
    for (int j=0; j< data.number().value(); ++j) {
      double dc=data.d_factor().value();
      cur.push_back(create_particle(this, data.radius().value(),
                                    angular_d_factor_, dc,
                                    color,
                                    type, type.get_string()));
      IMP_NEW(BodyStatisticsOptimizerState, os, (cur.back()));
      os->set_period(statistics_interval_);
      float_stats_.back().push_back(os);
      cur_root.add_child(atom::Hierarchy::setup_particle(cur.back()));
      if (data.interactions().value() >0) {
        int nsites= data.interactions().value();
        set_sites(type, nsites, data.radius().value());
      }
    }
    particles_[type]=cur;
    interaction_range_factors_[type]= data.interaction_range_factor().value();
    interaction_k_factors_[type]= data.interaction_k_factor().value();
  }
}

/**
   Adds the FG Nup chains to the model hierarchy,
   based on the settings in data
*/
void SimulationData::
create_fgs(const ::npctransport::Assignment_FGAssignment&data,
           core::ParticleType type) {
  if (data.number().value() > 0) {
    fgs_stats_.push_back(base::Vector<BodyStatisticsOptimizerStates>());
    chain_stats_.push_back(ChainStatisticsOptimizerStates());
    ParticlesTemp cur;
    atom::Hierarchy hi= atom::Hierarchy::setup_particle(new Particle(get_m()));
    hi->set_name(type.get_string());
    atom::Hierarchy(get_root()).add_child(hi);
    for (int j=0; j< data.number().value(); ++j) {
      double dc=data.d_factor().value();
      double rlf=data.rest_length_factor().value();
      atom::Hierarchy hc(create_chain(this,
                                      data.number_of_beads().value(),
                                      data.radius().value(),
                                      angular_d_factor_, dc, rlf,
                                      backbone_k_,
                                      display::Color(.3,.3,.3), type, "fg",
                                      bonds_));
      cur.push_back(hc);
      ParticlesTemp chain=hc.get_children();
      chain_stats_.back()
        .push_back(new ChainStatisticsOptimizerState(chain));
      chain_stats_.back().back()->set_period(statistics_interval_);
      fgs_stats_.back().push_back(BodyStatisticsOptimizerStates());
      for (unsigned int k=0; k < chain.size(); ++k) {
        fgs_stats_.back().back()
          .push_back(new BodyStatisticsOptimizerState(chain[k]));
        fgs_stats_.back().back().back()->set_period(statistics_interval_);
      }
      hi.add_child(atom::Hierarchy(cur.back()));
      if (data.interactions().value() > 0) {
        int nsites= data.interactions().value();
        set_sites(type, nsites, data.radius().value());
      }
    }
    particles_[type]=cur;
    interaction_range_factors_[type]= data.interaction_range_factor().value();
    interaction_k_factors_[type]= data.interaction_k_factor().value();
  }
}


/**
   Creates bounding volume restraints such as box restraint and slab restraints,
   based on the box_size_, slab_height_, slab_radius_, etc. class variables
*/
void SimulationData::create_bounding_box_restraint()
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

void SimulationData::create_slab_restraint() {
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


// initialize a writer that outputs the particles hierarchy
// using the name return by ::get_rmf_file_name()
//
// \exception RMF::IOException if couldn't open RMF file
rmf::SaveOptimizerState *SimulationData::get_rmf_writer() {
  if (!rmf_writer_) {
    IMP_LOG(TERSE, "Setting up dump" << std::endl);
    RMF::FileHandle fh=RMF::create_rmf_file
      (get_rmf_file_name());
    IMP_NEW(rmf::SaveOptimizerState, los,
            (fh));
    rmf_writer_=los;
    los->set_period(dump_interval_);
    add_hierarchy(fh, atom::Hierarchy(get_root()));
    IMP::rmf::add_restraints(fh, RestraintsTemp(1, get_predr()));
    IMP::rmf::add_restraints(fh, chain_restraints_);
    if(get_has_bounding_box()){
      IMP::rmf::add_restraints(fh, RestraintsTemp(1, box_restraint_));
      IMP_NEW(display::BoundingBoxGeometry, bbg, (get_box()));
      IMP::rmf::add_geometries(fh, display::Geometries(1, bbg));
    }
    if(get_has_slab()){
      IMP::rmf::add_restraints(fh, RestraintsTemp(1, slab_restraint_));
      IMP_NEW(SlabWireGeometry, slab_geometry,
              ( slab_thickness_ , tunnel_radius_ , box_side_ ) );
        IMP::rmf::add_static_geometries
          (fh, slab_geometry->get_components() );
        //      IMP_NEW(display::CylinderGeometry, cyl_geom, (get_cylinder()) );
        //      IMP::rmf::add_static_geometries
        //          (fh, display::Geometries(1, cyl_geom));
    }
  }
  return rmf_writer_;
}

void SimulationData::dump_geometry() {
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
  if(get_has_slab()){
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
  if (!rmf_writer_) return;
  get_rmf_writer();
  rmf_writer_->reset();
}

atom::BrownianDynamics *SimulationData::get_bd() {
  if (!bd_) {
    bd_=new atom::BrownianDynamics(m_);
    bd_->set_maximum_time_step(time_step_);
    bd_->set_maximum_move(range_/4);
    if (dump_interval_ > 0) {
      bd_->add_optimizer_state(get_rmf_writer());
    }
    RestraintsTemp rs= chain_restraints_;
    if(get_has_bounding_box()) rs.push_back(box_restraint_);
    if(get_has_slab()) rs.push_back(slab_restraint_);
    rs.push_back(predr_);
    IMP_NEW(core::RestraintsScoringFunction, rsf, (rs));
    bd_->set_scoring_function(rsf);
    for (unsigned int i=0; i< fgs_stats_.size(); ++i) {
      for (unsigned int j=0; j< fgs_stats_[i].size(); ++j) {
        bd_->add_optimizer_states(fgs_stats_[i][j]);
      }
    }
    for (unsigned int i=0; i< float_stats_.size(); ++i) {
      bd_->add_optimizer_states(float_stats_[i]);
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
    // returned by get_cpc()
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
(const ::npctransport::Assignment_InteractionAssignment& idata )
{
  // extract interaction params
  core::ParticleType type0(idata.type0());
  core::ParticleType type1(idata.type1());
  double base_k=interaction_k_;
  if (idata.has_interaction_k()) {
    base_k= idata.interaction_k().value();
  }
  double interaction_k= base_k
    * interaction_k_factors_.find(type0)->second
    * interaction_k_factors_.find(type1)->second;
  double base_range=interaction_range_;
  if (idata.has_interaction_range()) {
    base_range= idata.interaction_range().value();
  }
  double interaction_range= base_range
    * interaction_range_factors_.find(type0)->second
    * interaction_range_factors_.find(type1)->second;

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
  IMP_LOG( PROGRESS,
           "Interaction "
           << type0.get_string() << ", "
           << type1.get_string()
           << "  sizes: " << set0.size() << ", " << set1.size()
           << " statistics range: " << range_ << std::endl );
  if(set0.size() > 0 && set1.size() > 0) {
    InteractionType interaction_type = std::make_pair(type0,type1);
    IMP_NEW( BipartitePairsStatisticsOptimizerState,
             bpsos ,
             ( get_m(), interaction_type,
               set0, set1,  range_ /* distance thresh */ ) );
    bpsos->set_period(statistics_interval_);
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
  if(get_has_slab()){
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
  if(get_has_slab()){
    // IMP_NEW(SlabWireGeometry, slab_geometry,
    //         (1000 /* h */ , 100 /* r */, 300 /* w */) );
    // w->add_geometry(slab_geometry);
    IMP_NEW(display::CylinderGeometry, cyl_geom, (get_cylinder()) );
    w->add_geometry( cyl_geom);
  }

}

#define UPDATE(frame, message, field, newvalue)  \
  (message).set_##field(static_cast<double>(frame)/(frame+1)*(message).field() \
                    + 1.0/(frame+1)*newvalue)


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

void SimulationData::update_statistics(const boost::timer &timer) const {
  ::npctransport::Statistics stats;
  if (first_stats_) { // first initialization
    for (unsigned int i=0; i<type_of_fg.size(); ++i) {
      if (particles_.find(type_of_fg[i]) != particles_.end()) {
        stats.add_fgs();
        IMP_USAGE_CHECK(stats.fgs_size() ==static_cast<int>(i+1),
                        "Wrong size: " << stats.fgs_size());
      }
    }
    for (unsigned int i=0; i<type_of_float.size(); ++i) {
      if (particles_.find(type_of_float[i]) != particles_.end()) {
        stats.add_floaters();
      }
    }
    for(unsigned int i = 0; i < interactions_stats_.size(); i++){
      ::npctransport::Statistics_InteractionStats*
        pOutStats_i = stats.add_interactions();
      InteractionType itype = interactions_stats_[i]->get_interaction_type();
      pOutStats_i->set_type0( itype.first.get_string() );
      pOutStats_i->set_type1( itype.second.get_string() );
    }
    first_stats_=false;
  } else { // not first initialization
    std::ifstream inf(statistics_file_name_.c_str(), std::ios::binary);
    stats.ParseFromIstream(&inf);
    inf.close();
  }

  int nf= stats.number_of_frames();
  ParticlesTemp all;
  ParticlesTemps floaters;
  base::Vector<ParticlesTemps> fgs;
  for (unsigned int i=0; i<type_of_fg.size(); ++i) {
    if (particles_.find(type_of_fg[i]) != particles_.end()) {
      fgs.push_back(ParticlesTemps());
      ParticlesTemp cur= particles_.find(type_of_fg[i])->second;
      for (unsigned int j=0; j< cur.size(); ++j) {
        atom::Hierarchy h(cur[j]);
        ParticlesTemp chain
            = get_as<ParticlesTemp>(atom::get_leaves(h));
        all+=chain;
        fgs.back().push_back(chain);
#ifdef IMP_NPCTRANSPORT_USE_IMP_CGAL
        double volume= atom::get_volume(h);
        double radius_of_gyration= atom::get_radius_of_gyration(chain);
#else
        double volume= -1.;
        double radius_of_gyration= -1.;
#endif
        UPDATE(nf, *stats.mutable_fgs(i), volume, volume);
        double length= core::get_distance(core::XYZ(chain[0]),
                                          core::XYZ(chain.back()));
        UPDATE(nf, *stats.mutable_fgs(i), length, length);
        UPDATE(nf, *stats.mutable_fgs(i), radius_of_gyration,
               radius_of_gyration);
      }
    }
  }
  for (unsigned int i=0; i<float_stats_.size(); ++i) {
    for (unsigned int j=0; j<float_stats_[i].size(); ++j) {
      int cnf=(nf)*float_stats_[i].size()+j;
      UPDATE(cnf,
             *stats.mutable_floaters(i), diffusion_coefficient,
             float_stats_[i][j]->get_diffusion_coefficient());
      UPDATE(cnf,
             *stats.mutable_floaters(i), correlation_time,
             float_stats_[i][j]->get_correlation_time());
    }
  }
  for (unsigned int i=0; i<fgs_stats_.size(); ++i) {
    for (unsigned int j=0; j<fgs_stats_[i].size(); ++j) {
      for (unsigned int k=0; k< fgs_stats_[i][j].size(); ++k) {
        unsigned int n=fgs_stats_[i].size()*fgs_stats_[i][j].size();
        unsigned int cnf=(nf)*n+j*fgs_stats_[i][j].size()+k;;
        UPDATE(cnf,
               *stats.mutable_fgs(i),particle_correlation_time,
               fgs_stats_[i][j][k]->get_correlation_time());
        UPDATE(cnf,
               *stats.mutable_fgs(i),particle_diffusion_coefficient,
               fgs_stats_[i][j][k]->get_diffusion_coefficient());
      }
    }
  }
  for (unsigned int i=0; i<chain_stats_.size(); ++i) {
    for (unsigned int j=0; j<chain_stats_[i].size(); ++j) {
      unsigned int n=chain_stats_[i].size();
      unsigned int cnf=(nf)*n+j;
      UPDATE(cnf,
             *stats.mutable_fgs(i), chain_correlation_time,
             chain_stats_[i][j]->get_correlation_time());
      UPDATE(cnf,
             *stats.mutable_fgs(i), chain_diffusion_coefficient,
             chain_stats_[i][j]->get_diffusion_coefficient());
      Floats dfs=chain_stats_[i][j]->get_diffusion_coefficient();
      double df= std::accumulate(dfs.begin(), dfs.end(), 0.0)/dfs.size();
      UPDATE(cnf,
             *stats.mutable_fgs(i), local_diffusion_coefficient,
             df);
    }
  }

  for (unsigned int i=0; i<type_of_float.size(); ++i) {
    if (particles_.find(type_of_float[i]) != particles_.end()) {
      floaters.push_back(particles_.find(type_of_float[i])->second);
      all+=floaters.back();
      double interactions, interacting, bead_partners, chain_partners;
      boost::tie(interactions, interacting, bead_partners, chain_partners)
          =get_interactions_and_interacting(floaters.back(), fgs);
      UPDATE(nf, *stats.mutable_floaters(i), interactions,
             interactions/floaters.back().size());
      UPDATE(nf, *stats.mutable_floaters(i), interacting,
             interacting/floaters.back().size());
      UPDATE(nf, *stats.mutable_floaters(i), interaction_partner_chains,
             chain_partners/interacting);
      UPDATE(nf, *stats.mutable_floaters(i), interaction_partner_beads,
             bead_partners/interacting);
    }
  }

  // update statistics gathered on interaction rates
  for(unsigned int i = 0; i < interactions_stats_.size(); i++){
    IMP_LOG( PROGRESS, "adding interaction statistics # " << i << std::endl );
    ::npctransport::Statistics_InteractionStats*
      pOutStats_i = stats.mutable_interactions(i);
    IMP::Pointer<BipartitePairsStatisticsOptimizerState>
      pInStats_i = interactions_stats_[i];
    // verify correct interaction type is stored
    InteractionType itype = pInStats_i->get_interaction_type();
    std::string s_type0 = itype.first.get_string();
    std::string s_type1 = itype.second.get_string();
    if( pOutStats_i->type0() != s_type0  || pOutStats_i->type1() != s_type1) {
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

    UPDATE(nf, *pOutStats_i, avg_contacts_per_particle0,
           avg_contacts_num / n0);
    UPDATE(nf, *pOutStats_i, avg_contacts_per_particle1,
           avg_contacts_num / n1);
    UPDATE(nf, *pOutStats_i, avg_pct_bound_particles0,
           pInStats_i->get_average_percentage_bound_particles_1() );
    UPDATE(nf, *pOutStats_i, avg_pct_bound_particles1,
           pInStats_i->get_average_percentage_bound_particles_2() );
  }

  UPDATE(nf, stats, energy_per_particle, get_m()->evaluate(false)/all.size());

  UPDATE(nf, stats, seconds_per_iteration, timer.elapsed());

  stats.set_number_of_frames(nf+1);
  std::ofstream outf(statistics_file_name_.c_str(), std::ios::binary);
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
