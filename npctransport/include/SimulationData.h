/**
 *  \file simulation_data.h
 *  \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_SIMULATION_DATA_H
#define IMPNPCTRANSPORT_SIMULATION_DATA_H

#include "npctransport_config.h"
#include "io.h"
#include "statistics_optimizer_states.h"
#include <IMP/Model.h>
#include <IMP/atom/BrownianDynamics.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/container/ClosePairContainer.h>
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/core/pair_predicates.h>
#include <IMP/core/Typed.h>
#include <IMP/display/declare_Geometry.h>
#include <IMP/container/PairContainerSet.h>
#include <IMP/rmf/SaveOptimizerState.h>
#include "SlabSingletonScore.h"
#include <IMP/core/BoundingBox3DSingletonScore.h>
#include <IMP/container/PredicatePairsRestraint.h>
#include <IMP/Pointer.h>

#include <boost/timer.hpp>
#include <string>

//#ifndef SWIG
#ifdef IMP_NPC_GOOGLE
#include "third_party/npc/npctransport/data/npctransport.pb.h"
#else
#include "npctransport.pb.h"
#endif
//#endif
IMPNPCTRANSPORT_BEGIN_NAMESPACE


class IMPNPCTRANSPORTEXPORT SimulationData: public base::Object {
#ifndef SWIG
  template <class T>
    struct Parameter
    {
      T t_;
      bool init_;
    public:
    Parameter(): init_(false){}
      operator T() const {
        IMP_USAGE_CHECK(init_, "Not initialized");
        return t_;
      }
      void operator=(T t) {
        t_=t;
        init_=true;
      }
        bool is_init() { return init_; }
    };
  Parameter<double> interaction_k_;
  Parameter<double> interaction_range_;
  Parameter<double> backbone_k_;
  Parameter<double> box_side_;
  Parameter<double> tunnel_radius_;
  Parameter<double> slab_thickness_;
  Parameter<double> slack_;
  Parameter<int> number_of_trials_;
  Parameter<int> number_of_frames_;
  Parameter<int> dump_interval_;
  Parameter<double> nonspecific_k_;
  Parameter<double> nonspecific_range_;
  Parameter<double> angular_d_factor_;
  Parameter<int> statistics_interval_;
  Parameter<double> excluded_volume_k_;
  Parameter<double> range_;
  Parameter<double> time_step_;

  // Per type scaling factors for the interaction parameters
  compatibility::map<core::ParticleType, double> interaction_range_factors_;
  compatibility::map<core::ParticleType, double> interaction_k_factors_;
#endif
  // the model on which the simulation is run
  Pointer<Model> m_;

  // The BrownianDynamic simulator
  Pointer<atom::BrownianDynamics> bd_;

  Pointer<container::ListSingletonContainer> diffusers_;

  Pointer<container::ClosePairContainer> cpc_;

  // generates hash values ('predicates') for ordered types pairs
  // (e.g., pairs of ParticleTypes)
  Pointer<core::OrderedTypePairPredicate> otpp_;

  // contains all restraints between pairs of particle types
  Pointer<container::PredicatePairsRestraint> predr_;

  // a writer to an RMF (Rich Molecular Format) type file
  Pointer<rmf::SaveOptimizerState> rmf_writer_;

  // the root of the model hierarchy
  Pointer<Particle> root_;

  Pointer<display::Geometry> static_geom_;
  Pointer<container::PairContainerSet> bonds_;
  Pointer<Restraint> box_restraint_;
  Pointer<Restraint> slab_restraint_;
  Restraints chain_restraints_;
  compatibility::map<core::ParticleType, algebra::Vector3Ds> sites_;
  compatibility::map<core::ParticleType, double> ranges_;

  // the list of particles for each particle type
  // e.g., particles_[ ParticleType("fg0") ]
  compatibility::map<core::ParticleType, ParticlesTemp> particles_;

  mutable bool first_stats_;
  // the file to which simulation statistics are dumped:
  std::string statistics_file_name_;
  // the RMF format file to which simulation output is dumped:
  std::string rmf_file_name_;

  // statistics about all FG nups individual particles, for each FG type
  base::Vector<base::Vector<BodyStatisticsOptimizerStates> > fgs_stats_;

  // statistics about all Kaps and non-specific binders ("floats"),
  // for each floats type
  base::Vector<BodyStatisticsOptimizerStates> float_stats_;

  // statistics about entire FG chains, for each FG type
  base::Vector<ChainStatisticsOptimizerStates> chain_stats_;

  // statistics of pairs of interactions for each interaction type
  BipartitePairsStatisticsOptimizerStates interactions_stats_;

  boost::tuple<double,double,double,double>
      get_interactions_and_interacting(const ParticlesTemp &kaps,
                                       const base::Vector<ParticlesTemps> &fgs)
    const;

  /**
     Adds the FG Nup chains to the model hierarchy,
     based on the settings in data
   */
  void create_fgs(const  ::npctransport::Assignment_FGAssignment&data,
                  core::ParticleType type);

  /**
     Adds the 'floaters' (free diffusing particles) to the model hierarchy,
     based on the settings in data
   */
  void create_floaters(const ::npctransport::Assignment_FloaterAssignment&data,
                       core::ParticleType type, display::Color color);

  /**
     Creates bounding volume restraints such as box restraint and slab
     restraints, based on the box_size_, slab_height_ and slab_radius_
     class variables, etc.
   */
  void create_bounding_box_restraint();

  void create_slab_restraint();

 public:
  /**
     @param[out] assignment_file name of input assignment file for
                                 initializing the simulation params
     @param[out] statistics_file name of output statistics file
     @param[out] quick if true, perform a very short simulation,
                       typically for calibration or for initial testing
     @param[out] rmf_file_name RMF file to which simulation trajectory
                               is recorded
   */
  SimulationData(std::string assignment_file,
                 std::string statistics_file,
                 bool quick,
                 std::string rmf_file_name = "output.prmf");

  Model *get_m();
#ifndef SWIG
  Model *get_m() const  {return m_;}
#endif
  double get_range() const {return range_;}
  atom::BrownianDynamics *get_bd();
  container::ListSingletonContainer *get_diffusers();

  // a close pair container for all diffusers except diffusers that
  // appear consecutively within the model (e.g., fg repeats)
  container::ClosePairContainer* get_cpc();

  container::PredicatePairsRestraint* get_predr();
  // get the number of interactions between two particles
  int get_number_of_interactions(Particle *a, Particle *b) const;
  void set_sites(core::ParticleType t0,
                 unsigned int n, double r);
  algebra::Vector3Ds get_sites(core::ParticleType t0) const {
    if (sites_.find(t0) != sites_.end()) {
      return sites_.find(t0)->second;
    } else {
      return algebra::Vector3Ds();
    }
  }
  double get_site_radius(core::ParticleType) const {
    return range_/2;
  }
  PairContainer* get_bonds() const {return bonds_;}
  bool get_has_bounding_box() const { return box_restraint_; }
  bool get_has_slab() const { return slab_restraint_; }

  // swig doesn't equate the two protobuf types
#ifndef SWIG
  /**
     add a SitesPairScore restraint that applies to particles of
     types t0 and t1 (specified in idata) to the PredicatePairsRestraint
     object that is returned by get_predr().
     The restraint is symmetric, and applies to both (t0,t1) and (t1,t0).
     In additions, adds statistics optimizer state about pairs of this type.

     A SitesPairScore restraint means site-specific
     attractive forces between surface bidning sites on each particle +
     non-specific attraction and soft repulsion between the entire particles.

     @param idata the protobuf data about the interaction (particle types,
            interaction coefficients, etc.)
  */
  void add_interaction
    ( const ::npctransport::Assignment_InteractionAssignment& idata );
#endif

  algebra::BoundingBox3D get_box() const {
    return algebra::get_cube_d<3>(.5*box_side_);
  }
  algebra::Cylinder3D get_cylinder() const;
  double get_backbone_k() const {return backbone_k_;}
  double get_interaction_k() const {return interaction_k_;}

  /**
   Initialize a writer that outputs the particles hierarchy
   using the name returned by ::get_rmf_file_name()

   \exception RMF::IOException couldn't create RMF file
  */
  rmf::SaveOptimizerState *get_rmf_writer();

  void reset_rmf();
  void add_chain_restraint(Restraint *r) {
    chain_restraints_.push_back(r);
  }
  void write_geometry(std::string out);
  void dump_geometry();
  ParticlesTemp get_particles(core::ParticleType type) const {
    return particles_.find(type)->second;
  }
  void update_statistics(const boost::timer &timer) const;
  unsigned int get_number_of_frames() const {return number_of_frames_;}
  unsigned int get_number_of_trials() const {return number_of_trials_;}
  atom::Hierarchy get_root() const {return atom::Hierarchy(root_);}
  Restraints get_chain_restraints() const {return chain_restraints_;}
  Restraint* get_box_restraint() const {return box_restraint_;}
  Restraint* get_slab_restraint() const {return slab_restraint_;}
  display::Geometry* get_static_geometry();

  int get_rmf_dump_interval() const { return dump_interval_; }
  std::string get_rmf_file_name() const { return rmf_file_name_; }

  /**
     resets the name of the RMF file that records the simulation.
     the previous RMF writer is invalidated by this action (closed and flushed)
     TODO: make sure it is indeed closed and flushed
  */
  void set_rmf_file_name(const std::string& new_name)
  {
    rmf_file_name_ = new_name;
    rmf_writer_ = IMP_NULLPTR; // invalidate the existing writer
  }

  IMP_OBJECT_INLINE(SimulationData,IMP_UNUSED(out),);
};

inline IMPNPCTRANSPORTEXPORT
boost::timer create_boost_timer()
{
  return boost::timer();
}



inline WeakObjectKey get_simulation_data_key() {
  static WeakObjectKey simdata("simulation data");
  return simdata;
}

IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_SIMULATION_DATA_H */
