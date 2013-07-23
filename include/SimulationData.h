/**
 *  \file npctransport/SimulationData.h
 *  \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_SIMULATION_DATA_H
#define IMPNPCTRANSPORT_SIMULATION_DATA_H

#include "npctransport_config.h"
#include "io.h"
#include "BodyStatisticsOptimizerState.h"
#include "ParticleTransportStatisticsOptimizerState.h"
#include "ChainStatisticsOptimizerState.h"
#include "BipartitePairsStatisticsOptimizerState.h"
#include <IMP/Model.h>
#include <IMP/atom/BrownianDynamics.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/container/ClosePairContainer.h>
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/core/pair_predicates.h>
#include "linear_distance_pair_scores.h"
#include <IMP/core/Typed.h>
#include <IMP/base/map.h>
#include <IMP/display/declare_Geometry.h>
#include <IMP/container/PairContainerSet.h>
#include <IMP/rmf/SaveOptimizerState.h>
#include "SlabSingletonScore.h"
#include <IMP/core/BoundingBox3DSingletonScore.h>
#include <IMP/container/PredicatePairsRestraint.h>
#include <IMP/base/Pointer.h>

#include <boost/timer.hpp>
#include <string>

#ifdef SWIG
namespace boost {
struct timer {};
}
#endif
#ifndef SWIG
namespace npctransport_proto {
class Assignment_FGAssignment;
class Assignment_InteractionAssignment;
class Assignment_FloaterAssignment;
class Assignment_ObstacleAssignment;
}
#endif


IMPNPCTRANSPORT_BEGIN_NAMESPACE

class IMPNPCTRANSPORTEXPORT SimulationData : public base::Object {
 private:
#ifndef SWIG
  template <class T>
  struct Parameter {
    T t_;
    bool init_;

   public:
    Parameter() : init_(false) {}
    operator T() const {
      IMP_USAGE_CHECK(init_, "Not initialized");
      return t_;
    }
    void operator=(T t) {
      t_ = t;
      init_ = true;
    }
    bool is_init() { return init_; }
  };
  Parameter<double> interaction_k_;
  Parameter<double> interaction_range_;
  Parameter<double> backbone_k_;
  Parameter<double> box_side_;
  Parameter<double> tunnel_radius_;
  Parameter<double> slab_thickness_;
  Parameter<bool> box_is_on_;
  Parameter<bool> slab_is_on_;
  Parameter<double> slack_;
  Parameter<int> number_of_trials_;
  Parameter<int> number_of_frames_;
  Parameter<int> dump_interval_frames_;
  Parameter<double> nonspecific_k_;
  Parameter<double> nonspecific_range_;
  Parameter<double> angular_d_factor_;
  Parameter<int> statistics_interval_frames_;
  Parameter<double> excluded_volume_k_;
  Parameter<double> range_;
  Parameter<double> time_step_;
  Parameter<double> statistics_fraction_;
  Parameter<double> maximum_number_of_minutes_;

  // Per type scaling factors for the interaction parameters
  base::map<core::ParticleType, double> interaction_range_factors_;
  base::map<core::ParticleType, double> interaction_k_factors_;
#endif
  // the model on which the simulation is run
  base::Pointer<Model> m_;

  // The BrownianDynamic simulator
  base::Pointer<atom::BrownianDynamics> bd_;

  base::Pointer<container::ListSingletonContainer> diffusers_;

  bool diffusers_updated_; // keeps track of whether the diffusers list is
                           // up to date (can be invalid if particles added)
  bool obstacles_updated_; // keeps track of whether the obstacles list is
                           // up to date (can be invalid if particles added)
                           // TODO: possibly need be same as diffusers - matter
                           //       of efficiency

  base::Pointer<container::ClosePairContainer> cpc_;

  // generates hash values ('predicates') for ordered types pairs
  // (e.g., pairs of ParticleTypes)
  base::Pointer<core::OrderedTypePairPredicate> otpp_;

  // contains all restraints between pairs of particle types
  base::Pointer<container::PredicatePairsRestraint> predr_;

  // a writer to an RMF (Rich Molecular Format) type file
  base::Pointer<rmf::SaveOptimizerState> rmf_sos_writer_;

  // the root of the model hierarchy
  base::Pointer<Particle> root_;

  base::Pointer<display::Geometry> static_geom_;
  base::Pointer<Restraint> box_restraint_;
  base::Pointer<Restraint> slab_restraint_;
  Restraints chain_restraints_;
  base::map<core::ParticleType, algebra::Vector3Ds> sites_;
  base::map<core::ParticleType, double> ranges_;
  LinearWellPairScores backbone_scores_;

  // the list of particles for each particle type
  // e.g., particles_[ ParticleType("fg0") ]
  base::map<core::ParticleType, ParticlesTemp> particles_;

  // the file to which simulation statistics are dumped:
  std::string output_file_name_;
  // the RMF format file to which simulation output is dumped:
  std::string rmf_file_name_;

  // statistics about all FG nups individual particles, for each FG type
  base::Vector<base::Vector<BodyStatisticsOptimizerStates> > fgs_stats_;

  // statistics about all Kaps and non-specific binders ("floats"),
  // for each floats type
  base::Vector<BodyStatisticsOptimizerStates> float_stats_;

  // transport statistics about all Kaps and non-specific binders ("floats"),
  // for each floats type
  base::Vector<ParticleTransportStatisticsOptimizerStates>
      float_transport_stats_;

  // statistics about entire FG chains, for each FG type
  base::Vector<ChainStatisticsOptimizerStates> chain_stats_;

  // statistics of pairs of interactions for each interaction type
  BipartitePairsStatisticsOptimizerStates interactions_stats_;

  // true if statistics were recently reset, so that update_statistics
  // should restart averaging from 0 frames
  mutable bool is_stats_reset_;

 private:

  boost::tuple<double, double, double, double> get_interactions_and_interacting(
      const ParticlesTemp &kaps, const base::Vector<ParticlesTemps> &fgs) const;

  /**
     Adds the FG Nup chains to the model hierarchy,
     based on the settings in data

     @param data data for fgs of this type in protobuf format as specified
                 in data/npctransport.proto
     @param default_type for backward compatibility only, a default type of particle
                         if one is not specified in data.type
   */
  void create_fgs(const ::npctransport_proto::Assignment_FGAssignment &data,
                  core::ParticleType default_type);

  /**
     Adds the 'floaters' (free diffusing particles) to the model hierarchy,
     based on the settings in data

     @param data data for floater in protobuf format as specified
                 in data/npctransport.proto
     @param default_type for backward compatibility only, a default type
                         of particle if one is not specified in data.type
     @param color a color for all floaters of this type
   */
  void create_floaters(
      const ::npctransport_proto::Assignment_FloaterAssignment &data,
      core::ParticleType default_type, display::Color color);

  /**
     Adds the 'obstacles' (possibly static e.g. nups that make the pore)
     to the model hierarchy, based on the settings in data

     @param data data for obstacles in protobuf format as specified in
                 data/npctransport.proto
  */
  void create_obstacles
    ( const ::npctransport_proto::Assignment_ObstacleAssignment &data);

  /**
     Creates bounding box restraint based on the box_size_
     class variable, and apply it to all diffusers returned
     by get_diffusers()
   */
  void create_bounding_box_restraint_on_diffusers();

  /**
     Creates slab bounding volume restraint, based on the slab_thickness_
     and tunnel_radius_ class variables, and apply it to all diffusers
     returned by get_diffusers()
   */
  void create_slab_restraint_on_diffusers();

  /** Initializes the simulation based on the specified assignment file

      @param assignments_file the file with the current assignment, which
                              is also to be used for progress and statistics
                              output
      @param quick whether to perform a quick simulation with a small number
                   of iterations
  */
  void initialize(std::string assignments_file, bool quick);

 public:
  /**
     @param[out] output_file name of protobuf file that encapsulates the
     assignment
     data for initializing the simulation params and the output statistics file.
     The output file is assumed to already exist at this point and contain the
     assignment.
     If it contains an rmf_conformation or conformation field, it will be used
     to initialize particle positions,
     in that priority.
     @param[out] quick if true, perform a very short simulation,
                       typically for calibration or for initial testing
     @param[out] rmf_file_name RMF file to which simulation trajectory
                               is recorded (if empty, no RMF file is created
     upon construction)
                               \see set_rmf_file_name()
   */
  SimulationData(std::string output_file, bool quick,
                 std::string rmf_file_name = std::string());

  Model *get_m();
#ifndef SWIG
  Model *get_m() const { return m_; }
#endif
  double get_range() const { return range_; }
  atom::BrownianDynamics *get_bd();

  /**
     returns all diffusing particles that currently exist in
     this simulation data object (or an empty list if none exists)
   */
  container::ListSingletonContainer *get_diffusers();

  /**
     returns a reference to the collection of score functions for FG backbones
     (can be used to e.g. scale them up or down during optimization)
   */
  LinearWellPairScores get_backbone_scores() const { return backbone_scores_; }

  /**
     a close pair container for all diffusers except diffusers that
     appear consecutively within the model (e.g., fg repeats)
  */
  container::ClosePairContainer *get_cpc();

  /**
     Returns the container for restraints over pairs of particles. Different
     scores
     are used for particles of different (ordered) particle types.
     When called for the first time, returns a new PredicatePairsRestraints
     over all diffusing particles and sets a default linear repulsion restraint
     between all close pairs returned by get_cpc()
  */
  container::PredicatePairsRestraint *get_predr();

  /** get the number of interactions between two particles */
  int get_number_of_interactions(Particle *a, Particle *b) const;

  /**
      Create n interaction sites spread around the surface of a ball of
      radius r, and associate these sites with all particles of type t0
  */
  void set_sites(core::ParticleType t0, unsigned int n, double r);

  /**
      returns a list of 3D coordinates for the interaction sites associated
      with particles of type t0. The site coordinates are relative to the
      particles centers.
  */
  algebra::Vector3Ds get_sites(core::ParticleType t0) const {
    if (sites_.find(t0) != sites_.end()) {
      return sites_.find(t0)->second;
    } else {
      return algebra::Vector3Ds();
    }
  }

  /** Set the sites explicitly. */
  void set_sites(core::ParticleType t0, const algebra::Vector3Ds &sites) {
    sites_[t0] = sites;
  }

  /**
      Returns the effective interaction range radius of a
      site on a floater */
  // TODO: this is not the true value - the true one might depend on the
  //       specific range of each interaction type, so range_ is more
  //       like an upper bound
  double get_site_radius(core::ParticleType) const { return range_ / 2; }

  // returns true if a bounding box restraint has been defined */
  bool get_has_bounding_box() const { return box_restraint_; }

  /** returns true if a slab restraint has been defined */
  bool get_has_slab() const { return slab_restraint_; }

  double get_statistics_fraction() const { return statistics_fraction_; }

  //! Return the maximum number of minutes the simulation can run
  /** Or 0 for no limit. */
  double get_maximum_number_of_minutes() const {
    return maximum_number_of_minutes_;
  }

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
  void add_interaction(
      const ::npctransport_proto::Assignment_InteractionAssignment &idata);
#endif

  algebra::BoundingBox3D get_box() const {
    return algebra::get_cube_d<3>(.5 * box_side_);
  }

  algebra::Cylinder3D get_cylinder() const;

  double get_backbone_k() const { return backbone_k_; }

  double get_interaction_k() const { return interaction_k_; }

  /**
   Open the specified RMF file, links it to the hierarchies of this object, and
   copy
   the XYZ coordinates from its last frame into the particles of this object.

   @note this method assumes that the hierarchies that are stored in the RMF
   file
   was constructed in the same way as the hierarchies within this SimulationData
   object. The only major difference between the two is assumed to be the
   coordinates of
   the various particles. Otherwise, results might be unexpected, and it is not
   guaranteed
   that a mismatch would be detected at runtime.
   \see IMP::rmf::link_hierarchies()
   \see IMP::rmf::load_frame()

   The frame number is the frame in the file to use, if, eg, there are several
   initial frames to use. If it is -1, the last frame is used.

   @exception RMF::IOException if couldn't open RMF file, or unsupported file
   format
  */
  void initialize_positions_from_rmf(RMF::FileConstHandle fh,
                                     int frame_number = -1);

  /** Links a handle to an open rmf file with the
      hierarchy and constraints of this simulation data,
      so that writing frames into this handle will write
      the updated state of the simulation

      @return a handle to the file

      \exception RMF::IOException IO error with RMF file handle
  */
  void link_rmf_file_handle(RMF::FileHandle fh);

  /**
     Returns the internal periodic SaveOptimizerState writer that
     periodically outputs the particles hierarchy and restraints,
     using the file name returned by SimulationData::get_rmf_file_name(), which is
     assumed to be non empty.  If the writer does not exist, then
     it is constructed first.

     \exception RMF::IOException couldn't create RMF file
  */
  rmf::SaveOptimizerState *get_rmf_sos_writer();

  /**
     Resets the RMF file to which the internal periodic SaveOptimizerState
     writer dumps the output. If it does not exist, create it.

     \exception RMF::IOException couldn't create RMF file or other IO related
                                 problems with RMF
   */
  void reset_rmf();

  void reset_statistics_optimizer_states();

  void add_chain_restraint(Restraint *r) { chain_restraints_.push_back(r); }

  /**
     Write the geometry to the file path 'out'
  */
  void write_geometry(std::string out);

  void dump_geometry();

  /**
      Returns the list of particles in the simulation model
      of type 'type', or empty list if not found
  */
  ParticlesTemp get_particles(core::ParticleType type) const {
    base::map<core::ParticleType, ParticlesTemp>::const_iterator iter;
    iter = particles_.find(type);
    if (iter != particles_.end()) return iter->second;
    return ParticlesTemp();
  }

  void set_interrupted(bool tf);

  /**
      opens / creates statistics protobuf file, and update it
      with appropriate statistics, using statistics file
      originally specified in the constructor.

      @param timer the timer that was used to measure the time
                   that has elapsed for statistics
      @param nf_new the number of frames by which the statistics file
                    should be advanced. This is used to weight the
                    contribution of average statistics over time.

      @note this method is not const cause it may invoke e.g., energy evaluation
            though it does not substantially change anything in the state of the object
   */
  void update_statistics(const boost::timer &timer,
                         unsigned int nf_new = 1);

  unsigned int get_number_of_frames() const { return number_of_frames_; }

  unsigned int get_number_of_trials() const { return number_of_trials_; }

  atom::Hierarchy get_root() const { return atom::Hierarchy(root_); }

  Restraints get_chain_restraints() const { return chain_restraints_; }

  Restraint *get_box_restraint() const { return box_restraint_; }

  Restraint *get_slab_restraint() const { return slab_restraint_; }

  double get_slab_thickness() const { return slab_thickness_; }

  display::Geometry *get_static_geometry();

  int get_rmf_dump_interval_frames() const { return dump_interval_frames_; }

  std::string get_rmf_file_name() const { return rmf_file_name_; }

  /**
     resets the name of the RMF file that records the simulation.  If
     the new name is different than the old one, then the previous RMF
     writer is invalidated by this action (closed and flushed).  If
     the Brownian Dynamics object was already initialized, a new
     writer with the new name is added as an optimizer state to it
     instead of the existing one.
     TODO: make sure the old writer is indeed closed and flushed
  */
  void set_rmf_file_name(const std::string &new_name);

  IMP_OBJECT_METHODS(SimulationData);
};

inline IMPNPCTRANSPORTEXPORT boost::timer create_boost_timer() {
  return boost::timer();
}

inline WeakObjectKey get_simulation_data_key() {
  static WeakObjectKey simdata("simulation data");
  return simdata;
}

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SIMULATION_DATA_H */
