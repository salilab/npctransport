/**
 *  \file npctransport/SimulationData.h
 *  \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_SIMULATION_DATA_H
#define IMPNPCTRANSPORT_SIMULATION_DATA_H


#include "npctransport_config.h"
#include <IMP/Model.h>
#include <IMP/PairContainer.h>
#include <IMP/atom/BrownianDynamics.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/container/CloseBipartitePairContainer.h>
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/container/PairContainerSet.h>
#include <IMP/container/PredicatePairsRestraint.h>
#include <IMP/core/pair_predicates.h>
#include <IMP/core/BoundingBox3DSingletonScore.h>
#include <IMP/core/Typed.h>
#include <IMP/display/declare_Geometry.h>
#include <IMP/rmf/SaveOptimizerState.h>
#include <IMP/base/Pointer.h>
#include <IMP/base/map.h>
#include <IMP/base/set.h>
#include "io.h"
#include "Parameter.h"
#include "Scoring.h"
#include "Statistics.h"
#include "npctransport_proto.fwd.h"
#include <string>



IMPNPCTRANSPORT_BEGIN_NAMESPACE

class IMPNPCTRANSPORTEXPORT SimulationData : public base::Object {
 private:
  // params
  Parameter<double> box_side_;
  Parameter<double> tunnel_radius_;
  Parameter<double> slab_thickness_;
  Parameter<bool> box_is_on_;
  Parameter<bool> slab_is_on_;
  Parameter<int> number_of_trials_;
  Parameter<int> number_of_frames_;
  Parameter<int> dump_interval_frames_;
  Parameter<double> angular_d_factor_;
  Parameter<double> range_;
  Parameter<double> statistics_fraction_;
  Parameter<int> statistics_interval_frames_;
  Parameter<double> time_step_;
  Parameter<double> maximum_number_of_minutes_;
  Parameter<double> fg_anchor_inflate_factor_; // By how much to inflate static
                                               // anchors of FG nups.
  Parameter<double> are_floaters_on_one_slab_side_;

  // the model on which the simulation is run
  base::PointerMember<Model> m_;

  // The BrownianDynamic simulator
  base::PointerMember<atom::BrownianDynamics> bd_;

  // The scoring function wrapper for the simulation
  base::PointerMember< IMP::npctransport::Scoring >
    scoring_;

  // The statistics manger for this simulation
  base::PointerMember< IMP::npctransport::Statistics >
    statistics_;

  // keeps track of whether the diffusers list has
  // changed (can be invalid if particles added)
   bool diffusers_changed_;


  // keeps track of whether the obstacles list
  // chnaged (can be invalid if particles added)
  // TODO: possibly need be same as diffusers - matter
  //       of efficiency
   bool obstacles_changed_;

  base::PointerMember
    <IMP::container::ListSingletonContainer> optimizable_diffusers_;

   base::PointerMember<container::ListSingletonContainer> diffusers_;

  // a writer to an RMF (Rich Molecular Format) type file
  base::PointerMember<rmf::SaveOptimizerState> rmf_sos_writer_;

  // the root of the model hierarchy
  base::PointerMember<Particle> root_;

  // fg types  - a list of all fg/floater/obstacle types that were
  // added via create_fgs/floaters/obstacles(), so far
  //! a set of particle types
  ParticleTypeSet fg_types_;
  ParticleTypeSet floater_types_;
  ParticleTypeSet obstacle_types_;

  base::PointerMember<display::Geometry> static_geom_;

  base::map<core::ParticleType, algebra::Vector3Ds> sites_;

  base::map<core::ParticleType, double> ranges_;

  // the list of particles for each particle type
  // e.g., particles_[ ParticleType("fg0") ]
  base::map<core::ParticleType, ParticlesTemp> particles_;

  // the RMF format file to which simulation output is dumped:
  std::string rmf_file_name_;


 private:

  /**
     Adds the FG Nup chains to the model hierarchy,
     based on the settings in data

     @param fg_data data for fgs of this type in protobuf format as specified
                 in data/npctransport.proto
   */
  void create_fgs
    ( const ::npctransport_proto::Assignment_FGAssignment &fg_data );

  /**
     Adds the 'floaters' (free diffusing particles) to the model hierarchy,
     based on the settings in data

     @param f_data data for floater in protobuf format as specified
                 in data/npctransport.proto
     @param default_type for backward compatibility only, a default type
                         of particle if one is not specified in data.type
     @param color a color for all floaters of this type
   */
  void create_floaters(
      const ::npctransport_proto::Assignment_FloaterAssignment &f_data,
      core::ParticleType default_type, display::Color color);

  /**
     Adds the 'obstacles' (possibly static e.g. nups that make the pore)
     to the model hierarchy, based on the settings in data

     @param o_data data for obstacles in protobuf format as specified in
                 data/npctransport.proto
  */
  void create_obstacles
    ( const ::npctransport_proto::Assignment_ObstacleAssignment &o_data);

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

  Model *get_model();
#ifndef SWIG
  Model * get_model() const {
    IMP_USAGE_CHECK(m_, "model not initialized in get_model()");
    return m_;
  }
#endif

  /** returns a scoring object that is updated with the current list of
      diffusers and obstacles and associated with this sd
  */
 Scoring * get_scoring();

  /** returns a Statistics object that is updated by this simulation data
  */
 Statistics* get_statistics();


  /** gets the Brownian Dynamics object that is capable of simulating
      this data

      @param recreate if true, forces recreation of the bd object
  */
  atom::BrownianDynamics *get_bd(bool recreate = false);

  /** returns the requested fraction of time for taking statistics */
  double get_statistics_fraction() const { return statistics_fraction_; }

  /** returns the simulation angular d factor */
  double get_angular_d_factor() const
  { return angular_d_factor_; }

  bool get_are_floaters_on_one_slab_side() const
  { return are_floaters_on_one_slab_side_; }

  /** returns true if particle type is of fg type
      (that is, particle was added within create_fgs()
  */
  bool get_is_fg(ParticleIndex pi) const;


  /** returns true if particle type is of fg type
      (that is, it is one of the types added via create_fgs())
  */
  bool get_is_fg_type(core::ParticleType pt) const{
    return fg_types_.find(pt) != fg_types_.end();
  }


  /** return all the types of fgs that were added
      via create_fgs() */
  ParticleTypeSet const& get_fg_types() const {
    return fg_types_;
  }


  /** return all the types of floaters that were added
      via create_floaters() */
  ParticleTypeSet const& get_floater_types() const {
    return floater_types_;
  }

  /** return all the types of obstacles that were added
      via create_obstacles() */
  ParticleTypeSet const&  get_obstacle_types() const {
    return obstacle_types_;
  }

  /**
     retrieve fg chain hierarchies

     @return all the fg hierarchies in the simulation data object
  */
  atom::Hierarchies get_fg_chains() const
    {
      return get_fg_chains(get_root());
    }


  /**
     retrieve fg chains from root

     @param root the root under which to look for FGs

     @return all the fg hierarchies under root, which match any fg type
             that was added to this simulation previously
  */
  atom::Hierarchies get_fg_chains(atom::Hierarchy root) const;


  /** return all the obstacle particles */
  ParticlesTemp get_obstacle_particles() const;


  /**
     returns the container of all diffusing particles that currently exist in
     this simulation data object (or an empty list if none exists)
   */
  container::ListSingletonContainer *get_diffusers();

  /**
     returns all diffusing particles that currently exist in
     this simulation data object get_sd(), which are also optimizable
     (or an empty list if none exists)
   */
  container::ListSingletonContainer* get_optimizable_diffusers();

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

  double get_range() const { return range_; }

  /**
      Returns the effective interaction range radius of a
      site on a floater */
  // TODO: this is not the true value - the true one might depend on the
  //       specific range of each interaction type, so range_ is more
  //       like an upper bound
  double get_site_radius(core::ParticleType) const { return range_ / 2; }

  //! Return the maximum number of minutes the simulation can run
  /** Or 0 for no limit. */
  double get_maximum_number_of_minutes() const {
    return maximum_number_of_minutes_;
  }

// swig doesn't equate the two protobuf types
#ifndef SWIG
  /**
     add the specified interaction to the scoring.
     In additions, adds statistics optimizer state about pairs of this type.

     @param idata the protobuf data about the interaction (particle types,
            interaction coefficients, etc.)
  */
  void add_interaction(
      const ::npctransport_proto::Assignment_InteractionAssignment &idata);
#endif

  // get the bounding box for this simulation
  algebra::BoundingBox3D get_box() const;

  // get the cylinder in the slab for this simulation
  algebra::Cylinder3D get_cylinder() const;

  bool get_has_slab() const { return slab_is_on_; }

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


  /*
     temporarily suspend output to RMF file

     @param suspend if true, suspend an existing rmf save optimizer state,
            associated with the bd, otherwise restores an old one, if it
            exists or was created.
     @note the suspension will expire when calling get_bd(true)
  */
  void switch_suspend_rmf(bool suspend);

  /**
     Write the geometry to the file path 'out'
  */
  void write_geometry(std::string out);

  void dump_geometry();

  /**
      Returns the list of particles in the simulation model
      of type 'type', or empty list if not found.
      The particles returned are all sets of direct children of the
      main root (so the diffusing particles are their leaves)
  */
  ParticlesTemp get_particles_of_type(core::ParticleType type) const {
    base::map<core::ParticleType, ParticlesTemp>::const_iterator iter;
    iter = particles_.find(type);
    if (iter != particles_.end()) return iter->second;
    return ParticlesTemp();
  }

  unsigned int get_number_of_frames() const { return number_of_frames_; }

  unsigned int get_number_of_trials() const { return number_of_trials_; }

  atom::Hierarchy get_root() const { return atom::Hierarchy(root_); }

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
 private:
  /** mark that the list of diffusers may have changed recently */
  void set_diffusers_changed(bool is_changed){
    diffusers_changed_ = is_changed;
    if(diffusers_changed_){
      //      get_scoring()->update_particles();
    }
  }

  /** mark that the list of obstacles may have changed recently */
  void set_obstacles_changed(bool is_changed){
    obstacles_changed_ = is_changed;
    if(obstacles_changed_){
      //      get_scoring()->update_particles();
    }
  }

};


inline WeakObjectKey get_simulation_data_key() {
  static WeakObjectKey simdata("simulation data");
  return simdata;
}

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SIMULATION_DATA_H */
