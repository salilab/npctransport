/**
 *  \file npctransport/SimulationData.h
 *  \brief description
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_SIMULATION_DATA_H
#define IMPNPCTRANSPORT_SIMULATION_DATA_H


#include "npctransport_config.h"
#include <IMP/Model.h>
#include <IMP/PairContainer.h>
#include <IMP/atom/BrownianDynamics.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/SingletonContainer.h>
#include <IMP/container/CloseBipartitePairContainer.h>
//#include <IMP/container/DynamicListSingletonContainer.h>
#include <IMP/container/PairContainerSet.h>
#include <IMP/container/PredicatePairsRestraint.h>
#include <IMP/core/pair_predicates.h>
#include <IMP/core/BoundingBox3DSingletonScore.h>
#include <IMP/core/Typed.h>
#include <IMP/display/declare_Geometry.h>
#include <IMP/rmf/SaveOptimizerState.h>
#include <IMP/Pointer.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/algebra/Sphere3D.h>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include "io.h"
#include "Parameter.h"
#include "Scoring.h"
#include "Statistics.h"
#include "npctransport_proto.fwd.h"
#include <string>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

// Version 2.5 - turned slab into a particle
// Version 3.0 - added harmonic bond with proper tau
// Version 4.0 - more properly handle various suffixes of FG chains (=different types of beads)
#define IMPNPCTRANSPORT_VERSION 4.0

//! Store all parameters for a simulation
class IMPNPCTRANSPORTEXPORT SimulationData : public Object {
 private:
  // params
  Parameter<double> output_npctransport_version_;
  Parameter<double> box_side_;
  Parameter<double> tunnel_radius_; // note this is the initial pore radius (major radius if toroidal pore)
  Parameter<double> tunnel_radius_k_; // k for harmonic restraint on pore radius
  Parameter<double> pore_anchored_beads_k_; // k for harmonic restraint on beads anchored to pore
  Parameter<double> slab_thickness_; // note this is the initial thickness (also minor vertical radius if toroidal pore)
  Parameter<int> box_is_on_;
  Parameter<int> slab_is_on_;
  Parameter<int> number_of_trials_;
  Parameter<int> number_of_frames_;
  Parameter<int> dump_interval_frames_;
  Parameter<bool> is_backbone_harmonic_;
  Parameter<double> backbone_tau_ns_;
  Parameter<double> angular_d_factor_;
  Parameter<double> range_;
  Parameter<double> statistics_fraction_;
  Parameter<int> statistics_interval_frames_;
  Parameter<int> output_statistics_interval_frames_;
  Parameter<int> full_output_statistics_interval_factor_;
  Parameter<int> is_multiple_hdf5s_;
  Parameter<double> time_step_;
  Parameter<double> time_step_wave_factor_;
  Parameter<double> maximum_number_of_minutes_;
  Parameter<double> fg_anchor_inflate_factor_; // By how much to inflate static
                                               // anchors of FG nups.
  Parameter<int> is_exclude_floaters_from_slab_initially_;
  Parameter<double> are_floaters_on_one_slab_side_;
  Parameter<int> is_xyz_hist_stats_;
  Parameter<double> xyz_stats_crop_factor_;
  Parameter<double> xyz_stats_voxel_size_a_;
  Parameter<double> xyz_stats_max_box_size_a_;
 // time when simulation has started for this process
  Parameter<double> initial_simulation_time_ns_;
  Parameter<double> temperature_k_;

 public:
  double get_output_npctransport_version() const { return output_npctransport_version_; }

  /** returns the maximal interaction range between particles */
  double get_range() const { return range_; }

  int get_output_statistics_interval_frames() const
  { return output_statistics_interval_frames_; }

  int get_full_output_statistics_interval_factor() const
  { return full_output_statistics_interval_factor_; }

  // if true, hdf5s are generated every get_output_statistics_interval_frames() frames
  bool get_is_multiple_hdf5s()
  {
    return is_multiple_hdf5s_;
  }

  /** returns whether should exclude floaters from slab during initialization */
  bool get_is_exclude_floaters_from_slab_initially() const
  { return is_exclude_floaters_from_slab_initially_; }

  /** returns whether should exclude floaters from one slab side,
      if excluded at all
  */
  bool get_are_floaters_on_one_slab_side() const
  { return are_floaters_on_one_slab_side_;}

  //! whether xyz histogram statistics are on (into HDF5 file)
  bool get_is_xyz_hist_stats() const
  { return is_xyz_hist_stats_; }

  //! return the factor by which the XYZ histogram is cropped on each axis
  //! (symmetrically), up to the maximal size
  double get_xyz_stats_crop_factor() const
  { return xyz_stats_crop_factor_; }

  //! return the voxel size in angstroms for the XYZ histogram
  double get_xyz_stats_voxel_size_A() const
  { return xyz_stats_voxel_size_a_; }

  //! returns the maximal box size in angstroms for the XYZ histogram
  //! i.e. crop at most this many angstroms from each dimension
  double get_xyz_stats_max_box_size_A() const
  { return xyz_stats_max_box_size_a_; }

  /** returns the simulation angular d factor */
  double get_angular_d_factor() const
  { return angular_d_factor_; }


 private:
  // the model on which the simulation is run
  PointerMember<Model> m_;

  // The BrownianDynamic simulator
  PointerMember<atom::BrownianDynamics> bd_;

  // The scoring function wrapper for the simulation
  PointerMember <IMP::npctransport::Scoring >
    scoring_;

  // The statistics manager for this simulation
  PointerMember <IMP::npctransport::Statistics >
    statistics_;

  // all beads in the simulation (=fine-level particles)
  Particles beads_;

  // a writer to an RMF (Rich Molecular Format) type file
  PointerMember<rmf::SaveOptimizerState> rmf_sos_writer_;

  // the root of the model hierarchy
  PointerMember<Particle> root_;

  // Membrane slab, if exists (mutable but is expected to be accessed only from get_slab_particle()
  mutable PointerMember<Particle> slab_particle_;

  // fg bead types  - a set of all fg/floater/obstacle types that were
  // added via create_fgs/floaters/obstacles(), including
  // the suffixes of individual beads in a chain added to their chain root name
  ParticleTypeSet fg_bead_types_;
  // fg chain types - same as fg_bead_types but just for the root particles
  // without the suffixes
  ParticleTypeSet fg_chain_types_;
  ParticleTypeSet floater_types_;
  ParticleTypeSet obstacle_types_;

  PointerMember<display::Geometry> static_geom_;

  boost::unordered_map<core::ParticleType, algebra::Sphere3Ds> sites_;

  boost::unordered_map<core::ParticleType, double> ranges_;

  // the RMF format file to which simulation output is dumped:
  std::string rmf_file_name_;

  // whether to save restraints to rmf when
  // calling get_rmf_sos_writer()
  bool is_save_restraints_to_rmf_;

 private:

  //! initialize slab_particle_ based on current parameters
  void create_slab_particle();

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
     @param color a color for all floaters of this type
   */
  void create_floaters(
      const ::npctransport_proto::Assignment_FloaterAssignment &f_data,
      display::Color color);

  /**
     Adds the 'obstacles' (possibly static e.g. nups that make the pore)
     to the model hierarchy, based on the settings in data

     @param o_data data for obstacles in protobuf format as specified in
                 data/npctransport.proto
  */
  void create_obstacles
    ( const ::npctransport_proto::Assignment_ObstacleAssignment &o_data);

  /** Initializes the simulation based on the specified assignment file

      @param prev_output_file protobuf file with the current assignment, and
                              possibly previous progress and statistics info
     @param[in] new_output_file new protobuf file for subsequent output and
                                statistics.
                                If it is "", then prev_output_file is being
                                rewritten.
      @param quick whether to perform a quick simulation with a small number
                   of iterations
  */
  void initialize(std::string prev_output_file,
                  std::string new_output_file,
                  bool quick);

 public:
  /**
     @param[in] prev_output_file name of protobuf file that encapsulates the
                                  assignment data for initializing the
                                  simulation params and the output
                                  statistics file.  The output file is
                                  assumed to already exist at this
                                  point and contain the assignment.
                                  If it contains an rmf_conformation
                                  or conformation field, it will be
                                  used to initialize particle
                                  positions, in that priority.
     @param[in] quick if true, perform a very short simulation,
                       typically for calibration or for initial testing
     @param[out] rmf_file_name RMF file to which simulation trajectory
                               is recorded (if empty, no RMF file is created
                               upon construction)
                               \see set_rmf_file()
     @param[in] new_output_file new protobuf file for subsequent output and statistics
                                If it is "", then prev_output_file is being rewritten
                                (= default).
   */
  SimulationData(std::string prev_output_file, bool quick,
                 std::string rmf_file_name = std::string(),
                 std::string new_output_file = "");

  Model *get_model();
#ifndef SWIG
  Model * get_model() const {
    IMP_USAGE_CHECK(m_, "model not initialized in get_model()");
    return m_;
  }
#endif

  /** returns a scoring object that is updated with the current particles
      hierarch and associated with this sd
  */
 Scoring * get_scoring();

  /** returns a Statistics object that is updated by this simulation data
  */
 Statistics* get_statistics();

#ifndef SWIG
  /** returns a scoring object that is updated with the current particles
      hierarchy and associated with this sd
  */
 Scoring const* get_scoring() const;

  /** returns a Statistics object that is updated by this simulation data
   */
 Statistics const* get_statistics() const;

#endif

  /** gets the Brownian Dynamics object that is capable of simulating
      this data. Statistics are not yet active for the simulation.

      @param recreate if true, forces recreation of the bd object
  */
 atom::BrownianDynamics *get_bd(bool recreate = false);

  //! activates Brownian Dynamics statistics tracking
 //! by adding all appropriate optimizer states, if they weren't already
  void activate_statistics();

  /** returns the requested fraction of time for taking statistics */
  double get_statistics_fraction() const { return statistics_fraction_; }

  /** returns true if bead has an fg type
      that is, particle was added within create_fgs()
  */
  bool get_is_fg_bead(ParticleIndex pi) const;


  /** returns true if particle type is of fg type
      (that is, it is one of the types added via create_fgs())
  */
  bool get_is_fg_bead(core::ParticleType pt) const{
    return fg_bead_types_.find(pt) != fg_bead_types_.end();
  }

  /** returns true if particle has an fg chain type
      that is, particle has the same type as the root of an
      fg chain that was added via create_fgs()
  */
  bool get_is_fg_chain(ParticleIndex pi) const;


  /** returns true if particle type is an fg chain type
      (that is, it is one of the types of a chain added
      via create_fgs())
  */
  bool get_is_fg_chain(core::ParticleType pt) const{
    return fg_chain_types_.find(pt) != fg_chain_types_.end();
  }


  /** return all the types of fg beads that were added
      via create_fgs() including the suffix appended to their
      chain type name */
  ParticleTypeSet const& get_fg_bead_types() const {
    return fg_bead_types_;
  }

  /** return all the types of fg chains that were added
      via create_fgs(), i.e. the types of chain roots */
  ParticleTypeSet const& get_fg_chain_types() const {
    return fg_chain_types_;
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
             that stand for individual FG chains
  */
  atom::Hierarchies get_fg_chain_roots() const;

  /** \deprecated_at{2.2}, use get_fg_chain_roots() instead */
  IMPCORE_DEPRECATED_METHOD_DECL(2.2)
  atom::Hierarchies get_fg_chains() const
    { return get_fg_chain_roots(); }


  /** return all the obstacle particles */
  ParticlesTemp get_obstacle_particles() const;


#ifndef SWIG

  /**
     returns all diffusing fine-level beads in this
     simulation data object (or an empty list if none exists)
     @note efficient - returns an existing container by ref
   */
  Particles& get_beads_byref(){ return beads_; }
#endif

  /**
     returns all diffusing fine-level beads in this
     simulation data object (or an empty list if none exists)
   */
  ParticlesTemp get_beads(){ return beads_; }


  /**
     returns all fine-level beads that currently exist in
     this simulation data object get_sd(), which are also optimizable
     (or an empty list if no such bead exists)
   */
  ParticlesTemp get_optimizable_beads();

  /**
     returns all fine-level beads that currently exist in
     this simulation data object get_sd(), which are also non-optimizable
     (or an empty list if no such bead exists)
   */
  ParticlesTemp get_non_optimizable_beads();

  bool get_is_backbone_harmonic() const
  { return is_backbone_harmonic_; }

  double get_backbone_tau_ns() const {
    return backbone_tau_ns_;
  }

  double get_temperature_k() const
  { return temperature_k_; }

  /** get time from which simulation begins */
  double get_initial_simulation_time_ns() const
  {return initial_simulation_time_ns_; }

  /**
      Create n interaction sites spread around the surface of a ball of
      radius r, and associate these sites with all particles of type t0

      If n==-1, creates a single site, centered at the bead, and r is ignored

      @param t0 type of particle to associate with sites
      @param n number of sites
      @param r radius at which to position sites relative to t0 origin
      @param sr radius from site center within which site
                has maximal interaction energy
  */
  void set_sites(core::ParticleType t0, int n, double r, double sr);

  /**
      returns a list of 3D coordinates for the interaction sites associated
      with particles of type t0. The site coordinates are in the
      particles local reference frame
  */
  algebra::Vector3Ds get_site_centers(core::ParticleType t0) const {
    algebra::Vector3Ds ret;
    if (sites_.find(t0) != sites_.end()) {
      algebra::Sphere3Ds sites = sites_.find(t0)->second;
      for(unsigned int i = 0; i < sites.size(); i++){
        ret.push_back(sites[i].get_center());
      }
      return ret;
    } else {
      return algebra::Vector3Ds();
    }
  }

  /**
      returns a list of Spheres3D for the interaction sites associated
      with particles of type t0. The site coordinates are in the
      particles local reference frame.
  */
  algebra::Sphere3Ds get_sites(core::ParticleType t0) const {
    if (sites_.find(t0) != sites_.end()) {
      return sites_.find(t0)->second;
    } else {
      return algebra::Sphere3Ds();
    }
  }

  /** Set the sites explicitly. */
  void set_sites(core::ParticleType t0, const algebra::Sphere3Ds &sites) {
    IMP_USAGE_CHECK(sites.size()>0, "trying to set zero sites for particle type"
                    << t0.get_string() );
    sites_[t0] = sites;
  }

  /**
      Returns the effective interaction range radius of a
      site on a floater */
  // TODO: this is not the true value - the true one might depend on the
  //       specific range of each interaction type, so range_ is more
  //       like an upper bound
  double get_site_display_radius(core::ParticleType) const { return std::max(range_ / 2, 6.0); }

  //! Return the maximum number of minutes the simulation can run
  /** Or 0 for no limit. */
  double get_maximum_number_of_minutes() const {
    return maximum_number_of_minutes_;
  }

  //! remove particle of type pt from simulation, including all
  //! model particles, statistics and interactions associated with it
  //! \note this does not update the protobuf assignment or protobuf statistics
  //! \see remove_fgs_with_prefix
  void remove_particle_type(core::ParticleType pt);

  //! remove all FGs of specified type.
  /**
      Remove any FG chain if the prefix of its string representation matches s_fg_type
      (= begins with s_fg_type, followed by an empty string or a non-digit
       character) from the simulation, including model particles, statistics and
       associated interaction, as well as related assignment and statistics data
       from the output protobuf file

       @param s_fg_type The FG type prefix (chain name, to which a suffix could
                        be added for children beads)

       \see remove_particle_type
  */
  void remove_fgs_with_prefix(std::string s_fg_type);

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
  algebra::BoundingBox3D get_bounding_box() const;

  // get the bounding box for this simulation
  algebra::Sphere3D get_bounding_sphere() const;

  //! returns the length in A of one side of the box
  //! (valid only if get_has_bounding_box() is true)
  double get_bounding_box_size() const;

  //! sets the length in A of one side of the box
  //! (valid only if get_has_bounding_box() is true)
  void set_bounding_box_size(double box_size);

  //! get the radius of the bounding sphere
  //! (valid only if get_has_bounding_sphere() is true)
  double get_bounding_sphere_radius() const;

  //! set radius of bounding sphere in angstrom
  //! (valid only if get_has_bounding_sphere() is true)
  void set_bounding_sphere_radius(double sphere_radius);

  double get_bounding_volume() const;

  //! sets the volume of the bounding box or sphere in A^3
  //! (box side or sphere radius are adjusted accordingly)
  void set_bounding_volume(double volume_A3);

  //* returns true if a slab is defnied */
  bool get_has_slab() const { return slab_is_on_!=0; }

    /** returns true if a slab is defined and has a cylindrical pore */
  bool get_is_slab_with_cylindrical_pore() const
  { return slab_is_on_==1; }

  /** returns true if a slab is defined and has a toroidal pore */
  bool get_is_slab_with_toroidal_pore() const
  { return slab_is_on_==2; }

  /** returns the slab particle, initializing it based on current parameters
      if needed
   */
  Particle* get_slab_particle() const{
    if(!slab_particle_ && get_has_slab()){
      // Const cast below - cause the initialization of slab_particle_ is transparent when
      // using get_slab_particle(), and it should only be accessed from this method
      const_cast<SimulationData*>(this)->create_slab_particle();
    }
    return slab_particle_.get();
  }

  // get the cylinder in the slab for this simulation
  algebra::Cylinder3D get_cylinder() const;

  //! returns true if simulation has a bounding simulation box
  bool get_has_bounding_box() const {
    return box_is_on_==1;
  }

  //! returns true if simulation has a bounding simulation sphere ('cell-like')
  bool get_has_bounding_sphere() const {
    return box_is_on_==2;
  }

  //! returns true if simulation has any bounding volume (box or sphere are supported)
  bool get_has_bounding_volume() const {
    bool has_bounding_volume= (box_is_on_ != 0);
    IMP_IF_CHECK(USAGE)
      if(has_bounding_volume){
        IMP_USAGE_CHECK(get_has_bounding_box() || get_has_bounding_sphere(),
                        "Invalid box_is_on value (typically defined in protobuf)");
      }
    return has_bounding_volume;

  }

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

      @param fh the handle to be linked
      @param is_restraints if true, save restraints to RMF

      \exception RMF::IOException IO error with RMF file handle
  */
  void link_rmf_file_handle(RMF::FileHandle fh, bool is_restraints = true);

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
      Returns the root of all chains of type 'type'
  */
  atom::Hierarchy get_root_of_type(core::ParticleType type) const;

  unsigned int get_number_of_frames() const { return number_of_frames_; }

  unsigned int get_number_of_trials() const { return number_of_trials_; }

  atom::Hierarchy get_root() const { return atom::Hierarchy(root_); }

  double get_slab_thickness() const;

  //! returns the current tunnel radius
  double get_tunnel_radius() const;

  //! alias to get_tunnel_radius
    double get_pore_radius() const{
    return get_tunnel_radius();
  }

  //! returns the force coefficient for keeping the tunnel radius in
  //! equilibrium (radius is dynamic only if this coefficient is positive)
  double get_tunnel_radius_k() const {
    return tunnel_radius_k_;
  }

  //! alias to get_tunnel_radius_k
  double get_pore_radius_k() const {
    return get_tunnel_radius_k();
  }

  //! returns true if pore radius can change dynamically
  bool get_is_pore_radius_dynamic() const {
    return get_has_slab() && get_tunnel_radius_k()>0.0;
  }

  double get_pore_anchored_beads_k() const {
    return pore_anchored_beads_k_;
  }

  display::Geometry *get_static_geometry();

  int get_rmf_dump_interval_frames() const { return dump_interval_frames_; }

  std::string get_rmf_file_name() const { return rmf_file_name_; }

  /**
     resets the name of the RMF file that records the simulation.  If
     an old one exists from a previous call with a different name and
     different restraints settings, then the previous RMF
     writer is invalidated by this action (closed and flushed).  If
     the Brownian Dynamics object has already been initialized, a new
     writer with the new name is added as an optimizer state to it
     instead of the existing one.
     TODO: make sure the old writer is indeed closed and flushed

     @param new_name the new name of the rmf file
     @param is_save_restraints_to_rmf whether to save restraints to this
                                      rmf file
     @param is_force_reset if true, force reset even if new_name and
                        is_save_restraints_to_rmf are equal to their
                        old values.
  */
  void set_rmf_file(const std::string &new_name,
                    bool is_save_restraints_to_rmf = true,
                    bool is_force_restart= false);

  IMP_OBJECT_METHODS(SimulationData);

};


inline WeakObjectKey get_simulation_data_key() {
  static WeakObjectKey simdata("simulation data");
  return simdata;
}

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SIMULATION_DATA_H */
