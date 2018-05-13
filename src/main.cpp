/**
 *  \file io.cpp
 *  \brief description.
 *
 *  Copyright 2007-2018 IMP Inventors. All rights reserved.
 *
 */

// TODO: cleanup all headers - probably most are not needed

#include <IMP/npctransport/main.h>
#include <IMP/base_types.h>
#include <boost/timer.hpp>
#include <IMP/log.h>
#include <IMP/exception.h>
#include <IMP/check_macros.h>
//#include <IMP/benchmark/Profiler.h>
#include <IMP/npctransport/internal/npctransport.pb.h>
#include <IMP/npctransport/internal/boost_main.h>
#include <IMP/npctransport.h>
#include <IMP/npctransport/initialize_positions.h>
#include <IMP/npctransport/protobuf.h>
#include <IMP/npctransport/util.h>
#include <IMP/rmf/frames.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/atom/Diffusion.h>
#include <IMP/atom/estimates.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/container.h>
#include <IMP/CreateLogContext.h>
#include <IMP/random.h>
#include <IMP/nullptr.h>
#include <IMP/nullptr_macros.h>
#include <IMP/check_macros.h>
#include <IMP/exception.h>
#include <IMP/ScoringFunction.h>
#include <IMP/Model.h>

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/coded_stream.h>

#include <boost/cstdint.hpp>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <iostream>
#include <numeric>
#include <fcntl.h>
#if defined(_MSC_VER)
#include <io.h>
#else
#include <sys/types.h>
#include <sys/stat.h>
#endif


IMPNPCTRANSPORT_BEGIN_NAMESPACE
boost::int64_t work_unit = -1;
IMP::AddIntFlag work_unitadder
( "work_unit",
  "The work unit",
  &work_unit);
std::string configuration = "configuration.pb";
IMP::AddStringFlag configuration_adder
( "configuration",
  "input configuration file in protobuf format"
  " [default: %default]",
  &configuration);
std::string output = "output.pb";
IMP::AddStringFlag output_adder
( "output",
  "output assignments and statistics file in protobuf format,"
  " recording the assignment being executed"
  " [default: %default]",
  &output);
std::string restart = "";
IMP::AddStringFlag restart_adder
( "restart",
  "output file of a previous run, from which to restart"
  " this run (also initializing the final coordinates"
  " from this previous run, if they exist in the output"
  " file)"
  " [default: %default]",
  &restart);
std::string restart_fgs_only = "";
IMP::AddStringFlag restart_fgs_only_adder
( "restart_fgs_only",
  "output file of a previous run, from which to restart"
  " only the FGs (but not diffusers). The system then goes through equilibration"
  " with the new diffusers, but not the FGs"
  " [default: %default]",
  &restart_fgs_only);
std::string conformations = "";
IMP::AddStringFlag conformations_adder
( "conformations",
  "RMF file for recording the conforomations along the "
  " simulation [default: %default]",
  &conformations);
std::string init_rmffile = "";
IMP::AddStringFlag init_rmf_adder
( "init_rmffile",
  "[OBSOLETE]"
  " RMF file for initializing the simulation with its"
  " last frame (to continue a previous run). Note that"
  " this option overrides the coordinates specified in"
  " an older output file when using the --restart flag",
  &init_rmffile);
std::string final_conformations = "final_conformations.rmf";
IMP::AddStringFlag final_conformations_adder
( "final_conformations",
  "RMF file for recording the initial and final conformations "
  " [default: %default]",
  &final_conformations);
bool verbose = false;
IMP::AddBoolFlag verbose_adder
( "verbose",
  "Print more info during run",
  &verbose);
bool first_only = false;
IMP::AddBoolFlag first_only_adder
( "first_only",
  "Only do the first simulation block",
  &first_only);
bool initialize_only = false;
IMP::AddBoolFlag initialize_only_adder
( "initialize_only", "Run the initialization and then stop",
  &initialize_only);
bool show_steps = false;
IMP::AddBoolFlag show_steps_adder
( "show_steps", "Show the steps for each modified variable",
  &show_steps);
bool show_number_of_work_units = false;
IMP::AddBoolFlag show_work_units_adder
( "show_number_of_work_units",
  "Show the number of work units" ,
  &show_number_of_work_units);
double short_init_factor = 1.0;
AddFloatFlag short_init_adder
( "short_init_factor",
  "Run an abbreviated version of system initialization, which takes"
  " the specified fraction of a full initialization (or more if >1.0)"
  " [default=1.0]"
  " Note: initialization does not run in restart mode, unless --force_initialization_on_restart"
  " flag is specified explicitly",
  &short_init_factor);
bool is_force_initialization_on_restart= false;
IMP::AddBoolFlag is_force_initialization_on_restart_adder
( "force_initialization_on_restart",
  "Force re-initalization on restart, but without randomizing positions"
  " (essentially, relaxation from specified starting position)",
  &is_force_initialization_on_restart);
double short_sim_factor = 1.0;
AddFloatFlag short_sim_adder
( "short_sim_factor",
  "Run an abbreviated version of the simulation, which takes"
  " the specified fraction of a full simulation (or more if >1.0)"
  " [default=1.0]",
  &short_sim_factor);
bool no_save_restraints_to_rmf = false;
IMP::AddBoolFlag no_save_restraints_to_rmf_adder
( "no_save_restraints_to_rmf",
  "whether not to save restraints to rmf conformations and final rms,"
  " if either is applicable",
  &no_save_restraints_to_rmf);
bool is_inflate_kap28;
IMP::AddBoolFlag is_inflate_kap28_adder
( "inflate_kap28",
  "whether to inflate kap28 to radius of 150 nm + double the box size + switch it to having 16 interaction sites",
  &is_inflate_kap28);
double kap_interaction_k_factor = 1.0;
AddFloatFlag kap_interaction_k_factor_adder
( "kap_interaction_k_factor",
  "increase kap interaction factor by said ratio relative to input configuration/restart file, and update to output file",
  &kap_interaction_k_factor );

namespace {
  /*********************************** internal functions
   * *********************************/

  // TODO: move to H file?
  //! print this score for the current state of sd
  void print_score_and_positions(SimulationData *sd, bool print_positions = false,
                                 std::string header = "Score = ");

  void print_score_and_positions(SimulationData *sd, bool print_positions,
                                 std::string header) {
    IMP_OMP_PRAGMA(critical)
      std::cout << header
                << sd->get_bd()->get_scoring_function()->evaluate(false)
                << " ; PredicatePairsRestraint score = "
                << sd->get_scoring()
                    ->get_predicates_pair_restraint()->evaluate(false)
                << std::endl;
    if (print_positions) {
      Particles ps = sd->get_beads();
      for (unsigned int i = 0; i < ps.size(); i++) {
        IMP_OMP_PRAGMA(critical)
          std::cout << ps[i] << ", " << IMP::core::RigidBody(ps[i])
          .get_reference_frame() << std::endl;
      }
    }
  }

  void adjust_kap_stickiness
  ( ::npctransport_proto::Assignment* pb_assignment );

  inline void adjust_kap_stickiness
  ( ::npctransport_proto::Assignment* pb_assignment )
  {
    for (int i = 0; i < pb_assignment->floaters_size(); ++i) {
      ::npctransport_proto::Assignment_FloaterAssignment* f_data=
        pb_assignment->mutable_floaters(i);
      std::string prefix("kap");
      if( std::mismatch(prefix.begin(),prefix.end(), f_data->type().begin()).first ==
          prefix.end() ) { // = if "kap" is a prefix of f_data->type
        double new_k_factor=
          f_data->interaction_k_factor().value() * kap_interaction_k_factor;
        f_data->mutable_interaction_k_factor()->set_value( new_k_factor );
      }
    }
  }

  //! Use output file specified in ref_output_fname to load
  //! coordinates of FGs into target_sd. It is assumed that the FGs
  //! from ref_output_fname have the same topology as in target_sd
  void restart_fgs_from_reference_output_file
  ( IMP::npctransport::SimulationData* target_sd, std::string ref_output_fname )
  {
    std::cout << "Loading FGs from " << ref_output_fname << std::endl;
    IMP_NEW(SimulationData, reference_sd, (ref_output_fname, false));
    if(0){
      atom::Hierarchy h_fg0_target=
        target_sd->get_fg_chain_roots()[0].get_child(1);
      atom::Hierarchy h_fg0_reference=
        reference_sd->get_fg_chain_roots()[0].get_child(1);
      std::cout << "Coordinates of second fg bead in target and referefnce before: "
                << core::XYZ(h_fg0_target) << " "
                << core::XYZ(h_fg0_reference) << std::endl;
    }
    copy_FGs_coordinates(reference_sd.get(), target_sd);
    if(0){
      atom::Hierarchy h_fg0_target=
        target_sd->get_fg_chain_roots()[0].get_child(1);
      atom::Hierarchy h_fg0_reference=
        reference_sd->get_fg_chain_roots()[0].get_child(1);
      std::cout << "Coordinates of second fg bead in target and referefnce after: "
                << core::XYZ(h_fg0_target) << " "
                << core::XYZ(h_fg0_reference) << std::endl;
    }
  }

  /** writes the output assignment file based on the configuration parameters
      either assign a new work unit from a configuration file, or restart
      from a previous output file

      @param actual_seed the actual random seed used in the simulation,
      to be saved in the output file
  */
  inline void write_output_based_on_flags(boost::uint64_t actual_seed) {
    std::string prev_output_fname;
    if (restart.empty()) {
      int num = IMP::npctransport::assign_ranges(configuration, output, work_unit,
                                                 show_steps, actual_seed);
      if (show_number_of_work_units) {
        IMP_OMP_PRAGMA(critical)
          std::cout << "work units " << num << std::endl;
      }
      prev_output_fname= output;
    } else {  // resart.empty()
      prev_output_fname= restart;
      IMP_OMP_PRAGMA(critical)
        std::cout << "Restart simulation from " << restart << std::endl;
    }

    ::npctransport_proto::Output prev_output;
    // copy to new file to avoid modifying input file
    //std::ifstream file(restart.c_str(), std::ios::binary);
    //bool read = prev_output.ParseFromIstream(&file);
    bool read(false);
    int fd= IMP_C_OPEN(prev_output_fname.c_str(),
                       IMP_C_OPEN_FLAG(O_RDONLY) | IMP_C_OPEN_BINARY);
    if(fd!=-1){
      google::protobuf::io::FileInputStream fis(fd);
      google::protobuf::io::CodedInputStream cis(&fis);
      cis.SetTotalBytesLimit(500000000,200000000);
      read= prev_output.ParseFromCodedStream(&cis);
      IMP_C_CLOSE(fd);
    }
    IMP_ALWAYS_CHECK(read, "Couldn't read restart file " << restart << " file descriptor " << fd,
                     IMP::ValueException);
    ::npctransport_proto::Assignment* pb_assignment = prev_output.mutable_assignment();
    adjust_kap_stickiness(pb_assignment);
    if(!restart.empty()) {
      pb_assignment->set_random_seed(actual_seed);
    }
    std::ofstream outf(output.c_str(), std::ios::binary);
    prev_output.SerializeToOstream(&outf);
  }

  /**
     Run simulation <sd> for <number_of_frames> frames, in chunks of
     optimization that last <max_frames_per_chunks> frames each.
     Statistics are updated after each chunk of optimization, using
     <timer> to time the current simulation.

     @param sd simulation data used for optimization
     @param number_of_frames total number of simulation frames requires
     @param timer a timer that was reset before this simulation trial was
     initialized, to be used for tracking statistics
     @param total_time the total time that the simulation has spent.
                       The simulation would terminate at the end of an
                       optimization chunk in which this time has
                       elapsed, that is sd->get_maximum_number_of_minutes()
     @param silent_statistics if true, do not update statistics file (e.g.,
                              during equilibration)
     @param max_frames_per_chunk maximal number of frames to be simulated
                                 in a single optimization chunk

     @note if first_only==true, make a really short simulation

     @return true if succesful, false if terminated abnormally
  */
  bool run_it(SimulationData *sd, unsigned int number_of_frames,
              boost::timer &timer, boost::timer &total_time,
              bool silent_statistics = false,
              unsigned int max_frames_per_chunk = 100000) {
    // TODO: next line is a temporary hack - needed for some reason to
    // force the pair predicates to evaluate predicate pairs restraints
    sd->get_model()->update();
    do {
      unsigned int cur_nframes = std::min<unsigned int>
        ( first_only ? max_frames_per_chunk / 10 : max_frames_per_chunk,
          number_of_frames);
      // IMP_THREADS((sd, silent_statistics, cur_nframes),{
      std::cout << "Optimizing for " << cur_nframes << " frames in this iteration"
                << std::endl;
      sd->get_bd()->optimize(cur_nframes);
      print_score_and_positions(sd);
      //});
      if (sd->get_maximum_number_of_minutes() > 0 &&
          total_time.elapsed() / 60 > sd->get_maximum_number_of_minutes()) {
        sd->get_statistics()->set_interrupted(true);
        std::cout << "Terminating..." << std::endl;
        if (!silent_statistics) {
          sd->get_statistics()->update(timer, cur_nframes);
        }
        return false;
      }
      if (!silent_statistics) {
        sd->get_statistics()->update(timer, cur_nframes);
      }
      std::cout << "Done" << std::endl;
      number_of_frames -= cur_nframes;
    } while (number_of_frames > 0 && !first_only);
    return true;
  }

};
/****** END of anonymous namespace ********/

/********************* public functions **********************/

// initialize and return a simulation data object based on
// program command line parameters
IMP::npctransport::SimulationData *startup(int argc, char *argv[])
{
  IMP_NPC_PARSE_OPTIONS(argc, argv);
  IMP_ALWAYS_CHECK( short_init_factor > 0,
                    "short_init_factor must be positive",
                    IMP::ValueException );
  IMP_OMP_PRAGMA(critical)
    std::cout << "Random seed is " << IMP::get_random_seed() << std::endl;
  IMP::Pointer<IMP::npctransport::SimulationData> sd;
  write_output_based_on_flags(IMP::get_random_seed());
  sd = new IMP::npctransport::SimulationData(output, IMP::run_quick_test);
  if (!conformations.empty()) {
    sd->set_rmf_file(conformations,
                     !no_save_restraints_to_rmf);
  }
  if (!init_rmffile.empty()) {
    sd->initialize_positions_from_rmf
      ( RMF::open_rmf_file_read_only(init_rmffile), -1);
    IMP_OMP_PRAGMA(critical)
      std::cout
      << "Initialize coordinates from last frame of an existing RMF file "
      << init_rmffile << std::endl;
  }
  if (!restart_fgs_only.empty()) {
    std::cout
      << "Initialize coordinates of FGs only from output file "
      << restart_fgs_only << std::endl;
    restart_fgs_from_reference_output_file(sd, restart_fgs_only);
  }
  return sd.release();
}

//! remove nup42 and its anchors (also from obstacles)
void remove_Nup42(SimulationData* sd)
{
  // remove all diffusing Nup42 beads:
  {
    std::cout << "Removing Nup42 diffusing beads" << std::endl;
    ParticlesTemp beads=sd->get_beads();
    core::ParticleType Nup42_type("Nup42");
    for(unsigned int i=0; i<beads.size(); i++){
      core::Typed typed_pi(beads[i]);
      if(typed_pi.get_type() == Nup42_type){
        //        m->remove_particle(beads[i]->get_index());
        core::XYZ xyz(beads[i]);
        if(xyz.get_z()<500.0){
          xyz.set_z(xyz.get_z()+2000.0);
        }
      }
    }
  }
  // remove all obstacles close to center - these are Nup42 anchors:
  {
    ParticlesTemp obstacles= sd->get_obstacle_particles();
    for(unsigned int i=0; i<obstacles.size(); i++){
      core::XYZ xyz(obstacles[i]);
      double r2=std::pow(xyz.get_x(),2.0)+std::pow(xyz.get_y(),2.0);
      if(r2<200*200 && xyz.get_z()<500.0){
        std::cout << "Moving Nup42 obstacle at " << xyz.get_coordinates() << std::endl;
        //        m->remove_particle(obstacles[i]->get_index());
        xyz.set_z(xyz.get_z()+2000.0);
      }
    }
  }
}

//! inflate floater of specified type to new_radius
void inflate_floater
(SimulationData* sd, const std::string floater_name, const float new_radius){
  std::cout << "inflating " << floater_name << " to "
            << new_radius << " A" << std::endl;
  ParticlesTemp beads=sd->get_beads();
  core::XYZRs floaters;
  core::ParticleType floater_type(floater_name);
  double cur_radius(-1.0);
  for(unsigned int i=0; i<beads.size(); i++){
    core::Typed typed_pi(beads[i]);
    if(typed_pi.get_type() == floater_type){
      core::XYZR xyzr(beads[i]);
      floaters.push_back(xyzr);
      if(floaters.size()==1){
        cur_radius= xyzr.get_radius();
        std::cout << "Current radius = " << cur_radius << std::endl;
      } else {
        IMP_ALWAYS_CHECK(cur_radius==xyzr.get_radius(),
                         "inconsistent radii for same particle type",
                         IMP::UsageException);
      }
    }
  } // for i
  algebra::Sphere3Ds orig_sites=
    sd->get_sites(floater_type);
  algebra::Sphere3Ds inflated_sites=
    sd->get_sites(floater_type);
  ParticleTypeSet const& fg_bead_types = sd->get_fg_bead_types();
  std::cout << "Inflating begins" << std::endl;
  double delta;
  if(new_radius > cur_radius){
    delta= 0.05;
  }else{
    delta= -0.05;
  }
  for(double r=cur_radius;
      std::abs(r-new_radius)>0.5*std::abs(delta);
      r+= delta){
    std::cout << "Inflating sites to radius = " << r
              << " score "
              << sd->get_bd()->get_scoring_function()->evaluate(false)
              << std::endl;
    for(unsigned int j=0; j<inflated_sites.size(); j++){
      algebra::Vector3D orig_site=
        orig_sites[j].get_center();
      inflated_sites[j]._set_center(orig_site * r / cur_radius);
    }
    for(unsigned int i=0; i<floaters.size(); i++){
      floaters[i].set_radius(r);
      atom::RigidBodyDiffusion rbd(beads[i]);
      rbd.set_diffusion_coefficient // TODO: we assumed d_factor=1, need to adjust if it isn't
        (atom::get_einstein_diffusion_coefficient(r));
      rbd.set_rotational_diffusion_coefficient
        (atom::get_einstein_rotational_diffusion_coefficient(r)
         * sd->get_angular_d_factor());
      sd->set_sites(floater_type, inflated_sites);
      // Update all interactions involving type:
      // TODO: update score straight  from simulationData->set_sites() via Scoring,
      //       instead of from here
      for(ParticleTypeSet::const_iterator it=fg_bead_types.begin();
          it!=fg_bead_types.end(); it++){
        SitesPairScore* ps1=
          dynamic_cast<SitesPairScore*>(sd->get_scoring()->get_predicate_pair_score(floater_type, *it));
        IMP_ALWAYS_CHECK(ps1 != nullptr, "ps1 is not SitesPairScore", IMP::ValueException);
        ps1->set_sites0(inflated_sites);
        SitesPairScore* ps2=
          dynamic_cast<SitesPairScore*>(sd->get_scoring()->get_predicate_pair_score(*it, floater_type));
        IMP_ALWAYS_CHECK(ps2 != nullptr, "ps2 is not SitesPairScore", IMP::ValueException);
        ps1->set_sites1(inflated_sites);
      }
    } // i
    sd->get_bd()->optimize(50);
    // if(conformations_rmf_sos && (std::abs(r-std::round(r))<0.001)){
    //   conformations_rmf_sos->update_always("Inflating");
    // }
  } // r
  // Update output file (read, find floater, update radius):
  bool read(false);
  int fd= open(output.c_str(), O_RDONLY);
  google::protobuf::io::FileInputStream fis(fd);
  google::protobuf::io::CodedInputStream cis(&fis);
  cis.SetTotalBytesLimit(500000000,200000000);
  ::npctransport_proto::Output new_output;
  read= new_output.ParseFromCodedStream(&cis);
  close(fd);
  IMP_ALWAYS_CHECK(read, "Couldn't read output file " << output << " file descriptor " << fd,
                   IMP::ValueException);
  ::npctransport_proto::Assignment* pb_assignment = new_output.mutable_assignment();
  for (int i = 0; i < pb_assignment->floaters_size(); ++i) {
    ::npctransport_proto::Assignment_FloaterAssignment* f_data=
      pb_assignment->mutable_floaters(i);
    if(f_data->type() == floater_name){
      f_data->mutable_radius()->set_value(new_radius);
      break;
    }
  }
  std::ofstream outf(output.c_str(), std::ios::binary);
  new_output.SerializeToOstream(&outf);
}

void reset_box_size(SimulationData* sd, double box_size){
  sd->set_box_size(box_size);
  // Update output file:
  bool read(false);
  int fd= open(output.c_str(), O_RDONLY);
  google::protobuf::io::FileInputStream fis(fd);
  google::protobuf::io::CodedInputStream cis(&fis);
  cis.SetTotalBytesLimit(500000000,200000000);
  ::npctransport_proto::Output new_output;
  read= new_output.ParseFromCodedStream(&cis);
  close(fd);
  IMP_ALWAYS_CHECK(read, "Couldn't read output file " << output << " file descriptor " << fd,
                   IMP::ValueException);
  ::npctransport_proto::Assignment* pb_assignment = new_output.mutable_assignment();
  pb_assignment->mutable_box_side()->set_value(box_size);
  std::ofstream outf(output.c_str(), std::ios::binary);
  new_output.SerializeToOstream(&outf);
}

//!  Run simulation using preconstructed SimulationData object sd,
//!  with ad-hoc init restratins init_restraints
void do_main_loop(SimulationData *sd, const RestraintsTemp &init_restraints) {
  using namespace IMP;
  if(0){
    std::cout << "Coordinates of second fg bead in start of do_main_loop(): "
              << core::XYZ(sd->get_fg_chain_roots()[0].get_child(1))
              << std::endl;
  }

  sd->set_was_used( true );
  const int max_frames_per_chunk = sd->get_output_statistics_interval_frames();
  /** initial optimization and equilibration needed unless starting
      from another output file or rmf file */
  bool is_initial_optimization = (restart.empty() && init_rmffile.empty()) ||
    is_force_initialization_on_restart;
  bool is_BD_equilibration = is_initial_optimization;
  bool is_BD_full_run = !initialize_only;

  Pointer<rmf::SaveOptimizerState> conformations_rmf_sos =
      sd->get_rmf_sos_writer();
  RMF::FileHandle final_rmf_fh;
  if (!final_conformations.empty()) {
    final_rmf_fh = RMF::create_rmf_file(final_conformations);
    sd->link_rmf_file_handle(final_rmf_fh, !no_save_restraints_to_rmf);
  }
  boost::timer total_time;
  for (unsigned int i = 0; i < sd->get_number_of_trials(); ++i) {
    IMP::CreateLogContext clc("iteration");
    boost::timer timer;
    //    IMP::set_log_level(IMP::PROGRESS);
    std::cout << "Simulation trial " << i << " out of "
              << sd->get_number_of_trials() << std::endl;
    if (is_initial_optimization) {
      sd->switch_suspend_rmf(true);
      std::cout << "Doing initial coordinates optimization..." << std::endl;
      bool is_disable_randomize=
        is_force_initialization_on_restart || !restart_fgs_only.empty();
      initialize_positions(sd,
                           init_restraints,
                           verbose,
                           short_init_factor,
                           is_disable_randomize,
                           !restart_fgs_only.empty());
      if(0){
        std::cout << "Coordinates of second fg bead in start after initialization: "
                  << core::XYZ(sd->get_fg_chain_roots()[0].get_child(1))
                  << std::endl;
      }

      sd->get_bd()->set_current_time(0.0);
      sd->get_statistics()->reset_statistics_optimizer_states();
      sd->switch_suspend_rmf(false);
      //      conformation_rmf_sos->reset();
    }
    sd->activate_statistics();
    print_score_and_positions(sd, verbose, "Score right before BD = ");
    if (conformations_rmf_sos) {
      conformations_rmf_sos->update_always("Right before BD");
    }
    /*IMP::benchmark::Profiler p;
      if(i == 0)
      p.set("profiling.pprof");*/
    // sd->get_bd()->set_log_level(IMP::PROGRESS);
    unsigned int actual_nframes =
      (unsigned int)std::ceil(sd->get_number_of_frames() * short_sim_factor);
    unsigned int nframes_run = (unsigned int)(actual_nframes *
                                              sd->get_statistics_fraction());
    unsigned int nframes_equilibrate = actual_nframes - nframes_run;
    if (is_BD_equilibration) {
      std::cout << "Equilibrating for " << nframes_equilibrate << " frames..."
                << std::endl;
      bool ok = run_it(sd, nframes_equilibrate, timer, total_time,
                       true /* silent stats */, max_frames_per_chunk);
      if (!ok || first_only) return;
      // if(nframes_equilibrate > 0) {
      // if equilibrated, ignore equilibration stats
      // TODO: removed for now since this may be incosistent with
      //       consecutive runs
      // sd->reset_statistics_optimizer_states();
      // }
      sd->get_bd()->set_current_time(0.0);
      std::cout << "Equilibration finished succesfully" << std::endl;
    }
    if( is_inflate_kap28 ) {
      sd->switch_suspend_rmf(true);
      reset_box_size(sd, sd->get_box_size()*3.0);
      sd->get_bd()->set_scoring_function // update scoring function with new box size
      ( sd->get_scoring()->get_scoring_function(true) );
      remove_Nup42(sd);
      inflate_floater(sd, "kap16", 14.0);
      inflate_floater(sd, "kap18", 14.0);
      //      inflate_floater(sd, "kap20", 14.0);
      inflate_floater(sd, "kap22", 30.0);
      inflate_floater(sd, "kap24", 40.0);
      inflate_floater(sd, "kap26", 75.0);
      inflate_floater(sd, "kap28", 150.0);
      if(!conformations.empty()){
        sd->set_rmf_file("INFLATED_"+conformations,
                         !no_save_restraints_to_rmf);
      }
      sd->switch_suspend_rmf(false);
    }
    if (is_BD_full_run) {
      timer.restart();
      std::cout << "Running for " << nframes_run << " frames..." << std::endl;
      if (conformations_rmf_sos) {
        sd->get_bd()->get_scoring_function()->evaluate(false); // score before writing rmf so restraints and such are up to date
        conformations_rmf_sos->update_always(
            "Before running (post equilibration)");
      }
      // now start initial stats and run the rest of the sim
      //sd->get_statistics()->update(timer, 0);
      bool ok = run_it(sd, nframes_run, timer, total_time,
                       false /* silent stats */, max_frames_per_chunk);
      if (!ok) {
        return;
      }
      std::cout << "Run trial #" << i << " finished succesfully" << std::endl;
    }
    if (conformations_rmf_sos) {
      sd->get_bd()->get_scoring_function()->evaluate(false); // score before writing rmf so restraints and such are up to date
      conformations_rmf_sos->update_always("Final frame");
    }
    if (!final_conformations.empty()) {
      std::cout << "Printing last frame to " << final_conformations
                << std::endl;
      IMP::rmf::save_frame(final_rmf_fh);
    }
  }
  std::cout << "Entire run finished" << std::endl;
  print_score_and_positions(sd, verbose, "Final score = ");
}

IMPNPCTRANSPORT_END_NAMESPACE
