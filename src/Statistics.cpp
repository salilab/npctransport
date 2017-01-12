/**
 *  \file Statistics.cpp
 *  \brief statistics and order parameters about the simulations
 *         that is associated with a SimulationData object
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
*/

#include <IMP/npctransport/Statistics.h>
#include <IMP/npctransport/SimulationData.h>
#include <IMP/npctransport/FGChain.h>
#include <IMP/npctransport/protobuf.h>
#include <IMP/npctransport/enums.h>
#include <IMP/npctransport/io.h>
#include <IMP/npctransport/typedefs.h>
#include <IMP/npctransport/util.h>
#include <IMP/algebra/vector_generators.h>
#include <IMP/atom/estimates.h>
#include <IMP/atom/distance.h>
#include <IMP/atom/Diffusion.h>
#include <IMP/atom/Selection.h>
#include <IMP/log.h>
#include <IMP/check_macros.h>
#include <IMP/flags.h>
#include <IMP/core/pair_predicates.h>
#include <IMP/core/XYZR.h>
#include <IMP/core/generic.h>
#include <IMP/display/LogOptimizerState.h>
#include <IMP/rmf/frames.h>
#include <RMF/HDF5/File.h>
#include <RMF/HDF5/Group.h>
#include <RMF/HDF5/DataSetD.h>

#include <algorithm>
#include <numeric>
#include <set>
#include <math.h>
#include "boost/tuple/tuple.hpp"

#include <IMP/npctransport/internal/npctransport.pb.h>

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/coded_stream.h>
#include <fcntl.h>

// struct Int32TraitsBase : RMF::HDF5::IntTraitsBase {
//   // typedef int Type;
//   // typedef std::vector<int> Types;
//   // static const bool BatchOperations = true;
//   // static int get_index() { return 0; }
//   // static const Type& get_null_value() {
//   //   static Type null = std::numeric_limits<int>::max();
//   //   return null;
//   // }
//   // static bool get_is_null_value(Type t) { return t == get_null_value(); }
//   // static hid_t get_hdf5_fill_type() { return H5T_NATIVE_INT; }
//   static hid_t get_hdf5_disk_type() { return H5T_STD_I32LE; }
//   // static hid_t get_hdf5_memory_type() { return H5T_NATIVE_INT; }
//   // static const Type& get_fill_value() { return get_null_value(); }
//   // static std::string get_name() { return "int32"; }
// };

// struct Int32Traits : RMF::HDF5::SimpleTraits<Int32TraitsBase> {};

bool no_save_rmf_to_output = false;
IMP::AddBoolFlag  no_save_rmf_to_output_adder
( "no_save_rmf_to_output",
  "If true, save final rmf buffer to output file [default=false]",
  &no_save_rmf_to_output);

// TODO: turn into a template inline in unamed space?
/**
   updates (message).field() with a weighted average of its current
   value and new_value, giving weight n_old_frames, n_new_frames to each,
   respectively.
*/
#define UPDATE_AVG(n_frames, n_new_frames, message, field, new_value)   \
  if(n_new_frames>0) {                                                  \
    if((message).has_##field()){                                        \
      (message).set_##field(static_cast<double>(n_frames * (message).field() + \
                                                n_new_frames *new_value) /     \
                            (n_frames + n_new_frames));                 \
    } else {                                                            \
      (message).set_##field(new_value);                                 \
    }                                                                   \
  }



IMPNPCTRANSPORT_BEGIN_NAMESPACE

// ctr
Statistics::Statistics
( SimulationData* owner_sd,
  unsigned int statistics_interval_frames,
  std::string output_file_name)
: Object("Statistics%1%"),
  owner_sd_(owner_sd),
  statistics_interval_frames_(statistics_interval_frames),
  output_file_name_(output_file_name),
  is_stats_reset_(false)
{
}


//! add statistics about an FG chain
void Statistics::add_fg_chain_stats
( FGChain* fg_chain )
{
  ParticlesTemp chain_beads = fg_chain->get_beads();
  if(chain_beads.size() == 0)
    return;
  // chain stats
  IMP_NEW( ChainStatisticsOptimizerState, csos,
           (chain_beads, statistics_interval_frames_ ) );
  core::ParticleType chain_type = core::Typed(fg_chain->get_root()).get_type();
  chains_stats_map_[chain_type].push_back( csos );
  // stats for each chain particle
  std::set<core::ParticleType> p_types_encountered;
  for (unsigned int k = 0; k < chain_beads.size(); ++k) {
    Particle* p = chain_beads[k];
    core::ParticleType p_type = core::Typed(p).get_type();
    //    IMP_ALWAYS_CHECK( chain_type == p_type,
    //                  "All beads in chain must be of same type for now",
    //                  ValueException);
    if(p_types_encountered.find(p_type) == p_types_encountered.end()) {
      p_types_encountered.insert(p_type);
      fgs_bodies_stats_map_[p_type].push_back( BodyStatisticsOptimizerStates() );
    }
    IMP_NEW(BodyStatisticsOptimizerState, bsos,
            ( p, this,  statistics_interval_frames_ ) );
    fgs_bodies_stats_map_[p_type].back().push_back( bsos );
  }  // for k
}

//! add statistics about a floater particle
void Statistics::add_floater_stats
( Particle* p )
{
  core::ParticleType type = core::Typed(p).get_type();
  IMP_NEW(BodyStatisticsOptimizerState, bsos,
          (p, this, statistics_interval_frames_));
  floaters_stats_map_[type].push_back(bsos);
  if (get_sd()->get_has_slab() )
    {  // only if has tunnel
      IMP_NEW(ParticleTransportStatisticsOptimizerState, ptsos,
              (p,
               -0.5 * get_sd()->get_slab_thickness(),  // tunnel bottom
               0.5 * get_sd()->get_slab_thickness(),  // tunnel top
               this // statistics manager
               )
              );
      ptsos->set_period(statistics_interval_frames_);
      floaters_transport_stats_map_[type].push_back( ptsos );
    }
}

void Statistics::add_interaction_stats
( core::ParticleType type0, core::ParticleType type1)
{
    InteractionType interaction_type = std::make_pair(type0, type1);
    IMP_USAGE_CHECK( interaction_stats_map_.find( interaction_type ) ==
                     interaction_stats_map_.end(),
                     "cannot add stats about the same interaction type twice" );

    // add statistics about this interaction to interactions_stats_
    // between all diffusing particles
    ParticlesTemp set0, set1; // TODO: turn to a real set?!
    Particles& beads = get_sd()->get_beads_byref();
    for(unsigned int i = 0; i < beads.size(); i++){
          if (core::Typed(beads[i]).get_type() == type0)
            { set0.push_back( beads[i] ); }
          if (core::Typed(beads[i]).get_type() == type1)
            { set1.push_back( beads[i] ); }
    }
    bool include_site_site = true;
    bool include_non_specific = true;
    double range =
      get_sd()->get_scoring()->get_interaction_range_for
      ( type0, type1, include_site_site, include_non_specific);
    double slack=3.0; // TODO: param
    IMP_LOG(PROGRESS,
            "Interaction " << type0.get_string() << ", " << type1.get_string()
            << "  sizes: " << set0.size() << ", " << set1.size()
            << " statistics range: " << range << std::endl );
    if (set0.size() > 0 && set1.size() > 0)
      {
        IMP_NEW(BipartitePairsStatisticsOptimizerState, bpsos,
                (this, interaction_type, set0, set1,
                 range, slack));
        bpsos->set_period(statistics_interval_frames_);
        interaction_stats_map_[interaction_type] = bpsos;
      }
}

// get all statistics periodic optimizer states in one list
OptimizerStates Statistics::add_optimizer_states(Optimizer* o)
{
  if(o == nullptr) o = get_sd()->get_bd();
  IMP_ALWAYS_CHECK( o, "null optimizer in add_optimizer_states()"
                    " and get_sd()->get_bd() is invalid",
                    ValueException);
  OptimizerStates ret;
  for (FGsBodyStatisticsOSsMap::iterator iter = fgs_bodies_stats_map_.begin();
       iter != fgs_bodies_stats_map_.end(); iter++)
    {
      for (unsigned int j = 0; j < iter->second.size(); j++)
        {
          ret += iter->second[j];
        } // for j
    } // for iter
  for (BodyStatisticsOSsMap::iterator iter = floaters_stats_map_.begin();
       iter != floaters_stats_map_.end(); iter++)
    {
      ret += iter->second;
    } // for iter
  if ( get_sd()->get_has_slab() )
    {
      for (ParticleTransportStatisticsOSsMap::iterator
             iter = floaters_transport_stats_map_.begin();
           iter != floaters_transport_stats_map_.end(); iter++)
        {
          ret += iter->second;
          for(unsigned int j = 0; j < iter->second.size(); j++) {
            // TODO: this is problematic encapsulation wise
            //       perhaps needs to provide 'owner' as parameter,
            //       with default being get_sd()->get_bd()
            iter->second[j]->set_owner( get_sd()->get_bd() );
          } // for j
        } // for iter
    }
  for (ChainStatisticsOSsMap::iterator iter = chains_stats_map_.begin();
       iter != chains_stats_map_.end(); iter++)
    {
      ret += iter->second;
    } // for iter
  for (BipartitePairsStatisticsOSMap::iterator
         iter = interaction_stats_map_.begin();
       iter != interaction_stats_map_.end(); iter++)
    {
      ret.push_back( iter->second ) ;
    }
  o->add_optimizer_states(ret);
  return ret;
}


bool
Statistics::update_xyz_distribution_to_hdf5
(RMF::HDF5::Group hdf5_group,
 core::ParticleType p_type){
  // Recreate z-y-x histogram based on retrieved xyz_hist:
  ParticleTypeXYZDistributionMap::const_iterator ptxyzdm_it=
    particle_type_xyz_distribution_map_.find(p_type);
  if(ptxyzdm_it != particle_type_xyz_distribution_map_.end()) {
    ParticleTypeXYZDistributionMap::mapped_type xyz_hist=
      ptxyzdm_it->second;
    // retrieve or create dataset in hdf5
    RMF::HDF5::DataSetD<RMF::HDF5::IntTraits, 3> ds_xyz;
    std::string s_type= p_type.get_string();
    //    std::cout << "Outputing stats for " << s_type << std::endl;
    RMF::HDF5::IntTraits::Type fill_value(0);
    if(hdf5_group.get_has_child(s_type)) {
      ds_xyz= hdf5_group.get_child_data_set
        < RMF::HDF5::IntTraits,3 > (s_type);
      std::cout << "Warning - hist for " << s_type
                << " already exists but it was not expected to" << std::endl;
    }else{
      RMF::HDF5::DataSetCreationPropertiesD
        < RMF::HDF5::IntTraits,3 > dscp;
      RMF_HDF5_CALL(H5Pset_deflate(dscp.get_handle(), 1)); // compression
      dscp.set_custom_fill_value(&fill_value);
      ds_xyz= hdf5_group.add_child_data_set
        < RMF::HDF5::IntTraits,3 > (s_type, dscp);
    }
    // set size (override any existing settings):
    boost::uint_fast8_t& d0= xyz_distribution_sizes_.d0;
    boost::uint_fast8_t& d1= xyz_distribution_sizes_.d1;
    boost::uint_fast8_t& d2= xyz_distribution_sizes_.d2;
    ds_xyz.set_size(RMF::HDF5::DataSetIndexD<3>(d0,d1,d2));
    // set non-zero values:
    for(t_sparse_3d_matrix::iterator iter_ii= xyz_hist.begin();
        iter_ii != xyz_hist.end(); iter_ii++){
      boost::uint_fast8_t ii= iter_ii->first;
      for(t_sparse_3d_matrix::mapped_type::iterator iter_jj= iter_ii->second.begin();
          iter_jj != iter_ii->second.end(); iter_jj++){
        boost::uint_fast8_t jj= iter_jj->first;
        for(t_sparse_3d_matrix::mapped_type::mapped_type::iterator iter_kk= iter_jj->second.begin();
            iter_kk != iter_jj->second.end(); iter_kk++){
          boost::uint_fast8_t kk= iter_kk->first;
          RMF::HDF5::DataSetIndexD<3> iijjkk(ii,jj,kk);
          ds_xyz.set_value(iijjkk, iter_kk->second);
        } // iter_kk
      } // iter_jj
    } // iter_ii
    //    std::cout << "Finished outputing stats for " << s_type << std::endl;
    return true;
  } //if ptxyzdm_it
  return false;
}

void Statistics::update_fg_stats
( ::npctransport_proto::Statistics* stats,
  unsigned int nf_new,
  unsigned int zr_hist[4][3],
  RMF::HDF5::File hdf5_file)
{
  RMF::HDF5::Group hdf5_fg_xyz_hist_group;
  static const std::string  FG_XYZ_GROUP("fg_xyz_hist");
  if(hdf5_file.get_has_child(FG_XYZ_GROUP)){
    IMP_ALWAYS_CHECK(hdf5_file.get_child_is_group(FG_XYZ_GROUP),
                     FG_XYZ_GROUP << " is supposed to be an HDF5 group",
                     ValueException);
    hdf5_fg_xyz_hist_group= hdf5_file.get_child_group(FG_XYZ_GROUP);
  } else {
    hdf5_fg_xyz_hist_group= hdf5_file.add_child_group(FG_XYZ_GROUP);
  }
  int nf = stats->number_of_frames();
  double sim_time_ns = const_cast<SimulationData *>( get_sd() )
    ->get_bd()->get_current_time() / FS_IN_NS;

  // Go over all fg types:
  ParticleTypeSet const &fgt = get_sd()->get_fg_types();
  for ( ParticleTypeSet::const_iterator
          it = fgt.begin(); it != fgt.end(); it++ )
    {
      core::ParticleType type_i = *it;
      atom::Hierarchy root_i = get_sd()->get_root_of_type( type_i );
      ParticlesTemp chains_i = root_i.get_children();
      if(chains_i.size() == 0)
        continue;
      unsigned int i = find_or_add_fg_of_type( stats, type_i );
      npctransport_proto::Statistics_FGOrderParams*
        fgi_op = stats->mutable_fgs(i)->add_order_params();
      fgi_op->set_time_ns(sim_time_ns);

      // Average chain stats from optimizer state:
      {
        IMP_USAGE_CHECK(chains_stats_map_.find(type_i)
                        != chains_stats_map_.end(),
                        "type missing from stats");
        ChainStatisticsOptimizerStates& cs_i =
          chains_stats_map_.find(type_i)->second;
        double mean_radius_of_gyration=0.0;
        double mean_square_radius_of_gyration=0.0;
        double mean_end_to_end_distance=0.0;
        double mean_square_end_to_end_distance=0.0;
        for (unsigned int j = 0; j < cs_i.size(); ++j)
          {
            unsigned int cnf = nf * cs_i.size() + j;
            UPDATE_AVG(cnf, nf_new,
                       *stats->mutable_fgs(i), chain_correlation_time,
                       cs_i[j]->get_correlation_time());
            UPDATE_AVG(cnf, nf_new, *stats->mutable_fgs(i),
                       chain_diffusion_coefficient,
                       cs_i[j]->get_diffusion_coefficient());
            Floats df = cs_i[j]->get_diffusion_coefficients();
            UPDATE_AVG(cnf, nf_new, *stats->mutable_fgs(i),
                       local_diffusion_coefficient,
                       std::accumulate(df.begin(), df.end(), 0.0) / df.size());
            mean_radius_of_gyration+=
              cs_i[j]->get_mean_radius_of_gyration()/cs_i.size();
            mean_square_radius_of_gyration+=
              cs_i[j]->get_mean_square_radius_of_gyration()/cs_i.size();
            mean_end_to_end_distance+=
              cs_i[j]->get_mean_end_to_end_distance()/cs_i.size();
            mean_square_end_to_end_distance+=
              cs_i[j]->get_mean_square_end_to_end_distance()/cs_i.size();
            cs_i[j]->reset();
          } // for j (fg chain)
        UPDATE_AVG(nf, nf_new, *stats->mutable_fgs(i),
                   radius_of_gyration, mean_radius_of_gyration);
        UPDATE_AVG(nf, nf_new, *stats->mutable_fgs(i),
                   length, mean_end_to_end_distance);
        fgi_op->set_mean_radius_of_gyration
          (mean_radius_of_gyration);
        fgi_op->set_mean_square_radius_of_gyration
          (mean_square_radius_of_gyration);
        fgi_op->set_mean_end_to_end_distance
          (mean_end_to_end_distance);
        fgi_op->set_mean_square_end_to_end_distance
          (mean_square_end_to_end_distance);
      }

      // Additional chain stats from current snapshot:
      {
        double avg_volume = 0.0;
        for (unsigned int j = 0; j < chains_i.size(); ++j)
          {
            Pointer<FGChain> chain_ij= get_fg_chain(chains_i[j]);
            fill_in_zr_hist(zr_hist, chain_ij->get_beads());
#ifdef IMP_NPCTRANSPORT_USE_IMP_CGAL
            double volume_ij =
              atom::get_volume(chain_ij->get_root()); // how does the work with TAMD?
#else
            double volume_ij = -1.;
#endif
            UPDATE_AVG(nf, nf_new, *stats->mutable_fgs(i), volume, volume_ij);
            avg_volume += volume_ij / chains_i.size();
          } // for j (fg chain)
        fgi_op->set_volume(avg_volume);
      }

      // FG average body stats from optimizer state:
      {
        IMP_USAGE_CHECK(fgs_bodies_stats_map_.find(type_i) !=
                        fgs_bodies_stats_map_.end(),
                        "type missing from stats");
        FGsBodyStatisticsOSs& fbs_i = fgs_bodies_stats_map_.find(type_i)->second;
        for (unsigned int j = 0; j < fbs_i.size(); ++j)
          {
            BodyStatisticsOptimizerStates& fbs_ij = fbs_i[j];
            for (unsigned int k = 0; k < fbs_ij.size(); ++k)
              {
                BodyStatisticsOptimizerState* fbs_ijk = fbs_ij[k];
                unsigned int per_frame = fbs_i.size() * fbs_ij.size();
                unsigned int cnf = (nf) * per_frame + j * fbs_ij.size() + k;
                UPDATE_AVG(cnf, nf_new, *stats->mutable_fgs(i),
                           particle_correlation_time,
                           fbs_ijk->get_correlation_time());
                UPDATE_AVG(cnf, nf_new, *stats->mutable_fgs(i),
                           particle_diffusion_coefficient,
                           fbs_ijk->get_diffusion_coefficient());
                fbs_ijk->reset();
              } // for k
          } // for j
      }

      // Recreate z-r histogram based on retrieved zr_hist:
      if(get_sd()->get_is_xyz_hist_stats()){
        stats->mutable_fgs(i)->clear_xyz_hist();
        update_xyz_distribution_to_hdf5(hdf5_fg_xyz_hist_group,
                                        type_i);
      } else {
        ParticleTypeZRDistributionMap::const_iterator ptzrdm_it=
          particle_type_zr_distribution_map_.find(type_i);
        if(ptzrdm_it != particle_type_zr_distribution_map_.end())
          {
            ParticleTypeZRDistributionMap::mapped_type zr_hist=
              ptzrdm_it->second;
            stats->mutable_fgs(i)->clear_zr_hist();
            for(unsigned int ii=0; ii < zr_hist.size(); ii++)
              {
                ::npctransport_proto::Statistics_Ints* zii_r_hist=
                  stats->mutable_fgs(i)->mutable_zr_hist()->add_ints_list();
                for(unsigned int jj=0; jj < zr_hist[ii].size(); jj++)
                  {
                    zii_r_hist->add_ints(zr_hist[ii][jj]);
                  } // for jj
              } // for ii
          } // if ptzed_it
      } // if is_xyz_hist
    } // for it (fg type)
}

void Statistics
::update_particle_type_zr_distribution_map
( Particle* p )
{
  if ( !get_sd()->get_has_slab() || !get_sd()->get_has_bounding_box() ){
    return;
  }
  const float GRID_RESOLUTION_ANGSTROMS= 10.0; // resolution of zr grid
  bool is_z_symmetric=
    (get_sd()->get_output_npctransport_version() < 2.0);
  core::ParticleType pt( core::Typed(p).get_type() );
  std::pair<ParticleTypeZRDistributionMap::iterator, bool> it_pair;
  it_pair.first= particle_type_zr_distribution_map_.find(pt);
  // add distribution table if needed
  if(it_pair.first==particle_type_zr_distribution_map_.end()) {
    float z_max =  get_sd()->get_box_size() / 2.0; // get_z_distribution_top();
    float r_max =  get_sd()->get_box_size() / std::sqrt(2.0); // get_r_distribution_max();
    unsigned int nz=
      (1 + !is_z_symmetric) * (std::floor(z_max/GRID_RESOLUTION_ANGSTROMS) + 5); // +5 for slack
    unsigned int nr=
      std::floor(r_max/GRID_RESOLUTION_ANGSTROMS)+5; // +5 for slack
    ParticleTypeZRDistributionMap::value_type vt =
      std::make_pair(pt,
                     std::vector<std::vector<int> >
                     (nz, std::vector<int>(nr, 0)));
    it_pair=particle_type_zr_distribution_map_.insert(vt);
    IMP_USAGE_CHECK(it_pair.second, "type " << pt << " already exists");
  }
  ParticleTypeZRDistributionMap::iterator& it=it_pair.first;
  //update distribution
  core::XYZ xyz(p);
  float x= xyz.get_x();
  float y= xyz.get_x();
  float z= xyz.get_z();
  if(is_z_symmetric){
    z= std::abs(z);
  }
  float r= std::sqrt( x*x+y*y);
  unsigned int nz= it->second.size();
  unsigned int zz= std::floor(z/GRID_RESOLUTION_ANGSTROMS)
    + (nz/2)*(!is_z_symmetric);
  unsigned int rr= std::floor(r/GRID_RESOLUTION_ANGSTROMS);
  it->second[zz][rr]++;
}

void Statistics
::update_particle_type_xyz_distribution_map
( Particle* p )
{
  if ( !get_sd()->get_has_slab() || !get_sd()->get_has_bounding_box() ){
    return;
  }
  const float GRID_RESOLUTION_ANGSTROMS=10; // resolution of zr grid
  const float CROP_FACTOR=0.5; // crop 0.5*Crop x 2 on each dimension (e.g. for box size of 200, only include -50 to +50 and not -100 to +100 on each dimension
  const double MAX_CROP=1000.0;
  bool is_z_symmetric=
    (get_sd()->get_output_npctransport_version() < 2.0);
  core::ParticleType pt( core::Typed(p).get_type() );
  // retrieve, add distribution table if needed
  std::pair<ParticleTypeXYZDistributionMap::iterator, bool> it_pair;
  it_pair.first= particle_type_xyz_distribution_map_.find(pt);
  if(it_pair.first==particle_type_xyz_distribution_map_.end()) {
    ParticleTypeXYZDistributionMap::value_type vt =
      std::make_pair
      (pt, t_sparse_3d_matrix());
    it_pair=particle_type_xyz_distribution_map_.insert(vt);
    IMP_USAGE_CHECK(it_pair.second, "type " << pt << " already exists");
  }
  ParticleTypeXYZDistributionMap::iterator& it=it_pair.first;
  //update distribution
  float box_half =  std::min(get_sd()->get_box_size() / 2.0, MAX_CROP); // get_z_distribution_top();
  unsigned int half_n_max= std::floor(box_half/GRID_RESOLUTION_ANGSTROMS*CROP_FACTOR) + 5; // +5 for slack
  unsigned int nx= 2 * half_n_max;
  unsigned int ny= 2 * half_n_max;
  unsigned int nz= (1 + !is_z_symmetric) * half_n_max;
  xyz_distribution_sizes_.d0= nx; // TODO: cap with maximal value of d0, which is (or used to be) uint8
  xyz_distribution_sizes_.d1= ny;
  xyz_distribution_sizes_.d2= nz;
  core::XYZ xyz(p);
  float x= xyz.get_x();
  float y= xyz.get_y();
  float z= xyz.get_z();
  if(is_z_symmetric){
    z= std::abs(z);
  }
  unsigned int xx=
    std::floor(x/GRID_RESOLUTION_ANGSTROMS) + nx/2;
  unsigned int yy=
    std::floor(y/GRID_RESOLUTION_ANGSTROMS) + ny/2;
  unsigned int zz=
    std::floor(z/GRID_RESOLUTION_ANGSTROMS) + (nz/2)*(!is_z_symmetric);
  if(xx<nx && yy < ny && zz<nz){
    it->second[xx][yy][zz]++;
    //    std::cout<<"xx "<<xx<<" yy "<<yy<<" zz "<<zz<<" = "<<it->second[xx][yy][zz]<<std::endl;
  }
}



// @param nf_new number of new frames accounted for in this statistics update
void Statistics::update
( const boost::timer &timer,
  unsigned int nf_new)
{
  IMP_OBJECT_LOG;
  ::npctransport_proto::Output output;

  //  std::ifstream inf(output_file_name_.c_str(), std::ios::binary);
  //  output.ParseFromIstream(&inf);
  //  inf.close();
  bool read(false);
  int fd=open(output_file_name_.c_str(), O_RDONLY);
  if(fd!=-1){
    google::protobuf::io::FileInputStream fis(fd);
    google::protobuf::io::CodedInputStream cis(&fis);
    cis.SetTotalBytesLimit(500000000,200000000);
    read=output.ParseFromCodedStream(&cis);
    close(fd);
  }
  IMP_ALWAYS_CHECK(read,
                   "Failed updating statistics to " << output_file_name_.c_str() << std::endl,
                   IMP::IOException);
  RMF::HDF5::File hdf5_file= RMF::HDF5::create_file(output_file_name_ + ".hdf5");
  RMF::HDF5::Group hdf5_floater_xyz_hist_group;
  static const std::string  FLOATER_XYZ_GROUP("floater_xyz_hist");
  if(hdf5_file.get_has_child(FLOATER_XYZ_GROUP)){
    IMP_ALWAYS_CHECK(hdf5_file.get_child_is_group(FLOATER_XYZ_GROUP),
                     FLOATER_XYZ_GROUP << " is supposed to be an HDF5 group",
                     ValueException);
    hdf5_floater_xyz_hist_group= hdf5_file.get_child_group(FLOATER_XYZ_GROUP);
  } else {
    hdf5_floater_xyz_hist_group= hdf5_file.add_child_group(FLOATER_XYZ_GROUP);
  }


  ::npctransport_proto::Statistics* stats = output.mutable_statistics();
  int nf = stats->number_of_frames();
  if (is_stats_reset_) {  // stats was just reset
    // TODO: what's if multiple trials?
    nf = 0;
    is_stats_reset_ = false;
    for (int i = 0; i < stats->floaters().size(); i++) {
      (*stats->mutable_floaters(i)).clear_transport_time_points_ns();
    }
  }
  IMP_LOG(VERBOSE, "Updating statistics file " << output_file_name_
            << " that currently has " << nf << " frames, with " << nf_new
            << " additional frames" << std::endl);

  // gather the statistics one by one
  double sim_time_ns = const_cast<SimulationData *>( get_sd() )
    ->get_bd()->get_current_time() / FS_IN_NS;

  unsigned int zr_hist[4][3]={{0},{0},{0},{0}};
  update_fg_stats(stats, nf_new, zr_hist, hdf5_file);

  std::map<IMP::core::ParticleType, double> type_to_diffusion_coefficeint_map; // to be used for floater order params
  // Floaters general body stats
  for (BodyStatisticsOSsMap::iterator
         it = floaters_stats_map_.begin() ;
       it != floaters_stats_map_.end(); it++)
    {
      unsigned int i= find_or_add_floater_of_type( stats, it->first );
      BodyStatisticsOptimizerStates& bsos = it->second;
      unsigned int n_particles_type_i = bsos.size();
      type_to_diffusion_coefficeint_map[it->first]=0.0;
      int nf_weighted = nf * n_particles_type_i; // number of particle frames
      for (unsigned int j= 0; j < n_particles_type_i; j++)
        {
          double dc_j= bsos[j]->get_diffusion_coefficient();
          UPDATE_AVG(nf_weighted, nf_new,
                     *stats->mutable_floaters(i),
                     diffusion_coefficient, dc_j);
          type_to_diffusion_coefficeint_map[it->first]+=
            dc_j/n_particles_type_i;
          double ct_j= bsos[j]->get_correlation_time();
          UPDATE_AVG(nf_weighted, nf_new,
                     *stats->mutable_floaters(i),
                     correlation_time, ct_j);
          bsos[j]->reset();
          nf_weighted++;
        } // for j
      if(get_sd()->get_is_xyz_hist_stats()){ // TODO: floaters are disabled for xyz for now to save space - perhaps add it later
        update_xyz_distribution_to_hdf5(hdf5_floater_xyz_hist_group,
                                        it->first);
      } else { // if/else is_xyz_hist_stats
        // Recreate z-r histogram based on retrieved zr_hist:
        ParticleTypeZRDistributionMap::const_iterator ptzrdm_it=
          particle_type_zr_distribution_map_.find(it->first);
        if(ptzrdm_it != particle_type_zr_distribution_map_.end()) {
          ParticleTypeZRDistributionMap::mapped_type zr_hist=
            ptzrdm_it->second;
          stats->mutable_floaters(i)->clear_zr_hist();
          for(unsigned int ii=0; ii < zr_hist.size(); ii++) {
            ::npctransport_proto::Statistics_Ints* zii_r_hist=
              stats->mutable_floaters(i)->mutable_zr_hist()->add_ints_list();
            for(unsigned int jj=0; jj < zr_hist[ii].size(); jj++) {
              zii_r_hist->add_ints(zr_hist[ii][jj]);
            } // for jj
          } // for ii
        } //if ptzrdm_it
      } // if/else is_xyz_hist_stats
    } // for it

  // Floaters avg number of transports per particle
  if ( get_sd()->get_has_slab() ){
    for (ParticleTransportStatisticsOSsMap::iterator
           it1 = floaters_transport_stats_map_.begin();
         it1 != floaters_transport_stats_map_.end() ; it1++)
      {
        unsigned int i = find_or_add_floater_of_type( stats, it1->first );
        ParticleTransportStatisticsOptimizerStates& pts_i = it1->second;
        // fetch old ones from stats msg, add new ones and rewrite all:
        std::set<double> times_i
          ( stats->floaters(i).transport_time_points_ns().begin(),
            stats->floaters(i).transport_time_points_ns().end() );
        for (unsigned int j = 0; j < pts_i.size(); ++j)
          {
            Floats const &new_times_ij =
              pts_i[j]->get_transport_time_points_in_ns();
            for (unsigned int k = 0; k < new_times_ij.size(); k++)
              {
                times_i.insert(new_times_ij[k]);
              } // for k
          } // for j
        (*stats->mutable_floaters(i)).clear_transport_time_points_ns();
        for (std::set<double>::const_iterator it2 = times_i.begin();
             it2 != times_i.end(); it2++)
          {
            (*stats->mutable_floaters(i)).add_transport_time_points_ns(*it2);
          } // for it2
         // update avg too:
         double avg_n_transports_i = times_i.size() * 1.0 / pts_i.size();
         (*stats->mutable_floaters(i)).set_avg_n_transports( avg_n_transports_i );
       } // for it1
   }

   // FG-floaters interactions
   ParticleTypeSet const &ft = get_sd()->get_floater_types();
   atom::Hierarchies fg_chain_roots =  get_sd()->get_fg_chain_roots() ;
   for ( ParticleTypeSet::const_iterator
           it = ft.begin(); it != ft.end(); it++ )
       {
         ParticlesTemp cur_floaters =
           get_sd()->get_root_of_type( *it ).get_children();
         if(cur_floaters.size() == 0) // TODO: superfluous?
           continue;
         unsigned int i = find_or_add_floater_of_type( stats, *it );
         double site_site_pairs;
         double interacting_floaters;
         double bead_floater_pairs;
         double chain_floater_pairs;
         boost::tie(site_site_pairs, interacting_floaters,
                    bead_floater_pairs, chain_floater_pairs)
           = get_interactions_and_interacting(cur_floaters, fg_chain_roots);
        npctransport_proto::Statistics_FloaterOrderParams*
          frc = stats->mutable_floaters(i)->add_order_params();
        frc->set_time_ns(sim_time_ns);
        frc->set_site_interactions_per_floater
          (site_site_pairs / cur_floaters.size());
        frc->set_interacting_fraction
          (interacting_floaters / cur_floaters.size());
        frc->set_beads_per_interacting_floater
          (bead_floater_pairs / interacting_floaters);
        frc->set_chains_per_interacting_floater
          (chain_floater_pairs / interacting_floaters);
        if(get_sd()->get_has_slab())
          {
            int n_z0, n_z1, n_z2, n_z3;
            boost::tie(n_z0, n_z1, n_z2, n_z3) =
              get_z_distribution(cur_floaters);
            frc->set_n_z0(n_z0);
            frc->set_n_z1(n_z1);
            frc->set_n_z2(n_z2);
            frc->set_n_z3(n_z3);
          }
        frc->set_diffusion_coefficient
          (type_to_diffusion_coefficeint_map[*it]);
      } // for(it)

  // update statistics gathered on interaction rates
  for (BipartitePairsStatisticsOSMap::iterator
         it = interaction_stats_map_.begin();
       it != interaction_stats_map_.end(); it++)
    {
      IMP_LOG(PROGRESS, "adding interaction statistics "
              << it->first << std::endl);
      unsigned int i = find_or_add_interaction_of_type( stats, it->first );
      BipartitePairsStatisticsOptimizerState* bps_i = it->second;
      ::npctransport_proto::Statistics_InteractionStats *pOutStats_i =
          stats->mutable_interactions(i);
      // verify correct interaction type is stored
      InteractionType itype = bps_i->get_interaction_type();
      std::string s_type0 = itype.first.get_string();
      std::string s_type1 = itype.second.get_string();
    if (std::string(pOutStats_i->type0()) != s_type0 ||
        std::string(pOutStats_i->type1()) != s_type1) {
      IMP_THROW("Incompatible interaction types in pOutStats_i ["
                    << pOutStats_i->type0() << ", " << pOutStats_i->type1()
                    << "] and bps_i " << s_type0 << ", " << s_type1 << "]"
                    << std::endl,
                ValueException);
    }
    // some preparations
    Int n0 = bps_i->get_number_of_particles_1();
    Int n1 = bps_i->get_number_of_particles_2();
    Float avg_contacts_num = bps_i->get_average_number_of_contacts();
    // save the rest of the interactions info
    npctransport_proto::Statistics_InteractionOrderParams*
      siop = pOutStats_i->add_order_params();
    siop->set_time_ns(sim_time_ns);
    siop->set_avg_off_per_contact_per_ns
      ( bps_i->get_average_off_per_contact_per_ns() );
    siop->set_avg_off_per_bound_i_per_ns
      ( bps_i->get_average_off_per_bound_I_per_ns() );
    siop->set_avg_off_per_bound_ii_per_ns
      ( bps_i->get_average_off_per_bound_II_per_ns() );
    siop->set_off_stats_period_ns
      ( bps_i->get_off_stats_period_ns() );
    siop->set_off_i_stats_period_ns
      ( bps_i->get_off_I_stats_period_ns() );
    siop->set_off_ii_stats_period_ns
      ( bps_i->get_off_II_stats_period_ns() );
    siop->set_avg_on_per_missing_contact_per_ns
      ( bps_i->get_average_on_per_missing_contact_per_ns() );
    siop->set_on_stats_period_ns
      ( bps_i->get_on_stats_period_ns() );
    siop->set_avg_on_per_unbound_i_per_ns
      ( bps_i->get_average_on_per_unbound_I_per_ns() );
    siop->set_on_i_stats_period_ns
      ( bps_i->get_on_I_stats_period_ns() );
    siop->set_avg_on_per_unbound_ii_per_ns
      ( bps_i->get_average_on_per_unbound_II_per_ns() );
    siop->set_on_ii_stats_period_ns
      ( bps_i->get_on_II_stats_period_ns() );
    siop->set_avg_contacts_per_particle_i
      ( avg_contacts_num / n0 );
    siop->set_avg_contacts_per_particle_ii
      ( avg_contacts_num / n1 );
    siop->set_avg_fraction_bound_particles_i
      ( bps_i->get_average_percentage_bound_particles_1());
    siop->set_avg_fraction_bound_particles_ii
      ( bps_i->get_average_percentage_bound_particles_2());
    siop->set_misc_stats_period_ns
      ( bps_i->get_misc_stats_period_ns() );
    // reset till next udpate_statistics()
    bps_i->reset();
  }

  // GLOBALS:
  {
    // Todo: define better what we want of timer
    stats->set_seconds_per_iteration(timer.elapsed());
    stats->set_number_of_frames(nf + nf_new);
    stats->set_bd_simulation_time_ns( sim_time_ns );
    double total_energy  =
      get_sd()->get_bd()->get_scoring_function()->evaluate(false);
    double energy_per_bead =
      total_energy / get_sd()->get_beads().size();
    UPDATE_AVG(nf, nf_new, (*stats), energy_per_particle,  // TODO: reset?
               // TODO: remove static beads from stats?
               energy_per_bead );
    ::npctransport_proto::Statistics_GlobalOrderParams*
        sgop = stats->add_global_order_params();
    sgop->set_time_ns(sim_time_ns);
    sgop->set_energy(total_energy);
    if(get_sd()->get_has_slab()){
      for(int zz=0; zz < 4; zz++)
        {
          ::npctransport_proto::Statistics_Ints* sis=
            sgop->add_zr_hists();
          for(int rr=0; rr < 3; rr++)
            {
              sis->add_ints(zr_hist[zz][rr]);
            }
        }
    }
  }

  // TODO: disable this for now
  //  ::npctransport_proto::Conformation *conformation =
  //   output.mutable_conformation();
  //    save_pb_conformation(get_beads(), sites_, conformation);

  // save RMF for future restarts
  if(!no_save_rmf_to_output){
    RMF::BufferHandle buf;
    {
      RMF::FileHandle fh = RMF::create_rmf_buffer(buf);
      const_cast<SimulationData*>(get_sd())->link_rmf_file_handle(fh, false);
      rmf::save_frame(fh);
    }
    output.set_rmf_conformation(buf.get_string());
  }

  // dump to file
  std::ofstream outf(output_file_name_.c_str(), std::ios::binary);
  output.SerializeToOstream(&outf);
  outf.flush();
  }

void Statistics::reset_statistics_optimizer_states()
{
  is_stats_reset_ = true;  // indicate to update()

  for (FGsBodyStatisticsOSsMap::iterator iter = fgs_bodies_stats_map_.begin();
       iter != fgs_bodies_stats_map_.end(); iter++)
    {
      for (unsigned int j = 0; j < iter->second.size(); j++)
        {
          for (unsigned int k = 0; k < iter->second[j].size(); ++k) {
            iter->second[j][k]->reset();
          } // for k
        } // for j
    } // for iter

  for (BodyStatisticsOSsMap::iterator iter = floaters_stats_map_.begin();
       iter != floaters_stats_map_.end(); iter++)
    {
      for (unsigned int j = 0; j < iter->second.size(); j++)
        {
          iter->second[j]->reset();
        } // for jb
    } // for iter

  if (get_sd()->get_has_slab()) {
    for (ParticleTransportStatisticsOSsMap::iterator
           iter = floaters_transport_stats_map_.begin();
         iter != floaters_transport_stats_map_.end(); iter++)
      {
        for (unsigned int j = 0; j < iter->second.size(); j++)
          {
            iter->second[j]->reset();
          } // for j
      } // for iter
    particle_type_zr_distribution_map_.clear();
  }

  for (ChainStatisticsOSsMap::iterator iter = chains_stats_map_.begin();
       iter != chains_stats_map_.end(); iter++)
    {
      for (unsigned int j = 0; j < iter->second.size(); j++) {
        iter->second[j]->reset();
      } // for j
    } // for iter

  //double time_ns = get_sd()->get_bd()->get_current_time() / FS_IN_NS;
  for (BipartitePairsStatisticsOSMap::iterator
         iter = interaction_stats_map_.begin();
       iter != interaction_stats_map_.end(); iter++)
    {
      iter->second->reset( );
    }
}

void Statistics::set_interrupted(bool tf) {
  ::npctransport_proto::Output output;
  std::ifstream inf(output_file_name_.c_str(), std::ios::binary);
  output.ParseFromIstream(&inf);
  inf.close();
  ::npctransport_proto::Statistics* stats = output.mutable_statistics();
  stats->set_interrupted(tf ? 1 : 0);
  std::ofstream outf(output_file_name_.c_str(), std::ios::binary);
  stats->SerializeToOstream(&outf);
}



/************ private utility methods ****************/

// see doc in .h file
boost::tuple<double, double, double, double>
Statistics::get_interactions_and_interacting
( const ParticlesTemp &floaters, const atom::Hierarchies &chain_roots) const
{
  double site_site_pairs = 0,   interacting_floaters = 0,
    bead_floater_pairs = 0,     chain_floater_pairs = 0;
  for (unsigned int i = 0; i < floaters.size(); ++i) {
    bool floater_found = false;
    for (unsigned int j = 0; j < chain_roots.size(); ++j) {
      bool chain_found = false;
      Pointer<FGChain> cur_chain= get_fg_chain(chain_roots[j]);
      Particles const& chain_particles = cur_chain->get_beads();
      for (unsigned int k = 0; k < chain_particles.size(); ++k) {
        int num = get_sd()->get_scoring()
          ->get_number_of_site_site_interactions
          (floaters[i], chain_particles[k] );
        if (num > 0) {
          IMP_LOG(VERBOSE, "Found " << num
            << " site-site interactions between floater " << i <<
            " and bead " << j << "." << k << std::endl);
          site_site_pairs += num;
          ++bead_floater_pairs;
          if (!floater_found) ++interacting_floaters;
          floater_found = true;
          if (!chain_found) ++chain_floater_pairs;
          chain_found = true;
        } // num > 0
      } // k (particle in chain)
    }// j ( chains)
  }// i (floaters)
  // std::cout << "Return " << interactions << "," << interacting << ","
  //           << bead_partners << "," << chain_partners << std::endl;
  return boost::make_tuple(site_site_pairs,
                           interacting_floaters,
                           bead_floater_pairs,
                           chain_floater_pairs);
}

double Statistics::get_z_distribution_top() const{
  if(get_sd()->get_has_slab())
    {
      return  get_sd()->get_slab_thickness() / 2.0;
    }
  if(get_sd()->get_has_bounding_box())
    {
      return get_sd()->get_box_size() / 4.0;
    }
  return 1.0;
}

double Statistics::get_r_distribution_max() const{
  if(get_sd()->get_has_slab())
    {
      return  get_sd()->get_tunnel_radius();
    }
  if(get_sd()->get_has_bounding_box())
    {
      return get_sd()->get_box_size() / 4.0; // diameter = 50% box edge
    }
  return 1.0;
}

// distribution of particles along z axis
boost::tuple<int, int, int, int>
Statistics::get_z_distribution(const ParticlesTemp& ps) const{
  int zh[4] = {0, 0, 0, 0};
  double top = get_z_distribution_top();
  for (unsigned int i = 0; i < ps.size(); i++) {
    core::XYZ d(ps[i]);
    double z = d.get_z();
    if(z>top) {
      zh[0]++;
    } else if (z>0) {
      zh[1]++;
    } else if (z>-top) {
      zh[2]++;
    } else {
      zh[3]++;
    }
  }
  return boost::make_tuple(zh[0],zh[1],zh[2],zh[3]);
}

void Statistics::fill_in_zr_hist(unsigned int zr_hist[4][3],
                                 ParticlesTemp ps) const{
  double top = get_z_distribution_top();
  double r_max = get_r_distribution_max();
  unsigned int zz, rr;
  for(unsigned int i = 0; i < ps.size() ; i++){
    core::XYZ d(ps[i]);
    double z = d.get_z();
    double r = std::sqrt( std::pow(d.get_x(),2) + std::pow(d.get_y(),2) );
    if(z>top) {
      zz=0;
    } else if (z>0) {
      zz=1;
    } else if (z>-top) {
      zz=2;
    } else {
      zz=3;
    }
    rr = (unsigned int)(floor(2*r/r_max));
    if(rr>2) {
      rr=2;
    }
    zr_hist[zz][rr]++;
  } // i (particles)
}




/************************************************************/
/************* various simple getters and setters *******************/
/************************************************************/



/** returns the model associated with the owned SimulationData */
Model* Statistics::get_model() {
  return get_sd()->get_model();
}


/** returns the model associated with the owned SimulationData */
Model * Statistics::get_model() const {
  return get_sd()->get_model();
}


#undef GET_ASSIGNMENT
#undef GET_VALUE

IMPNPCTRANSPORT_END_NAMESPACE
