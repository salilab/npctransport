/**
 *  \file Statistics.cpp
 *  \brief statistics and order parameters about the simulations
 *         that is associated with a SimulationData object
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/Statistics.h>
#include <IMP/npctransport/SimulationData.h>
#include <IMP/npctransport/particle_types.h>
#include <IMP/npctransport/protobuf.h>
#include <IMP/npctransport/util.h>
#include <IMP/npctransport/io.h>
#include <IMP/npctransport/typedefs.h>
#include <IMP/npctransport/util.h>
#include <IMP/algebra/vector_generators.h>
#include <IMP/atom/estimates.h>
#include <IMP/atom/distance.h>
#include <IMP/atom/Diffusion.h>
#include <IMP/atom/Selection.h>
#include <IMP/base/log.h>
#include <IMP/core/pair_predicates.h>
#include <IMP/core/XYZR.h>
#include <IMP/core/generic.h>
#include <IMP/display/LogOptimizerState.h>
#include <numeric>
#include <set>

#ifdef IMP_NPC_GOOGLE
IMP_GCC_PUSH_POP(diagnostic push)
IMP_GCC_PRAGMA(diagnostic ignored "-Wsign-compare")
#include "third_party/npc/npctransport/data/npctransport.pb.h"
IMP_GCC_PUSH_POP(diagnostic pop)
#else
#include <IMP/npctransport/internal/npctransport.pb.h>
#endif


// TODO: turn into a template inline in unamed space?
/**
   updates (message).field() with a weighted average of its current
   value and new_value, giving weight n_old_frames, n_new_frames to each,
   respectively.
*/
#define UPDATE_AVG(n_frames, n_new_frames, message, field, new_value)     \
  (message).set_##field(static_cast<double>(n_frames *(message).field() + \
                                            n_new_frames *new_value) /    \
                        (n_frames + n_new_frames));


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
( ParticlesTemp chain_beads )
{
  if(chain_beads.size() == 0)
    return;
  core::ParticleType type = core::Typed(chain_beads[0]).get_type();
  // chain stats
  IMP_NEW( ChainStatisticsOptimizerState, csos,
           (chain_beads, statistics_interval_frames_ ) );
  chains_stats_map_[type].push_back( csos );
  // stats for each chain particle
  fgs_bodies_stats_map_[type].push_back( BodyStatisticsOptimizerStates() );
  for (unsigned int k = 0; k < chain_beads.size(); ++k) {
    Particle* p = chain_beads[k];
    IMP_ALWAYS_CHECK( core::Typed(p).get_type() == type,
                      "All beads in chain must be of same type for now",
                      base::ValueException);
    IMP_NEW(BodyStatisticsOptimizerState, bsos,
            ( p, statistics_interval_frames_ ) );
    fgs_bodies_stats_map_[type].back().push_back( bsos );
  }  // for k
}

//! add statistics about a floater particle
void Statistics::add_floater_stats
( Particle* p )
{
  core::ParticleType type = core::Typed(p).get_type();
  IMP_NEW(BodyStatisticsOptimizerState, bsos,
          (p, statistics_interval_frames_));
  floaters_stats_map_[type].push_back(bsos);
  if (get_sd()->get_has_slab())
    {  // only if has tunnel
      IMP_NEW(ParticleTransportStatisticsOptimizerState, ptsos,
              (p,
               -0.5 * get_sd()->get_slab_thickness(),  // tunnel bottom
               0.5 * get_sd()->get_slab_thickness() )   // tunnel top
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
    IMP_CONTAINER_FOREACH // _1 is the particle index
      ( SingletonContainer,
        get_sd()->get_diffusers(),
        {
          if (core::Typed(get_model(), _1).get_type() == type0)
            { set0.push_back( get_model()->get_particle(_1) ); }
          if (core::Typed(get_model(), _1).get_type() == type1)
            { set1.push_back( get_model()->get_particle(_1) ); }
        }
        );
    // TODO: make params
    double stats_contact_range = 2;// (dist + 2 * slack) for ClosePairContainers
    double stats_slack=1;
    double stats_dist = stats_contact_range - stats_slack * 2;
  IMP_LOG(PROGRESS,
          "Interaction " << type0.get_string() << ", " << type1.get_string()
          << "  sizes: " << set0.size() << ", " << set1.size()
          << " statistics range: " << stats_contact_range
          << std::endl);
  if (set0.size() > 0 && set1.size() > 0)
    {
      IMP_NEW(BipartitePairsStatisticsOptimizerState, bpsos,
              (get_model(), interaction_type, set0, set1,
               stats_dist, stats_slack));
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
                    base::ValueException);
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
  if (get_sd()->get_has_slab())
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


// @param nf_new number of new frames accounted for in this statistics update
void Statistics::update
( const boost::timer &timer,
  unsigned int nf_new)
{
  IMP_OBJECT_LOG;
  ::npctransport_proto::Output output;

  std::ifstream inf(output_file_name_.c_str(), std::ios::binary);
  output.ParseFromIstream(&inf);
  inf.close();
  ::npctransport_proto::Statistics &stats = *output.mutable_statistics();
  int nf = stats.number_of_frames();
  if (is_stats_reset_) {  // stats was just reset
    // TODO: what's if multiple trials?
    nf = 0;
    is_stats_reset_ = false;
    for (int i = 0; i < stats.floaters().size(); i++) {
      (*stats.mutable_floaters(i)).clear_transport_time_points_ns();
    }
  }
  std::cout << "Updating statistics file " << output_file_name_
            << " that currently has " << nf << " frames, with " << nf_new
            << " additional frames" << std::endl;

  // gather the statistics one by one
  ParticlesTemp all;
  base::Vector<ParticlesTemps> fgs;

  // General FG statistics:
  {
    ParticleTypeSet const &fgt = get_sd()->get_fg_types();
    for ( ParticleTypeSet::const_iterator
            it = fgt.begin(); it != fgt.end(); it++ )
      {
        ParticlesTemp ps = get_sd()->get_particles_of_type( *it );
        if(ps.size() == 0) continue; // TODO: add warning?
                                     // can happen only if dynamic changes...
        fgs.push_back(ParticlesTemps());
        unsigned int i = find_or_add_fg_of_type( &stats, *it );
        for (unsigned int j = 0; j < ps.size(); ++j)
          {
            atom::Hierarchy hj(ps[j]);
            ParticlesTemp chainj = get_as<ParticlesTemp>(atom::get_leaves(hj));
            all += chainj;
            fgs.back().push_back(chainj);
#ifdef IMP_NPCTRANSPORT_USE_IMP_CGAL
            double volume = atom::get_volume(hj);
            double radius_of_gyration = atom::get_radius_of_gyration(chainj);
#else
            double volume = -1.;
            double radius_of_gyration = -1.;
#endif
            UPDATE_AVG(nf, nf_new, *stats.mutable_fgs(i), volume, volume);
            double length =
              core::get_distance(core::XYZ(chainj[0]), core::XYZ(chainj.back()));
            UPDATE_AVG(nf, nf_new, *stats.mutable_fgs(i), length, length);
            UPDATE_AVG(nf, nf_new, *stats.mutable_fgs(i), radius_of_gyration,
                       radius_of_gyration);
          } // for j
      } // for it
  }

  // FG body stats
  for (FGsBodyStatisticsOSsMap::iterator
         it = fgs_bodies_stats_map_.begin();
       it != fgs_bodies_stats_map_.end(); it++)
    {
      unsigned int i = find_or_add_fg_of_type( &stats, it->first );
      FGsBodyStatisticsOSs& fbs_i = it->second;
      for (unsigned int j = 0; j < fbs_i.size(); ++j)
        {
          BodyStatisticsOptimizerStates& fbs_ij = fbs_i[j];
        for (unsigned int k = 0; k < fbs_ij.size(); ++k)
          {
            BodyStatisticsOptimizerState* fbs_ijk = fbs_ij[k];
            unsigned int per_frame = fbs_i.size() * fbs_ij.size();
            unsigned int cnf = (nf) * per_frame + j * fbs_ij.size() + k;
            UPDATE_AVG(cnf, nf_new, *stats.mutable_fgs(i),
                       particle_correlation_time,
                       fbs_ijk->get_correlation_time());
            UPDATE_AVG(cnf, nf_new, *stats.mutable_fgs(i),
                       particle_diffusion_coefficient,
                       fbs_ijk->get_diffusion_coefficient());
            fbs_ijk->reset();
          } // for k
        } // for j
    } // for it

  // FG entire chain stats
  for (ChainStatisticsOSsMap::iterator
         it = chains_stats_map_.begin(); it != chains_stats_map_.end(); it++)
    {
      unsigned int i = find_or_add_fg_of_type( &stats, it->first );
      ChainStatisticsOptimizerStates& cs_i = it->second;
      for (unsigned int j = 0; j < cs_i.size(); ++j)
        {
          unsigned int cnf = nf * cs_i.size() + j;
          UPDATE_AVG(cnf, nf_new,
                     *stats.mutable_fgs(i), chain_correlation_time,
                     cs_i[j]->get_correlation_time());
          UPDATE_AVG(cnf, nf_new, *stats.mutable_fgs(i),
                     chain_diffusion_coefficient,
                     cs_i[j]->get_diffusion_coefficient());
          Floats df = cs_i[j]->get_diffusion_coefficients();
          UPDATE_AVG(cnf, nf_new, *stats.mutable_fgs(i),
                     local_diffusion_coefficient,
                     std::accumulate(df.begin(), df.end(), 0.0) / df.size());
          cs_i[j]->reset();
        } // for j
    } // for it


  // Floaters general body stats
  for (BodyStatisticsOSsMap::iterator
         it = floaters_stats_map_.begin() ;
       it != floaters_stats_map_.end(); it++)
    {
      unsigned int i = find_or_add_floater_of_type( &stats, it->first );
      BodyStatisticsOptimizerStates& bsos = it->second;
      for (unsigned int j = 0; j < bsos.size(); j++)
        {
          int cnf = (nf) * bsos.size() + j; // each old frame is for size particles
          UPDATE_AVG(cnf, nf_new,  // TODO: is nf_new correct? I think so
                     *stats.mutable_floaters(i), diffusion_coefficient,
                     bsos[j]->get_diffusion_coefficient());
          UPDATE_AVG(cnf, nf_new,  // TODO: is nf_new correct? I think so
                     *stats.mutable_floaters(i), correlation_time,
                     bsos[j]->get_correlation_time());
          bsos[j]->reset();
        } // for j
    } // for it

  // Floaters avg number of transports per particle
  if (get_sd()->get_has_slab()){
    for (ParticleTransportStatisticsOSsMap::iterator
           it1 = floaters_transport_stats_map_.begin();
         it1 != floaters_transport_stats_map_.end() ; it1++)
      {
        unsigned int i = find_or_add_floater_of_type( &stats, it1->first );
        ParticleTransportStatisticsOptimizerStates& pts_i = it1->second;
        // fetch old ones from stats msg, add new ones and rewrite all:
        std::set<double> times_i
          ( stats.floaters(i).transport_time_points_ns().begin(),
            stats.floaters(i).transport_time_points_ns().end() );
        for (unsigned int j = 0; j < pts_i.size(); ++j)
          {
            Floats const &new_times_ij =
              pts_i[j]->get_transport_time_points_in_ns();
            for (unsigned int k = 0; k < new_times_ij.size(); k++)
              {
                times_i.insert(new_times_ij[k]);
              } // for k
          } // for j
        (*stats.mutable_floaters(i)).clear_transport_time_points_ns();
        for (std::set<double>::const_iterator it2 = times_i.begin();
             it2 != times_i.end(); it2++)
          {
            (*stats.mutable_floaters(i)).add_transport_time_points_ns(*it2);
          } // for it2
        // update avg too:
        double avg_n_transports_i = times_i.size() * 1.0 / pts_i.size();
        (*stats.mutable_floaters(i)).set_avg_n_transports( avg_n_transports_i );
      } // for it1
  }

  // FG-floaters interactions
  ParticlesTemps floaters_list;
  ParticleTypeSet const &ft = get_sd()->get_floater_types();
  for ( ParticleTypeSet::const_iterator
          it = ft.begin(); it != ft.end(); it++ )
      {
        ParticlesTemp ps = get_sd()->get_particles_of_type( *it );
        if(ps.size() == 0) // TODO: makes sense that this should happen?
          continue;
        floaters_list.push_back(ps);
        all += ps;
        unsigned int i = find_or_add_floater_of_type( &stats, *it );
        double interactions, interacting, bead_partners, chain_partners;
        boost::tie(interactions, interacting, bead_partners, chain_partners) =
          get_interactions_and_interacting(ps, fgs);
        UPDATE_AVG(nf, nf_new, *stats.mutable_floaters(i), interactions,
                   interactions / ps.size());
        UPDATE_AVG(nf, nf_new, *stats.mutable_floaters(i), interacting,
                 interacting / ps.size());
        UPDATE_AVG(nf, nf_new, *stats.mutable_floaters(i),
                 interaction_partner_chains, chain_partners / interacting);
        UPDATE_AVG(nf, nf_new, *stats.mutable_floaters(i),
                   interaction_partner_beads, bead_partners / interacting);
      } // for it

  // update statistics gathered on interaction rates
  for (BipartitePairsStatisticsOSMap::iterator
         it = interaction_stats_map_.begin();
       it != interaction_stats_map_.end(); it++)
    {
      IMP_LOG(PROGRESS, "adding interaction statistics "
              << it->first << std::endl);
      unsigned int i = find_or_add_interaction_of_type( &stats, it->first );
      BipartitePairsStatisticsOptimizerState* bps_i = it->second;
      ::npctransport_proto::Statistics_InteractionStats *pOutStats_i =
          stats.mutable_interactions(i);
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
                base::ValueException);
    }
    // save the rest of the interactions info
    Int n0 = bps_i->get_number_of_particles_1();
    Int n1 = bps_i->get_number_of_particles_2();
    Float avg_contacts_num = bps_i->get_average_number_of_contacts();

    UPDATE_AVG(nf, nf_new, *pOutStats_i, avg_contacts_per_particle0,
               avg_contacts_num / n0);
    UPDATE_AVG(nf, nf_new, *pOutStats_i, avg_contacts_per_particle1,
               avg_contacts_num / n1);
    UPDATE_AVG(nf, nf_new, *pOutStats_i, avg_pct_bound_particles0,
               bps_i->get_average_percentage_bound_particles_1());
    UPDATE_AVG(nf, nf_new, *pOutStats_i, avg_pct_bound_particles1,
               bps_i->get_average_percentage_bound_particles_2());
    bps_i->reset();
  }

  // MISC:
  {
    // Todo: define better what we want of timer
    stats.set_seconds_per_iteration(timer.elapsed());
    stats.set_number_of_frames(nf + nf_new);
    const double FS_IN_NS = 1.0E+6;
    double sim_time_ns = const_cast<SimulationData *>( get_sd() )
      ->get_bd()->get_current_time() / FS_IN_NS;
    stats.set_bd_simulation_time_ns( sim_time_ns );
    double total_energy  =
      get_sd()->get_bd()->get_scoring_function()->evaluate(false);
    UPDATE_AVG(nf, nf_new, stats, energy_per_particle,  // TODO: reset?
               total_energy / all.size());
    ::npctransport_proto::Statistics_EnergyByTime*
        sebt = stats.add_energies_by_time();
    sebt->set_time_ns(sim_time_ns);
    sebt->set_energy(total_energy);
  }

  // TODO: disable this for now
  //  ::npctransport_proto::Conformation *conformation =
  //   output.mutable_conformation();
  //    save_pb_conformation(get_diffusers(), sites_, conformation);

  // save RMF for future restarts
  {
    std::string buf;
    RMF::FileHandle fh = RMF::create_rmf_buffer(buf);
    const_cast<SimulationData *>( get_sd() )->link_rmf_file_handle(fh);
    rmf::save_frame(fh, 0);
    output.set_rmf_conformation(buf);
  }

  // dump to file
  std::ofstream outf(output_file_name_.c_str(), std::ios::binary);
  output.SerializeToOstream(&outf);
}

void Statistics::reset_statistics_optimizer_states()
{
  is_stats_reset_ = true;  // indicate to update()
  get_sd()->get_bd()->set_current_time(0.0);

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
        } // for j
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
  }

  for (ChainStatisticsOSsMap::iterator iter = chains_stats_map_.begin();
       iter != chains_stats_map_.end(); iter++)
    {
      for (unsigned int j = 0; j < iter->second.size(); j++) {
        iter->second[j]->reset();
      } // for j
    } // for iter

  for (BipartitePairsStatisticsOSMap::iterator
         iter = interaction_stats_map_.begin();
       iter != interaction_stats_map_.end(); iter++)
    {
      iter->second->reset();
    }
}

void Statistics::set_interrupted(bool tf) {
  ::npctransport_proto::Output output;
  std::ifstream inf(output_file_name_.c_str(), std::ios::binary);
  output.ParseFromIstream(&inf);
  inf.close();
  ::npctransport_proto::Statistics &stats = *output.mutable_statistics();
  stats.set_interrupted(tf ? 1 : 0);
  std::ofstream outf(output_file_name_.c_str(), std::ios::binary);
  stats.SerializeToOstream(&outf);
}



/************ private utility methods ****************/

// number of site-site interactions between a and b
int Statistics::get_number_of_interactions(Particle *a, Particle *b) const {
  using namespace IMP::core;
  double range = get_sd()->get_range();
  if ( get_distance(XYZR(a), XYZR(b)) > range )
    return 0;
  const algebra::Vector3Ds &sa =
    get_sd()->get_sites( Typed(a).get_type() );
  const algebra::Vector3Ds  &sb =
    get_sd()->get_sites( Typed(b).get_type() );
  int ct = 0; // TODO: shouldn't sites be transformed ?!
  for (unsigned int i = 0; i < sa.size(); ++i) {
    for (unsigned int j = 0; j < sb.size(); ++j) {
      if (algebra::get_distance(sa[i], sb[j]) < range) {
        ++ct;
      }
    }
  }
  return ct;
}

// see doc in .h file
boost::tuple<double, double, double, double>
Statistics::get_interactions_and_interacting(
    const ParticlesTemp &kaps, const base::Vector<ParticlesTemps> &fgs) const {
  double interactions = 0, interacting = 0, bead_partners = 0,
         chain_partners = 0;
  for (unsigned int i = 0; i < kaps.size(); ++i) {
    bool kap_found = false;
    for (unsigned int j = 0; j < fgs.size(); ++j) {
      for (unsigned int k = 0; k < fgs[j].size(); ++k) {
        bool chain_found = false;
        for (unsigned int l = 0; l < fgs[j][k].size(); ++l) {
          int num = get_number_of_interactions(kaps[i], fgs[j][k][l]);
          if (num > 0) {
            interactions += num;
            ++bead_partners;
            if (!kap_found) ++interacting;
            kap_found = true;
            if (!chain_found) ++chain_partners;
            chain_found = true;
          }
        }
      }
    }
  }
  return boost::make_tuple(interactions, interacting, bead_partners,
                           chain_partners);
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