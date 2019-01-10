/**
 *  \file BipartitePairsStatisticsOptimizerState.cpp
 *  \brief description.
 *
 *  Copyright 2007-2019 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/BipartitePairsStatisticsOptimizerState.h>
#include <IMP/npctransport/Statistics.h>
#include <IMP/npctransport/SimulationData.h>
#include <IMP/npctransport/Scoring.h>
#include <IMP/npctransport/enums.h>
#include <IMP/npctransport/util.h>
#include <IMP/container/ClosePairContainer.h>
#include <IMP/atom/Simulator.h>
#include <IMP/pair_macros.h>
#include <boost/unordered_set.hpp>

#include <boost/tuple/tuple.hpp>
#include <algorithm>
#include <iterator>
#include <numeric>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

namespace {

  // gets weighted average of d1 and d2 with weights w1 and w2, resp.
  inline double get_weighted_average
  (double d1, double d2, double w1, double w2)
  {
    IMP_USAGE_CHECK(w1>=0.0 && w2 >= 0.0 && w1+w2 > 0.0,
                    "Invalid weights " << w1 << "," << w2);
    return (d1 * w1 + d2 * w2 ) / (w1 + w2);
  }

  // update avg with value new_val, with weights old_w for the
  // old average and delta_w for new_val
  inline void update_weighted_average
  (double &avg, double new_val, double old_w, double delta_w)
  {
    avg = get_weighted_average(avg, new_val, old_w, delta_w);
  }

}



BipartitePairsStatisticsOptimizerState::BipartitePairsStatisticsOptimizerState
(  IMP::npctransport::Statistics* statistics_manager,
   InteractionType interaction_type,
   const ParticlesTemp& particlesI, const ParticlesTemp& particlesII,
   double contact_range, double slack
   )
  : P(statistics_manager ? statistics_manager->get_model() : nullptr,
      "BipartitePairsStatisticsOptimizerState%1%"),
    statistics_manager_(statistics_manager),
    is_reset_(true),
    interaction_type_(interaction_type),
    n_particles_I_(particlesI.size()),
    n_particles_II_(particlesII.size())
{
  close_bipartite_pair_container_ =
    new IMP::container::CloseBipartitePairContainer
    (particlesI, particlesII, contact_range, slack);
  range_ = contact_range;

  n_possible_contacts_ =
    get_maximal_number_of_unordered_pairs(particlesI, particlesII);
  n_sites_I_= n_particles_I_ *
    statistics_manager_->get_sd()->get_sites(interaction_type_.first).size();
  n_sites_II_= n_particles_II_ *
    statistics_manager_->get_sd()->get_sites(interaction_type_.second).size();

  // TODO: do we want to add consecutive pair filter for fg chains?
  reset(); // make sure all counters are 0
}

void BipartitePairsStatisticsOptimizerState::reset() {
  is_reset_ = true; // indicate, for next round of update
  if(statistics_manager_.get() != nullptr) {
    if(statistics_manager_->get_is_activated()) {
      update_always(); // make sure all state variables are updated properly for next round
    }
  } else {
    time_ns_ = 0.0;
  }
  P::reset();
  n_updates_ = 0;
  stats_time_ns_ = 0;
  off_stats_time_ns_ = 0.0;
  off_I_stats_time_ns_ = 0.0;
  off_II_stats_time_ns_ = 0.0;
  on_stats_time_ns_ = 0.0;
  on_I_stats_time_ns_ = 0.0;
  on_II_stats_time_ns_ = 0.0;
}

namespace {
  typedef std::map<ParticleIndex, std::vector<unsigned int> >
    t_bound_sites_by_pi_map;

  //! updates a map of bound sites by map with new counts
  //! of bound sites for the specified particle index pi
  void accumulate_bound_sites_by_pi
    (t_bound_sites_by_pi_map& bound_sites_by_pi,
     ParticleIndex pi,
     std::vector<unsigned int> new_bound_sites)
  {
    if(bound_sites_by_pi.find(pi) == bound_sites_by_pi.end()) {
      bound_sites_by_pi[pi]= new_bound_sites;
    }
    else {
      std::vector<unsigned int>& bound_sites= bound_sites_by_pi[pi];
      IMP_USAGE_CHECK(bound_sites.size()==new_bound_sites.size(),
                      "Expected bound_sites and new_bound_sites to have identical sizes");
      std::transform(bound_sites.begin(),
                     bound_sites.end(),
                     new_bound_sites.begin(),
                     bound_sites.begin(),
                     std::plus<double>());
    }
  }

  //! returns total number of bound sites (>=1 interactions)
  //! for all particles in bound_sites_by_pi
  unsigned int get_number_of_bound_sites
  (t_bound_sites_by_pi_map const& bound_sites_by_pi){
    unsigned int ret_value= 0;
    for(t_bound_sites_by_pi_map::const_iterator iter= bound_sites_by_pi.begin();
        iter!= bound_sites_by_pi.end();
        iter++){
      for(unsigned int i=0; i<iter->second.size(); i++){
        ret_value+= (iter->second[i] > 0);
      }
    }
    return ret_value;
  }
} //namespace {


// count all the pairs that are currently in contact
// and update stats
void BipartitePairsStatisticsOptimizerState::do_update(unsigned int call_num)
{
  // Get simulation time and reset if needed:
  atom::Simulator* simulator =
    dynamic_cast< atom::Simulator* >( get_optimizer() );
  IMP_USAGE_CHECK( simulator, "Optimizer must be a simulator in order to use "
                   "BipartitePairsStatisticsOptimizerState, for time stats" );
  double new_time_ns = simulator->get_current_time() / FS_IN_NS;
  if(is_reset_ == true){
    IMP_LOG(PROGRESS, "Starting bpsos with simulation time " << new_time_ns);
    time_ns_ = new_time_ns;
    is_reset_ = false;
  }
  double elapsed_time_ns = new_time_ns - time_ns_;
  IMP_LOG(PROGRESS,
          "Bipartite stats new-time/old-time/elapsed-time "
          << std::setprecision(3) << new_time_ns << " "
          << std::setprecision(3) << time_ns_ << " "
          << std::setprecision(3) << elapsed_time_ns
          << std::endl);
  IMP_LOG(PROGRESS, "Bipartite stats time span before current update:"
          << " misc " << std::setprecision(3) << stats_time_ns_
          << " / off " << std::setprecision(3) << off_stats_time_ns_
          << " / on " << std::setprecision(3) << on_stats_time_ns_
          << " / onI " << std::setprecision(3) << on_I_stats_time_ns_
          << " / onII " << std::setprecision(3) << on_II_stats_time_ns_
          << std::endl);

  // Update the lists of bound particles and their interactions,
  // for all bipartite pairs of distinct particles
  close_bipartite_pair_container_->do_score_state_before_evaluate(); // refresh
  t_particle_index_ordered_set new_bounds_I, new_bounds_II;
  t_particle_index_pair_ordered_set new_contacts; // more efficient if ordered set
  std::map<ParticleIndex, std::vector<unsigned int> >
    bound_sites_I_by_pi;
  std::map<ParticleIndex, std::vector<unsigned int> >
    bound_sites_II_by_pi;
  unsigned int n_nonspecific_interactions(0);
  IMP_CONTAINER_FOREACH(IMP::container::CloseBipartitePairContainer,
                        close_bipartite_pair_container_,
                        {
                          ParticleIndexPair const& pip = _1;
                          unsigned int n_site_site_contacts;
                          std::vector<unsigned int> bound_sites_I;
                          std::vector<unsigned int> bound_sites_II;
                          bool is_nonspecific_interaction;
                          boost::tie(n_site_site_contacts,
                                     bound_sites_I,
                                     bound_sites_II,
                                     is_nonspecific_interaction)=
                            statistics_manager_->get_sd()->get_scoring()
                            ->get_site_interactions_statistics(pip[0], pip[1]);
                          if(n_site_site_contacts>0){
                            new_bounds_I.insert(pip[0]);
                            new_bounds_II.insert(pip[1]);
                            new_contacts.insert
                              ( make_unordered_particle_index_pair( _1 ) );
                            if(call_num % 10 == 0) {
                              //                              statistics_manager->get_fgs_markov_states->update_contact(pip[0],
                              //                                                                                        pip[1],
                              //                                                                                        new_time_ns)
                            }
                          }
                          accumulate_bound_sites_by_pi(bound_sites_I_by_pi,
                                                          pip[0],
                                                          bound_sites_I);
                          accumulate_bound_sites_by_pi(bound_sites_II_by_pi,
                                                          pip[1],
                                                          bound_sites_II);
                          if(is_nonspecific_interaction){
                            n_nonspecific_interactions++;
                          }
                        });
  unsigned int n_bound_sites_I= get_number_of_bound_sites(bound_sites_I_by_pi);
  unsigned int n_bound_sites_II= get_number_of_bound_sites(bound_sites_II_by_pi);
  double fraction_bound_sites_I= (n_sites_I_>0) ? n_bound_sites_I/ (n_sites_I_+.0) : 0.0;
  double fraction_bound_sites_II= (n_sites_II_>0) ? n_bound_sites_II/ (n_sites_II_+.0) : 0.0;
  double fraction_nonspecific_I= n_nonspecific_interactions / (n_particles_I_+0.);
  double fraction_nonspecific_II= n_nonspecific_interactions / (n_particles_II_+0.);

  IMP_LOG(PROGRESS,
          new_bounds_I.size() << "/" << n_particles_I_
          << " bound-I(" << interaction_type_.first<< "); "
          << new_bounds_II.size() << "/" << n_particles_II_
          << " bound-II" << interaction_type_.second << "); "
          << new_contacts.size() << " contacts" << std::endl);

  // Update avg_ncontacts_ and fraction_bound_sites_I/II:
  if(elapsed_time_ns>0)
    {
      update_weighted_average(avg_ncontacts_, // old
                              new_contacts.size(), // new
                              stats_time_ns_,
                              elapsed_time_ns);
      update_weighted_average(avg_fraction_bound_sites_I_, // old
                              fraction_bound_sites_I, // new
                              stats_time_ns_,
                              elapsed_time_ns);
      update_weighted_average(avg_fraction_bound_sites_II_, // old
                              fraction_bound_sites_II, // new
                              stats_time_ns_,
                              elapsed_time_ns);
      update_weighted_average(avg_fraction_nonspecific_I_,
                              fraction_nonspecific_I,
                              stats_time_ns_,
                              elapsed_time_ns);
      update_weighted_average(avg_fraction_nonspecific_II_,
                              fraction_nonspecific_II,
                              stats_time_ns_,
                              elapsed_time_ns);
    }

  // Update on-off rates:
  if(elapsed_time_ns > 0)
    {
      // All doubles for division:
      double n_contacts_lost, n_contacts_gained; // double for divide
      double n_bounds_I_lost, n_bounds_I_gained,
        n_bounds_II_lost, n_bounds_II_gained; // double for divide
      boost::tie(n_contacts_lost, n_contacts_gained) =
        get_n_lost_and_gained( contacts_, new_contacts);
      boost::tie(n_bounds_I_lost, n_bounds_I_gained) =
        get_n_lost_and_gained( bounds_I_, new_bounds_I);
      boost::tie(n_bounds_II_lost, n_bounds_II_gained) =
        get_n_lost_and_gained( bounds_II_, new_bounds_II);
      double n_contacts_before = contacts_.size();
      double n_bounds_I_before = bounds_I_.size();
      double n_bounds_II_before = bounds_II_.size();
      double n_unbounds_I_before = n_particles_I_ - n_bounds_I_before;
      double n_unbounds_II_before = n_particles_II_ - n_bounds_II_before;
      IMP_LOG(PROGRESS,
                 "lostI " << n_bounds_I_lost
              << " gainedI " << n_bounds_I_gained
              << " lostII " << n_bounds_II_lost
              << " gainedII " << n_bounds_II_gained
              << std::endl);
      // OFF
      if(  n_contacts_before > 0 )
        {
          IMP_USAGE_CHECK( n_bounds_I_before > 0 && n_bounds_II_before > 0,
                           "positive contacts but no bounds type I or II");
          {
            double weighted_time_ns = ( n_contacts_before * elapsed_time_ns);
            double off_per_contact_per_ns =
              n_contacts_lost / weighted_time_ns;
            update_weighted_average(avg_off_per_contact_per_ns_,
                                    off_per_contact_per_ns,
                                    off_stats_time_ns_,
                                    weighted_time_ns);
            off_stats_time_ns_ += weighted_time_ns;
          }
          {
            double weighted_time_ns=( n_bounds_I_before * elapsed_time_ns);
            double off_per_bound_I_per_ns =
              n_bounds_I_lost / weighted_time_ns;
            update_weighted_average(avg_off_per_bound_I_per_ns_,
                                    off_per_bound_I_per_ns,
                                    off_I_stats_time_ns_,
                                    weighted_time_ns);
            off_I_stats_time_ns_ += weighted_time_ns;
          }
          {
            double weighted_time_ns=( n_bounds_II_before * elapsed_time_ns);
            double off_per_bound_II_per_ns =
              n_bounds_II_lost / weighted_time_ns;
            update_weighted_average(avg_off_per_bound_II_per_ns_,
                                    off_per_bound_II_per_ns,
                                    off_II_stats_time_ns_,
                                    weighted_time_ns);
            off_II_stats_time_ns_ += weighted_time_ns;
          }
        }
      if( n_unbounds_I_before > 0) {
        double weighted_time_ns =  n_unbounds_I_before * elapsed_time_ns;
        double on_per_unbound_I_per_ns =
          n_bounds_I_gained / weighted_time_ns;
        update_weighted_average(avg_on_per_unbound_I_per_ns_,
                                on_per_unbound_I_per_ns,
                                on_I_stats_time_ns_,
                                weighted_time_ns);
        on_I_stats_time_ns_ += weighted_time_ns;
      }
      if( n_unbounds_II_before > 0) {
        double weighted_time_ns =  n_unbounds_II_before * elapsed_time_ns;
        double on_per_unbound_II_per_ns =
          n_bounds_II_gained / weighted_time_ns;
        update_weighted_average(avg_on_per_unbound_II_per_ns_,
                                on_per_unbound_II_per_ns,
                                on_II_stats_time_ns_,
                                weighted_time_ns);
        on_II_stats_time_ns_ += weighted_time_ns;
      }
      int n_missing_contacts_before= // simplifying assumption: each pair of particles can form at most a single site-site contact - if not true, this measure is skewed
        (n_sites_I_ - n_contacts_before)
        * (n_sites_II_ - n_contacts_before);
      if( n_missing_contacts_before > 0 )
        {
          double weighted_time_ns =
            n_missing_contacts_before * elapsed_time_ns;
          double on_per_missing_contact_per_ns =
            n_contacts_gained / weighted_time_ns;
          update_weighted_average(avg_on_per_missing_contact_per_ns_,
                                  on_per_missing_contact_per_ns,
                                  on_stats_time_ns_,
                                  weighted_time_ns);
          on_stats_time_ns_ += weighted_time_ns;
      }

      IMP_IF_LOG(PROGRESS) {
        if( n_contacts_before > 0 || n_unbounds_I_before > 0 )
          {
            IMP_LOG(PROGRESS, "prev_contacts " << contacts_.size()
                    << " new_contacts " << new_contacts.size()
                    << " contacts_lost " << n_contacts_lost
                    << " contacts_gained " << n_contacts_gained
                    << std::endl);
            IMP_LOG(PROGRESS, "Stats time:"
                    << " misc " << std::setprecision(3) << stats_time_ns_
                  << " / off " << std::setprecision(3) << off_stats_time_ns_
                  << " / offI " << std::setprecision(3) << off_I_stats_time_ns_
                  << " / offII " << std::setprecision(3) << off_II_stats_time_ns_
                    << " / on " << std::setprecision(3) << on_stats_time_ns_
                    << " / onI " << std::setprecision(3) << on_I_stats_time_ns_
                    << " / onII " << std::setprecision(3) << on_II_stats_time_ns_
                    << std::endl);
            IMP_LOG(PROGRESS, "Average OFF per ns: "
                    << " per_contact " << avg_off_per_contact_per_ns_
                    << " per bound I " << avg_off_per_bound_I_per_ns_
                    << " per bound II " << avg_off_per_bound_II_per_ns_
                    << std::endl);
            IMP_LOG(PROGRESS, "Average ON per ns: "
                    << " per_missing contact "
                    << avg_on_per_missing_contact_per_ns_
                    << " per unbound I " << avg_on_per_unbound_I_per_ns_
                    << " per unbound II " << avg_on_per_unbound_II_per_ns_
                    << std::endl);
          } // if
      } // IMP_IF_LOG
    } // ON/OFF rate update

  if(elapsed_time_ns>0)
    {
      // TODO: next lines - n_particles_XX_ not dynamic
      // TODO: pct is misleading - it is fraction
      double pct_bound_particles_I =
        (new_bounds_I.size() + 0.0) / n_particles_I_;
      update_weighted_average( avg_pct_bound_particles_I_,
                               pct_bound_particles_I,
                               stats_time_ns_,
                               elapsed_time_ns);
      double pct_bound_particles_II =
        (new_bounds_II.size() + 0.0) / n_particles_II_;
      update_weighted_average( avg_pct_bound_particles_II_,
                               pct_bound_particles_II,
                               stats_time_ns_,
                               elapsed_time_ns);
    }

  // update records
  n_updates_++;
  bounds_I_ = new_bounds_I;
  bounds_II_ = new_bounds_II;
  contacts_ = new_contacts;
  stats_time_ns_ += elapsed_time_ns;
  time_ns_ = new_time_ns;
}


IMPNPCTRANSPORT_END_NAMESPACE
