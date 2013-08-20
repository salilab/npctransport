/**
 *  \file npctransport/Statistics.h
 *  \brief statistics and order parameters about the simulations
 *         that is associated with a SimulationData object
 *
 *  \author Barak Raveh, Daniel Russell
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_STATISTICS_H
#define IMPNPCTRANSPORT_STATISTICS_H


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
#include "BodyStatisticsOptimizerState.h"
#include "ParticleTransportStatisticsOptimizerState.h"
#include "ChainStatisticsOptimizerState.h"
#include "BipartitePairsStatisticsOptimizerState.h"
#include "Parameter.h"
#include "Scoring.h"
#include "typedefs.h"

#include <boost/timer.hpp>
#include <string>

#ifdef SWIG
namespace boost {
struct timer {};
}
#endif


IMPNPCTRANSPORT_BEGIN_NAMESPACE

class SimulationData; // fwd incomplete declaration

class IMPNPCTRANSPORTEXPORT Statistics : public base::Object {
 private:
  base::UncheckedWeakPointer<SimulationData> owner_sd_;

  // interval of simulation frames for gathering stats
  Parameter<int> statistics_interval_frames_;

  // the file to which simulation statistics are dumped:
  std::string output_file_name_;

  // statistics about all fgs, per particle, per chain, per particle type
  typedef std::vector< BodyStatisticsOptimizerStates >
    FGsBodyStatisticsOSs;
  typedef base::map<core::ParticleType, FGsBodyStatisticsOSs>
    FGsBodyStatisticsOSsMap;
  FGsBodyStatisticsOSsMap fgs_bodies_stats_map_;

  // statistics about all floaters (kaps etc.), per particle type
  typedef base::map<core::ParticleType, BodyStatisticsOptimizerStates>
    BodyStatisticsOSsMap;
  BodyStatisticsOSsMap floaters_stats_map_;

  // transport statistics about all floaters (kaps etc.) per particle type
  typedef base::map< core::ParticleType,
    ParticleTransportStatisticsOptimizerStates>
    ParticleTransportStatisticsOSsMap;
  ParticleTransportStatisticsOSsMap floaters_transport_stats_map_;

  // statistics about entire FG chains, for each FG type
  typedef base::map<core::ParticleType, ChainStatisticsOptimizerStates>
    ChainStatisticsOSsMap;
  ChainStatisticsOSsMap chains_stats_map_;

  // statistics of pairs of interactions for each interaction type
  typedef base::map< npctransport::InteractionType,
    base::PointerMember<BipartitePairsStatisticsOptimizerState> >
    BipartitePairsStatisticsOSMap;
  BipartitePairsStatisticsOSMap interaction_stats_map_;

  // true if statistics were recently reset, so that update()
  // should restart averaging from 0 frames
  mutable bool is_stats_reset_;


 public:
  /**
     @param sd the sd that owns and uses this statistics object
     @param statistics_interval_frames the interval of simulation frames for gathering
                                       statistics
     @param output_file_name name of output file to which to dump statistics (or update
                             if it already exists) when calling update()
   */
  Statistics(SimulationData* sd,
             unsigned int statistics_interval_frames,
             std::string output_file_name);

  //! add statistics about an FG chain
  /** add statistics about an FG chains

      @param chain_beads the beads in the FG chain in consecutive order
  */
  void add_fg_chain_stats(IMP::ParticlesTemp chain_beads);

  //! add statistics about a floater particle
  /** add statistics about a floater particle

      @param p the particle
  */
  void add_floater_stats(IMP::Particle* p);

  //! add statistics about interactions between particles of type 0 and 1
  /** add statistics about interactions between particles of type 0 and 1
      (order does not matter)

      @param type0 type of first intracting particles
      @param type1 type of other interacting particles
  */
  void add_interaction_stats
    ( core::ParticleType type0, core::ParticleType type1);

  //! add all statistics-related optimizer states to o
  /**
      @param o optimizer to which optimizer states are added,
               use get_sd()->get_bd() if nullptr
      @return the list of optimizer states that were addedxb

      @note If called for more than one optimizer, only the last
            optimizer will be guaranteed to work well
  */
  OptimizerStates add_optimizer_states(Optimizer* o = nullptr);

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
  void update(const boost::timer &timer,
                         unsigned int nf_new = 1);

  /** resets all the counters of any statistics counters,
      and the simulation time to zero */
  void reset_statistics_optimizer_states();

  //! Loads the stats file and set the interrupted flag to true
  void set_interrupted(bool tf);

  /************************************************************/
  /************* various simple getters and setters *******************/
  /************************************************************/

  /** returns the model associated with the owned SimulationData */
  Model* get_model();

#ifndef SWIG
  /** returns the model associated with the owned SimulationData */
  Model* get_model() const;
#endif

  /** return the SimulationData object that owns this ScoringFunction */
  SimulationData* get_sd() {
    return owner_sd_;
  }

#ifndef SWIG
  /** return the SimulationData object that owns this ScoringFunction */
  SimulationData const* get_sd() const{
    return owner_sd_;
  }
#endif


 private:

  // TODO: move to util.h, possibly internal
  // return a 4-tuple (1,2,3,4) with:
  // 1 - total # of individual site-site interactions between the specified kaps
  //     and fgs
  // 2 - total # of bead pairs that make any site-site interaction
  // 3 - total # of kaps that contact any FG
  // 4 - sum of # of individual fg chains contacted by each kap
  boost::tuple<double, double, double, double> get_interactions_and_interacting(
      const ParticlesTemp &kaps, const base::Vector<ParticlesTemps> &fgs) const;


  // TODO: move to util.h, possibly internal
  /** get the number of site-site interactions between two particles */
  int get_number_of_interactions(Particle *a, Particle *b) const;

  /**
     for particles ps, returns the distribution along the z axis, in regions
     z0 = [top...) ;
     z1 = [0..top) ;
     z2 = [-top..top) ;
     z3 = (...top)
     with top being top of slab if slab exists,
     or 25% of box size if box exists but not slab, or 1.0 otherwise.


     @param ps the particles

     @return a tuple <z0,z1,z2,z3> with counts of particles from ps
             in z0, z1, z2 and z3 regions
   */
  boost::tuple<int, int, int, int>
    get_z_distribution(const ParticlesTemp& ps) const;


 public:
  IMP_OBJECT_METHODS(Statistics);


};

inline IMPNPCTRANSPORTEXPORT boost::timer create_boost_timer() {
  return boost::timer();
}

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_STATISTICS_H */
