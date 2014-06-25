/**
 * \file FGChain.h
 * \brief creating TAMD chains
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

//TODO: turn FGChain, TAMDChain onto a Decorator? seems to make sense
//TODO: use "Representation" decorator?

#ifndef IMPNPCTRANSPORT_FG_CHAIN_H
#define IMPNPCTRANSPORT_FG_CHAIN_H

#include "npctransport_config.h"
#include "linear_distance_pair_scores.h"
#include "npctransport_proto.fwd.h"
//#include "internal/npctransport.pb.h"

#include <IMP/atom/Hierarchy.h>
#include <IMP/base/Object.h>
#include <IMP/base/nullptr.h>
#include <IMP/container/ConsecutivePairContainer.h>
#include <IMP/display/Color.h>


IMPNPCTRANSPORT_BEGIN_NAMESPACE

class SimulationData;

/**
   a chain with a root

    root - root of hierarchy
    beads - fine chain particles
*/
class IMPNPCTRANSPORTEXPORT FGChain : public IMP::base::Object {
public:
 private:
  // the root particle in the chain hierarchy
  base::PointerMember<Particle> root_;

  // the restraint on the chain bonds
  base::PointerMember<Restraint> bonds_restraint_;

  // the score acting on consecutive chain beads in bonds_restraint_
  base::PointerMember<LinearWellPairScore> bonds_score_;

  // container for consecutive (bonded) chain beads
  base::PointerMember<container::ExclusiveConsecutivePairContainer>
    bead_pairs_;


 private:

  // TODO: this currently cannot work for more than two calls cause of
  // ExclusiveConsecutivePairContainer - will need to switch to
  // ConsecutivePairContainer or make a different design to solve this
    /** recreate the bonds restraint for the chain beads based on the
        current chain topology
    */
  void update_bonds_restraint();

 public:

  /** initialized an FG chain with given backbone k and rest length factor,
      and a specified root whose leaves are the beads of the chain.

      @param root root of the chain (leaves are assumed beads).
                  can be null = to be added later
      @param backbone_k force constant between consecutive chain beads
      @param rest_length_factor equilibrium distance factor between
                        consecutive beads relative to the sum of their radii
      @param name chain object name

      \see LinearWellPairScore
  */
 FGChain(IMP::Particle* root,
         double backbone_k = 0.0,
         double rest_length_factor = 1.0,
         std::string name = "chain %1%")
   : base::Object(name),
    root_(root),
    bonds_restraint_(nullptr),
    bead_pairs_(nullptr)
      {
        IMP_USAGE_CHECK(rest_length_factor>0.0, "bonds rest length factor" <<
                        " should be positive");
        bonds_score_ = new LinearWellPairScore(rest_length_factor, backbone_k);
      }

    atom::Hierarchy get_root() const
      { return atom::Hierarchy(root_); }

 protected:
    /** set the root of the chain to this particle
        with the beads being the leaves of Hierarchy(p)

        @note it is assumed that p is decorated as atom::Hiererachy
    */
    void set_root(Particle* p){
      IMP_USAGE_CHECK(atom::Hierarchy::get_is_setup(p),
                      "root must be Hierarchy decorated");
      root_ = p;
    }

    /** set the root of the chain to this particle
        with the beads being the leaves of Hierarchy(p)
    */
    void set_root(atom::Hierarchy root){
      root_ = root.get_particle();
    }



 public:

  //! get the beads of the chain (assume valid root)
  IMP::ParticlesTemp get_beads() const
    { return core::get_leaves(get_root()); }

  //! get the i'th bead in the chain (assume valid root)
  IMP::Particle* get_bead(unsigned int i) const
    { return get_beads()[i]; }

  //! get the i'th bead index in the chain (assume valid root)
  IMP::ParticleIndex get_bead_index(unsigned int i) const
  { return get_beads()[i]->get_index(); }

  //! get the number of beads in the chain (assume valid root)
  unsigned int get_number_of_beads() const
  { return get_beads().size(); }

  /**
      Returns a restraint accosciated with internal interactions by this chain.
      Once this method has been called once, it is assumed that the
      chain topology remains static (the behavior of the restraint
      if the chain topology changes is undefined, e.g. if beads are added)

      @note this restraint is affected by future calls to set_rest_length_factor()
      and set_backbone_k(). However, it applies only to the topology of the chain
      at the time of call, if the chain topology changes, this methods should be
      called again.

      @note assumes that the chain has a valid root whose leaves are beads
  */
  virtual Restraints get_chain_restraints() {
    // TODO: add support for a dynamic chain?
    IMP_USAGE_CHECK(root_, "Chain not initialized");
    if(!bonds_restraint_){ // TODO: fix for dynamic chain topology?
      update_bonds_restraint();
    }
    return Restraints(1,bonds_restraint_);
  }

    /** set the equilibrium distance factor between consecutive beads
        relative to the sum of their radii

        @note This affects also restraints previously returned by
              get_chain_restraints()

        \see LinearWellPairScore
    */
  void set_rest_length_factor(double rlf){
    IMP_USAGE_CHECK(rlf>0.0, "bonds rest length factor should be positive");
      bonds_score_->set_rest_length_factor(rlf);
  }

    /** set the force constant between consecutive chain beads

        @note This affects also restraints previously returned by
              get_chain_restraints()

        \see LinearWellPairScore
     */
  void set_backbone_k(double k) {
      bonds_score_->set_k(k);
  }

  //! get the equilibrium distance factor between consecutive beads relative
  //! to the sum of their radii
  double get_rest_length_factor(){
    return bonds_score_->get_rest_length_factor();
  }

  //! get the force constant between consecutive chain beads
  double get_backbone_k() {
    return bonds_score_->get_k();
  }

  IMP_OBJECT_METHODS(FGChain);
};

IMP_OBJECTS(FGChain, FGChains);

/******************  utility methods ***************/

/**
   Create a chain particle hierarchy, associated with the model of sd,
   with restraint bonding consecutive particles added to sd, according to the
   parameters specified in fg_data.

   Notes:

   All nodes in the hierarchy share the same type, based on the type
   string fg_data.type().

   The rest length between two consecutive chain beads is
   fg_data.radius() * 2.0 * fg_data.rest_length_factor() and the
   spring constant is the simulation backbone_k parameter.

   If fg_data.is_tamd() is true, created a TAMD hierarchy, otherwise
   a simple parent + beads structure. In the TAMD case, the custom restraint
   are added to sd->get_scoring() and the tamd images are added to sd->root()

   @param[in,out] sd the simulation data whose model is associated with the
                  new chain. The chain is also added to the simulation data
                  scoring object.
   @param parent parent hierarchy to which chain is added
   @param[in] fg_data data about the FG chain
   @param[in] c        color of chain particles

   @return chain structure (with root and chain beads)

 */
FGChain* create_fg_chain
( IMP::npctransport::SimulationData *sd,
  atom::Hierarchy parent,
  const ::npctransport_proto::Assignment_FGAssignment &fg_data,
  display::Color c );


/**
   gets a chain structure from a root of an FG nup
   (by adding its ordered leaves)
*/
 FGChain* get_fg_chain(atom::Hierarchy root);

 /**
    gets a chain structure from a root of an FG nup
    (by adding its ordered leaves)

    @param p_root a particle that is assumed to be
                  Hierarchy decorated
 */
 FGChain* get_fg_chain(Particle* p_root);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_FG_CHAIN_H */
