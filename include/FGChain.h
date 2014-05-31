/**
 * \file FGChain.h
 * \brief creating TAMD chains
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

//TODO: turn FGChain, TAMDChain onto a Decorator? seems to make sense
//TODO: use "Representation" decorator?

#ifndef IMPNPCTRANSPORT_CREATING_TAMD_PARTICLES_H
#define IMPNPCTRANSPORT_CREATING_TAMD_PARTICLES_H

#include "npctransport_config.h"
#include "SimulationData.h"
#include "linear_distance_pair_scores.h"
#include "internal/npctransport.pb.h"
#include <IMP/atom/Hierarchy.h>
#include <IMP/base/Object.h>
#include <IMP/base/nullptr.h>
#include <IMP/container/ConsecutivePairContainer.h>
#include <IMP/display/Color.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/**
   a chain with a root

    root - root of hierarchy
    beads - fine chain particles
*/
class IMPNPCTRANSPORTEXPORT FGChain : public IMP::base::Object {
public:
 private:
  base::PointerMember<Particle> root_;
  base::PointerMember<Restraint> bonds_restraint_;
  base::PointerMember<LinearWellPairScore> bonds_score_;
  double backbone_k_;
  double rest_length_factor_;

 private:

    /** create the bonds restraint for the chain beads based on the
        internal values of backbone k and bond rest length factor
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
      bonds_score_(nullptr),
      backbone_k_(backbone_k),
      rest_length_factor_(rest_length_factor)
      {
        IMP_USAGE_CHECK(rest_length_factor_>0.0, "bonds rest length factor should be positive");
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
      return the restraint accosciated with internal interactions by this chain

      @param assumes that the chain has a valid root whose leaves are beads
  */
  virtual Restraints get_chain_restraints() {
    IMP_USAGE_CHECK(root_, "Chain not initialized");
    update_bonds_restraint();
    return Restraints(1,bonds_restraint_);
  }

    /** set the equilibrium distance factor between consecutive beads
        relative to the sum of their radii
        \see LinearWellPairScore
    */
  void set_rest_length_factor(double brlf){
    IMP_USAGE_CHECK(brlf>0.0, "bonds rest length factor should be positive");
    rest_length_factor_ = brlf;
    if(bonds_restraint_ && bonds_score_){
      // quicker if score exists
      bonds_score_->set_rest_length_factor(brlf);
    } else {
      if(root_){
        update_bonds_restraint();
      }
    }
  }

    /** set the force constant between consecutive chain beads
        \see LinearWellPairScore
     */
  void set_backbone_k(double k) {
    backbone_k_ = k;
    if(bonds_restraint_ && bonds_score_){
      // quicker if score exists
      bonds_score_.get()->set_k(k);
    } else {
      if(root_){
        update_bonds_restraint();
      }
    }
  }

    double get_rest_length_factor(){ return rest_length_factor_; }

      double get_backbone_k() { return backbone_k_; }

  IMP_OBJECT_METHODS(FGChain);
};


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
                  new chain. A chain restraint is added to the simulation data
                  scoring object, and the particle is added to the simulation data
                  diffusers list.
   @param parent parent hierarchy to which chain is added
   @param[in] fg_data data about the FG chain
   @param[in] c        color of chain particles

   @return chain structure (with root and chain beads)

 */
FGChain* create_fg_chain
( SimulationData *sd,
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

#endif /* IMPNPCTRANSPORT_CREATING_TAMD_PARTICLES_H */
