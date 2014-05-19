/**
 *  \file FGChain.cpp
 *  \brief description.
 *
 *  Copyright 2007-2014 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/FGChain.h>
#include <IMP/npctransport/SimulationData.h>
#include <IMP/npctransport/internal/TAMDChain.h>
#include <IMP/npctransport/internal/npctransport.pb.h>

#include <IMP/atom/Diffusion.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/CenterOfMass.h>
#include <IMP/atom/Mass.h>
#include <IMP/atom/TAMDParticle.h>
#include <IMP/Pointer.h>
#include <IMP/nullptr.h>
#include <IMP/core/ChildrenRefiner.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/core/Typed.h>
#include <IMP/core/XYZR.h>
#include <IMP/display/Colored.h>

#include <boost/tuple/tuple.hpp>
#include <cmath>
#include <cstdio>
#include <sstream>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

// type - variable type to be declared e.g. double
// var - var name to be declared
// name - name in protobuf
// default_value - the default
#define GET_ASSIGNMENT_DEF(type, var, name, default_value)       \
  type var;                                                      \
  {                                                              \
    if(fg_data.has_##name() )                                    \
      { var = fg_data.name().value(); }                    \
    else                                                         \
      { var = default_value; }                                   \
  }

/***************** FGChain methods ************/

//!  create the bonds restraint for the chain beads
void FGChain::update_bonds_restraint()
{
  // TODO: this currently cannot work for more than two calls cause of
  // ExclusiveConsecutivePairContainer - will need to switch to
  // ConsecutivePairContainer or make a different design to solve this
  std::string name = get_root().get_particle()->get_name();
  if(bonds_restraint_){
     IMP_USAGE_CHECK(!bonds_restraint_->get_is_shared(),
                    "bonds restraint is supposed to become invalidated so it"
                  " shouldn't be owned by anyone else at this point");
      bonds_restraint_ = nullptr; // invalidate to release old bead_pairs
      IMP_USAGE_CHECK(!bead_pairs_->get_is_shared(),
                    "bead pairs consecutive pair container must be destroyed"
                    " before it is updated so it cannot be owned by others");
      bead_pairs_ = nullptr; // invalidate to destruct
  }
  bead_pairs_ = new IMP::container::ExclusiveConsecutivePairContainer
    (this->get_beads(), "Bonds %1% " + name + " consecutive pairs");
  bonds_restraint_ = container::create_restraint
    ( bonds_score_.get(), bead_pairs_.get(),  "Bonds " + name  );
}


/******************  internal utility methods ***************/

namespace {

  void add_chain_to_sd
    (FGChain const* chain,
     SimulationData* sd,
     atom::Hierarchy parent,
     const ::npctransport_proto::Assignment_FGAssignment &fg_data)
  {
    parent.add_child( chain->get_root() );
    if(fg_data.is_tamd()) {
      internal::TAMDChain const* tamd_chain =
        dynamic_cast<internal::TAMDChain const*>(chain);
      IMP_USAGE_CHECK(tamd_chain, "chain is expected to be TAMD chain*");
      Particles tamd_images = tamd_chain->get_tamd_images();
      for(unsigned int k = 0; k < tamd_images.size(); k++) {
        atom::Hierarchy h_k( tamd_images[k]);
        sd->get_root().add_child( h_k );
      }
    }
  }

} // anonymous namespace


/******************  utility methods ***************/

FGChain* create_fg_chain
( SimulationData* sd,
  atom::Hierarchy parent,
  const ::npctransport_proto::Assignment_FGAssignment &fg_data,
  display::Color c )
{
  base::Pointer<FGChain> ret_chain = nullptr;

  // set up factory for chain particles:
  core::ParticleType type(fg_data.type());
  int n = fg_data.number_of_beads().value();
  double radius = fg_data.radius().value();
  double D_factor = fg_data.d_factor().value();
  double angular_D_factor = sd->get_angular_d_factor();
  IMP_NEW( ParticleFactory, pf,
           (sd, radius, D_factor, angular_D_factor, c, type) );

  // create TAMD or non-TAMD particles P, within hierarchy of root:
  if(fg_data.is_tamd()) {
    int d = 2; // outdegree in TAMD hierarchy, TODO: parametrize
    int n_levels = ceil(log(n)/log(d)); // log_d{n}
    std::vector<double> T_factors(n_levels); // temperature scaling
    std::vector<double> F_factors(n_levels); // friction scaling
    std::vector<double> Ks(n_levels); // TAMD spring constant
    for(int i=0; i < n_levels; i++) { // i ~ increasing depth from root
      int level = n_levels - i; // level above leaves
      GET_ASSIGNMENT_DEF(double, t_coeff, tamd_t_factor_coeff, 1.0);
      GET_ASSIGNMENT_DEF(double, t_base, tamd_t_factor_base, 1.0);
      GET_ASSIGNMENT_DEF(double, f_coeff, tamd_f_factor_coeff, 1.0);
      GET_ASSIGNMENT_DEF(double, f_base, tamd_f_factor_base, 1.0);
      GET_ASSIGNMENT_DEF(double, k, tamd_k, 1.0);
      std::cout << "TAMD params: "
                << t_coeff << ", "  << t_base << ", "
                << f_coeff << ", "  << f_base << ", "
                << k << std::endl;
      T_factors[i] = t_coeff * pow(t_base,level-1);
      F_factors[i] = f_coeff * pow(f_base,level-1);
      Ks[i] = k;
    }
    ret_chain=
      internal::create_tamd_chain(pf, n, d, T_factors, F_factors, Ks);
  } else {
    Particles P;
    for (int i = 0; i < n; ++i) {
      P.push_back( pf->create() );
    }
    std::string root_name = type.get_string() + " chain_root";
    Particle* root = atom::Hierarchy::setup_particle
      ( new Particle( sd->get_model(), root_name ), P );
    core::Typed::setup_particle( root, type );
    ret_chain = new FGChain(root);
  }

  // add chain beads, etc. to sd under parent
  add_chain_to_sd(ret_chain, sd, parent, fg_data);

  // add scoring info
  ret_chain->set_rest_length_factor( fg_data.rest_length_factor().value() );
  ret_chain->set_backbone_k( sd->get_scoring()->get_default_backbone_k() );
  sd->get_scoring()->add_chain_restraints
    ( ret_chain );

  return ret_chain.release();
}
#endif

// gets a chain structure from a root of an FG nup
// (by adding its ordered leaves)
FGChain* get_fg_chain(atom::Hierarchy root){
  IMP_NEW(FGChain, ret, ( root.get_particle() ) );
  return ret.release();
}


FGChain* get_fg_chain(Particle* p_root)
{
  IMP_USAGE_CHECK(atom::Hierarchy::get_is_setup(p_root),
                  "p_root must be of type hierarchy");
  return get_fg_chain(atom::Hierarchy(p_root));
};


#undef GET_ASSIGNMENT_DEF

IMPNPCTRANSPORT_END_NAMESPACE
