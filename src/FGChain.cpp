/**
 *  \file FGChain.cpp
 *  \brief description.
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/FGChain.h>
#include <IMP/npctransport/HarmonicSpringSingletonScore.h>
#include <IMP/npctransport/linear_distance_pair_scores.h>
#include <IMP/npctransport/RelaxingSpring.h>
#include <IMP/npctransport/Scoring.h>
#include <IMP/npctransport/SimulationData.h>
#include <IMP/npctransport/internal/TAMDChain.h>
#include <IMP/npctransport/internal/npctransport.pb.h>

#include <IMP/atom/constants.h>
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
      {                                                          \
        var = fg_data.name().value();                            \
      }                                                          \
    else                                                         \
      { var = default_value; }                                   \
  }

/***************** FGChain methods ************/

//!  create the bonds restraint for the chain beads
void FGChain::update_bonds_restraint(Scoring const* scoring_manager)
{
  // TODO: this currently cannot work for more than two calls cause of
  // ExclusiveConsecutivePairContainer - will need to switch to
  // ConsecutivePairContainer or make a different design to solve this
  std::string name = get_root().get_particle()->get_name();
  if(bonds_restraint_){
     IMP_USAGE_CHECK(!bonds_restraint_->get_is_shared(),
                    "bonds restraint is supposed to become invalidated so it"
                  " shouldn't be owned by anyone else at this point");
      bonds_restraint_= nullptr; // invalidate to release old bead_pairs
      bonds_score_= nullptr;
  }
  boost::tie(bonds_restraint_, bonds_score_) =
    scoring_manager->create_backbone_restraint(rest_length_factor_,
					       backbone_k_,
					       this->get_beads(),
                                               name);
}

Restraints
FGChain::get_chain_restraints
(Scoring const* scoring_manager)
{
  // TODO: add support for a dynamic chain?
  IMP_USAGE_CHECK(root_, "Chain not initialized");
  if(!bonds_restraint_){ // TODO: fix for dynamic chain topology?
    update_bonds_restraint(scoring_manager);
  }
  return Restraints(1,bonds_restraint_);
}

//! update the rest length of the chain with rlf
void
FGChain::set_rest_length_factor
(double rlf)
{
  IMP_USAGE_CHECK(rlf>0.0, "bonds rest length factor should be positive");
  rest_length_factor_= rlf;
  if(bonds_score_){
    // this is ugly but just to support the old LinearWellPairScore properly
    LinearWellPairScore* lwps=
      dynamic_cast<LinearWellPairScore*>(bonds_score_.get());
    if(lwps != nullptr) {
      lwps->set_rest_length_factor(rlf);
    }
  }
  if(get_beads().size()==0) {
    return;
  }
  if(!RelaxingSpring::get_is_setup(get_bead(0))) {
    return;
  }
  // If RelaxingSpring decorator (e.g. for harmonic spring score) -
  // update decorator values - note this will affect even an existing scoring function
  // unlike with the older LinearWellPairScore, which stored the rest length factor internally
  unsigned int n(get_number_of_beads());
  for(unsigned int i=0; i < n-1; i++){
    IMP_USAGE_CHECK(RelaxingSpring::get_is_setup(get_bead(i)),
                    "If first bead in chain is decorated with a relaxing"
                    " spring, then all beads except last are");
    //    std::cout << "FGChain::set_rest_length_factor setting rest length factor to " << rlf << std::endl;
    RelaxingSpring rs_i(get_bead(i));
    rs_i.set_equilibrium_rest_length_factor(rlf);
  }
}


//! set the force constant between consecutive chain beads
void FGChain::set_backbone_k(double k)
{
  backbone_k_= k;
  // If score is already set, update it
  if(bonds_score_){
    // this is ugly but just to support the old LinearWellPairScore and the
    // new HarmonicSpringSingletonScore scores properly - in a perfect world
    // they derive from the same class with common set_backbone_k method
    {
      LinearWellPairScore* lwps=
        dynamic_cast<LinearWellPairScore*>(bonds_score_.get());
      if (lwps != nullptr) {
        lwps->set_k(k);
        return;
      }
    }
    {
      HarmonicSpringSingletonScore* hsss=
        dynamic_cast<HarmonicSpringSingletonScore*>(bonds_score_.get());
      if (hsss != nullptr) {
        // k is interpreted as k2 (force to bring spring rest length
        // to equilibrium) while maintaining the k1/k2 ratio constant
        // where k1 is the actual force constant on particles away from
        // rest length (which we generally assume to be tight)
        IMP_USAGE_CHECK(hsss->get_k2()*hsss->get_k1() != 0,
                        "expected the spring to already be set by now");
        double k1_to_k2= hsss->get_k1()/hsss->get_k2();
        hsss->set_k1(k1_to_k2 * k);
        hsss->set_k2(k);
        return;
      }
    }
    IMP_USAGE_CHECK(false, "bonds_score_ type couldn't be determined");
  }
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

  // append appropriate suffixes for type fields of chain particle, if specified
  // also update name, if current name is equal to the type
  void update_chain_particle_type_suffixes
  ( FGChain* fgc,
    const ::npctransport_proto::Assignment_FGAssignment &fg_data )
  {
    unsigned int n1 = fg_data.type_suffix_list_size();
    unsigned int n2 = fgc->get_number_of_beads();
    if ( n1 == 0 ) {
      return;
    }
    //    std::cout << "n1,n2: " << n1 << "," << n2 << std::endl;
    IMP_ALWAYS_CHECK(n1==n2,
                     "Size of list of type suffixes in FG chain assignment"
                     " should be equal to either zero or to the chain size",
                     IMP::ValueException);

    std::string type_prefix=fg_data.type();
    for(unsigned int i=0; i<n1; i++){
      std::string type_suffix=fg_data.type_suffix_list(i);
      core::ParticleType type(type_prefix + type_suffix);
      IMP::Particle* p=fgc->get_bead(i);
      IMP_USAGE_CHECK(core::Typed::get_is_setup(p),
                      "FG chain particle is expected to be typed");
      core::Typed(p).set_type(type);
      if(p->get_name()==type_prefix){
        p->set_name(type.get_string());
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
  Pointer<FGChain> ret_chain = nullptr;

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
    int n_levels = ceil(log(static_cast<float>(n))
                        /log(static_cast<float>(d))); // log_d{n}
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
    // Non-TAMD:
    Particles P;
    for (int i= 0; i<n; i++) {
      P.push_back( pf->create() );
    }
    std::string root_name = type.get_string(); // + " chain_root";
    Particle* root = atom::Hierarchy::setup_particle
      ( new Particle( sd->get_model(), root_name ), P );
    core::Typed::setup_particle( root, type );
    ret_chain = new FGChain(root);
    if(sd->get_is_backbone_harmonic()){
      // set springs between consecutive chain beads from "from" beads
      double rest_length_factor(fg_data.rest_length_factor().value());
      double tau_fs(sd->get_backbone_tau_ns()*(1e+6));
      double rest_length_diffusion_coefficient=
        atom::get_kt(sd->get_temperature_k())/(tau_fs*sd->get_scoring()->get_default_backbone_k()); // diffuse by kT/K per tau
      for (int i= 0; i<n-1; i++) {
	RelaxingSpring::setup_particle
	  (P[i],
	   P[i]->get_index(),
	   P[i+1]->get_index(),
	   rest_length_factor,
	   rest_length_diffusion_coefficient); // TODO: set this parameter properly
      }
    }
  }

  // update individual particle types if specified
  update_chain_particle_type_suffixes(ret_chain, fg_data);

  // add chain beads, etc. to sd under parent
  add_chain_to_sd(ret_chain, sd, parent, fg_data);

  // add scoring info
  ret_chain->set_rest_length_factor( fg_data.rest_length_factor().value() );
  ret_chain->set_backbone_k( sd->get_scoring()->get_default_backbone_k() );
  sd->get_scoring()->add_chain_restraints
    ( ret_chain );

  return ret_chain.release();
}


// create a newly allocated chain structure from a root of an FG nup
// (by adding its ordered leaves)
FGChain* get_fg_chain(atom::Hierarchy root){
  IMP_NEW(FGChain, ret, ( root.get_particle(), 0.0, 1.0, "chain copy %1%" ) );
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
