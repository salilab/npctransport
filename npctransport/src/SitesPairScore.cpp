/**
 *  \file DistancePairScore.cpp
 *  \brief A Score on the distance between a pair of particles.
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
 */

#include <IMP/npctransport/SitesPairScore.h>
#include <IMP/core/XYZR.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/generic.h>
#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/to_tuple.hpp>
#include <boost/preprocessor/facilities/apply.hpp>
#include <boost/preprocessor/seq/elem.hpp>
#include <boost/preprocessor/tuple/to_seq.hpp>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/**
   A score between two spherical particles that contain a fixed set
   of interaction sites,  sites0 and sites1 (for first and second particle,
   resp.).

   Note that for a specific pair of particles, each particle might have
   a different reference frame (rigid body translation and rotation),
   which is applied to the sites list upon score evaluation.

   @param range          range of site specific attraction
   @param k_attraction   site specific attraction coefficient
   @param range_nonspec_attraction range for non-specific attraction
                                   between particles that contain the sites
   @param k_nonspec_attraction     non-specific attraction coefficient
   @param k_repulsion    repulsion coefficient between penetrating particles
   @param sites0         list of sites on the first particle
   @param sites1         list of sites on the second particle
*/
SitesPairScore::SitesPairScore(double range, double k_attraction,
                               double range_nonspec_attraction,
                               double k_nonspec_attraction,
                               double k_repulsion,
                               const algebra::Vector3Ds &sites0,
                               const algebra::Vector3Ds &sites1)
  : P(k_repulsion, range_nonspec_attraction, k_nonspec_attraction, "Sites %1%"),
    range_(range), k_(k_attraction) {
  // store the big set of sizes in nnsites_
  if (sites0.size() > sites1.size()) {
    sites_first_=false;
    sites_=sites1;
    nnsites_=sites0;
  } else {
    sites_first_=true;
    sites_=sites0;
    nnsites_=sites1;
  }
  nn_= new algebra::NearestNeighbor3D(nnsites_);
}

void SitesPairScore::do_show(std::ostream & o) const {
  P::do_show( o );
  o << " Site specific range and attraction "
    << range_ << ", " << k_ << std::endl;
  o << "Sites 0: " << (sites_first_ ? sites_ : nnsites_) << std::endl;
  o << "Sites 1: " << (sites_first_ ? nnsites_ : sites_) << std::endl;
}



ParticlesTemp SitesPairScore
::get_input_particles(Particle*p) const {
  return ParticlesTemp(1, p);
}

ContainersTemp SitesPairScore
::get_input_containers(Particle*) const {
  return ContainersTemp();
}

// the sites of each particle are transformed to a common frame of reference
// (using the reference frame of each particle), and the site-specific
// attraction, and the inter-particle non specific attraction and repulsion
// are evaluated and summed.
inline double SitesPairScore
::evaluate_index(Model *m, const ParticleIndexPair& pp,
                 DerivativeAccumulator *da) const {
  IMP_OBJECT_LOG;

  // I. evaluate non-specific attraction and repulsion between
  //    parent particles before computing for specific sites :
  double soft=P::evaluate_index(m, pp, da);
  LinearInteractionPairScore::EvaluationCache
    lips_cache = this->get_evaluation_cache();

  // II. Return if parent particles are out of site-specific interaction range
  //     using cache to avoid some redundant calcs
  double particles_delta = std::sqrt( lips_cache.particles_delta_squared );
  double sum_radii = lips_cache.sum_particles_radii;
  double surface_delta = particles_delta - sum_radii; // between balls surface
  IMP_LOG(TERSE, "Surface_delta " << surface_delta << " ; range "
          << range_ << std::endl ) ; // TODO: VERBOSE
  if(surface_delta > range_){
    IMP_LOG(TERSE, "Sites contribution is 0.0 and soft sphere is "
            << soft << std::endl); // TODO: VERBOSE
    return soft; // can still contribute if non-specific range is longer
  }

  // III. evaluate site-specific contributions :
  // bring sites_ to the frame of reference of nn_sites_ and nn_
  core::RigidBody rb0(m, pp[0]); // first rigid particle
  core::RigidBody rb1(m, pp[1]); // second rigid-body particle
  // make sure rb0 is associated with sites_, and rb1 with nn_sites_:
  if (!sites_first_)
    std::swap(rb0, rb1);
  algebra::Transformation3D tr0
      = rb0.get_reference_frame().get_transformation_to();
  algebra::Transformation3D tr1
      = rb1.get_reference_frame().get_transformation_to();
  algebra::Rotation3D itr0
      = tr0.get_rotation().get_inverse();
  algebra::Rotation3D itr1
      = tr1.get_rotation().get_inverse();
  // bring sites_ to the correct orientation relative to nn_sites_[j]
  algebra::Transformation3D
    relative= tr1.get_inverse() *tr0;
  // sum over specific interactions between all pairs of sites:
  double sum= 0;
  for (unsigned int i=0; i< sites_.size(); ++i) {
    // filter to evaluate only sites within range of attraction:
    algebra::Vector3D trp=relative.get_transformed(sites_[i]);
    Ints nn= nn_->get_in_ball(trp, range_);
    for (unsigned int j=0; j < nn.size(); ++j) {
      //double d2=algebra::get_squared(range_);
      sum+= internal::evaluate_one_site_2(k_, range_,
                                        rb0, rb1,
                                        tr0, tr1, itr0, itr1, sites_[i],
                                        nnsites_[nn[j]], da);
    }
  }
  IMP_LOG(TERSE, "Sites contribution is " << sum << " and soft sphere is "
          << soft << std::endl); // TODO: VERBOSE
  return sum+ soft; // specific + non-specific
}


IMP::display::Geometries SitesGeometry::get_components() const {
  display::Geometries ret;
  algebra::ReferenceFrame3D rf=
    core::RigidBody(get_particle()).get_reference_frame();
  double r= .1*core::XYZR(get_particle()).get_radius();
  for (unsigned int i=0; i< sites_.size(); ++i) {
    IMP_NEW(display::SphereGeometry, g,
            (algebra::Sphere3D(rf.get_transformation_to()
                               .get_transformed(sites_[i]),
                               r)));
    g->set_color(display::Color(1,0,0));
    ret.push_back(g);
  }
  return ret+core::XYZRGeometry::get_components();
}

void SitesGeometry::do_show(std::ostream &out) const {
  out << " sites: " << Showable(sites_) << std::endl;
}



IMP::display::Geometries TypedSitesGeometry::get_components() const {
  display::Geometries ret;
  IMP_FOREACH_SINGLETON(get_container(), {
      core::ParticleType t= core::Typed(_1).get_type();
      IMP_NEW(SitesGeometry, g, (_1, sites_.find(t)->second));
      ret.push_back(g);
    });
  return ret;
}

void TypedSitesGeometry::do_show(std::ostream &) const {
  //SingletonsGeometry::do_show(out);
}


#define IMP_HS(na, nb)                                  \
  else if (sites0.size()==na                            \
      && sites1.size()==nb) {                           \
    typedef TemplateSitesPairScore<na, nb, true> TSPS;  \
    IMP_NEW(TSPS, ps, (range, ka, kr, sites0, sites1)); \
    ppr->set_score(value, ps.get());                    \
  }

#if IMP_BUILD==IMP_FAST
#define IMP_S0 (1)(2)(3)(4)(5)(6)(7)(8)(9)(10)(11)(12)(13)\
  (14)(15)(16)(17)(18)(19)(20)
#else
#define IMP_S0 (1)(2)(3)(4)(5)
#endif
#define IMP_CALL(a,bb) a(BOOST_PP_SEQ_ELEM(0, bb),BOOST_PP_SEQ_ELEM(1, bb))

#define IMP_MACRO(r, product) IMP_CALL(IMP_HS,product)

/*
// #define IMP_ADD_CASE(na,nb)  else if (sites0.size()==na \
//                                       && sites1.size()==nb) {   \
//     typedef TemplateSitesPairScore<na, nb, true> TSPS;          \
//     IMP_NEW(TSPS, ps, (range, ka, rangena, kna, kr, sites0, sites1));   \
//     ppr->set_score(value, ps.get());                            \
//   }
*/

/**
   @param rangea   range of specific attraction
   @param ka       coefficient for specific attraction
   @param rangena  range for non-specific attraction
   @param kna      coefficient for non-specific attraction
   @param kr       coefficient for repulsion between penetrating particles
   @param sites0
   @param sites1
   @param value
   @param ppr
*/
void set_sites_score(double rangea, double ka,
                     double rangena, double kna,
                     double kr,
                     const algebra::Vector3Ds &sites0,
                     const algebra::Vector3Ds &sites1,
                     int value,
                     container::PredicatePairsRestraint *ppr) {
  // use the point-location based one if appropriate
  if ((sites0.size()>=4 && sites1.size() >=4
       && (sites0.size()>.5*sites1.size()
           || sites1.size() > .5*sites0.size()))) {
    IMP_NEW(SitesPairScore, ps, (rangea, ka, rangena, kna, kr,
                                 sites0, sites1));
    ppr->set_score(value, ps.get());
  }
  // too slow (TODO: what is too slow?)
  // TODO: check if the approach of the next clause (originally IMP_ADD_CASE)
  //       can accelerate computations by having a template for spefified
  //       site sizes
  //BOOST_PP_SEQ_FOR_EACH_PRODUCT(IMP_MACRO, (IMP_S0)(IMP_S0))
  //  IMP_ADD_CASE(1,15)
  else if (sites0.size()==1
           && sites1.size()==15) {
    typedef TemplateSitesPairScore<1, 15, true> TSPS;
    IMP_NEW(TSPS, ps, (rangea, ka, rangena, kna, kr, sites0, sites1));
    ppr->set_score(value, ps.get());
  }  else {
    IMP_NEW(SitesPairScore, ps, (rangea, ka, rangena, kna, kr, sites0, sites1));
    ppr->set_score(value, ps.get());
  }
}


Restraints SitesPairScore
::create_current_decomposition(Model *m,
                               const ParticleIndexPair &pi) const {
  Restraints ret;
  if (evaluate_index(m, pi, nullptr) < 0) {
    return Restraints(1, IMP::internal::create_tuple_restraint(this, m, pi));
  } else {
    return Restraints();
  }
}
IMPNPCTRANSPORT_END_NAMESPACE
