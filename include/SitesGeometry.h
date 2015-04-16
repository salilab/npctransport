/**
 *  \file SitesGeometry.h
 *  \brief Geometry of sites on particle surfaces
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_SITES_GEOMETRY_H
#define IMPNPCTRANSPORT_SITES_GEOMETRY_H

#include "npctransport_config.h"
#include <IMP/Pointer.h>
#include <IMP/core/XYZR.h>
#include <IMP/core/Typed.h>
#include <IMP/display/particle_geometry.h>
#include <IMP/generic.h>
#include <IMP/algebra/vector_search.h>
#include <IMP/atom/estimates.h>
#include <boost/unordered_set.hpp>
#include "internal/sites.h"

#include <boost/array.hpp>

IMPNPCTRANSPORT_BEGIN_NAMESPACE


/** Show the sites. */
class IMPNPCTRANSPORTEXPORT SitesGeometry : public core::XYZRGeometry {
  algebra::Vector3Ds sites_;

 public:
  SitesGeometry(Particle *p, algebra::Vector3Ds sites)
      : core::XYZRGeometry(p), sites_(sites) {}
  virtual IMP::display::Geometries get_components() const IMP_OVERRIDE;
  IMP_OBJECT_METHODS(SitesGeometry);
};

/** Show the sites. */
class IMPNPCTRANSPORTEXPORT TypedSitesGeometry
    : public display::SingletonsGeometry {
  boost::unordered_map<core::ParticleType, algebra::Vector3Ds> sites_;

 public:
  TypedSitesGeometry(SingletonContainer *sc)
      : display::SingletonsGeometry(sc) {}
  void set_sites(core::ParticleType t, algebra::Vector3Ds s) {
    // std::cout << t << " gets " << s.size() << std::endl;
    sites_[t] = s;
  }
  virtual IMP::display::Geometries get_components() const IMP_OVERRIDE;
  IMP_OBJECT_METHODS(TypedSitesGeometry);
};

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SITES_GEOMETRY_H */
