/**
 *  \file SitesGeometry.cpp
 *  \brief Geometry for sites on particle surfaces
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#include <IMP/npctransport/npctransport_config.h>
#include <IMP/npctransport/SitesGeometry.h>
#include <IMP/core/XYZR.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/generic.h>
#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/to_tuple.hpp>
#include <boost/preprocessor/facilities/apply.hpp>
#include <boost/preprocessor/seq/elem.hpp>
#include <boost/preprocessor/tuple/to_seq.hpp>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

IMP::display::Geometries SitesGeometry::get_components() const {
  display::Geometries ret;
  algebra::ReferenceFrame3D rf =
      core::RigidBody(get_particle()).get_reference_frame();
  double r = .1 * core::XYZR(get_particle()).get_radius();
  for (unsigned int i = 0; i < sites_.size(); ++i) {
    IMP_NEW(display::SphereGeometry, g,
            (algebra::Sphere3D
             ( rf.get_transformation_to().get_transformed
               (sites_[i].get_center()), r)
             )
            );
    g->set_color(display::Color(1, 0, 0));
    ret.push_back(g);
  }
  return ret + core::XYZRGeometry::get_components();
}

IMP::display::Geometries TypedSitesGeometry::get_components() const {
  display::Geometries ret;
  Model *m = get_container()->get_model();
  IMP_CONTAINER_FOREACH(SingletonContainer, get_container(), {
    core::ParticleType t = core::Typed(m, _1).get_type();
    IMP_NEW(SitesGeometry, g, (m->get_particle(_1), sites_.find(t)->second));
    ret.push_back(g);
  });
  return ret;
}

IMPNPCTRANSPORT_END_NAMESPACE
