/**
 *  \file TunnelGeometry.cpp
 *  \brief XXXX.
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 *
 */

#include "IMP/npctransport/SlabGeometry.h"
#include <IMP/display/geometry.h>
#include <IMP/algebra/Cylinder3D.h>
#include <IMP/constants.h>
#include <IMP/algebra/Plane3D.h>
#ifdef IMP_NPCTRANSPORT_USE_IMP_CGAL
#include <IMP/cgal/internal/polyhedrons.h>
#endif
#include <boost/tuple/tuple.hpp>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/*
namespace {
  algebra::Vector3D flip(int l, algebra::Vector3D pt) {
    if (l==0) {}
    else if (l==1) {
      pt[0]=-pt[0];
    } else if (l==2) {
      pt[0]=-pt[0];
      pt[1]=-pt[1];
    } else {
      pt[1]=-pt[1];
    }
    return pt;
  }
}
*/
#ifdef IMP_NPCTRANSPORT_USE_IMP_CGAL
namespace {
const std::pair<algebra::Vector3Ds, Ints> get_mesh(double height, double radius,
                                                   double width) {
  Vector<algebra::Plane3D> slab;
  slab.push_back(algebra::Plane3D(height, algebra::Vector3D(0, 0, 1)));
  slab.push_back(algebra::Plane3D(height, algebra::Vector3D(0, 0, -1)));
  Vector<algebra::Plane3D> tunnel;
  unsigned int npoints = 40;
  for (unsigned int i = 0; i < npoints; ++i) {
    double f = 2.0 * i * PI / (npoints - 1.0);
    double x = std::cos(f);
    double y = std::sin(f);
    tunnel.push_back(algebra::Plane3D((radius), algebra::Vector3D(x, y, 0)));
  }

  algebra::Plane3Ds empty;
  // empty.push_back(algebra::Plane3D(-1, algebra::Vector3D(1, 0, 0)));
  // empty.push_back(algebra::Plane3D(1, algebra::Vector3D(-1, 0, 0)));

  algebra::Vector3D c = algebra::get_ones_vector_d<3>() * width / 2.0;
  algebra::BoundingBox3D bb(-c, c);
  return cgal::internal::get_polyhedron_indexed_facets(bb, slab, tunnel);
}
}

SlabSurfaceGeometry::SlabSurfaceGeometry(double height, double radius,
                                         double width)
    : SurfaceMeshGeometry(get_mesh(height, radius, width), "Slab") {}

#else
SlabSurfaceGeometry::SlabSurfaceGeometry(double height, double radius,
                                         double width)
    : SurfaceMeshGeometry(std::pair<algebra::Vector3Ds, Ints>(), "Slab") {
  IMP_UNUSED(height);
  IMP_UNUSED(radius);
  IMP_UNUSED(width);
}
#endif

SlabWireGeometry::SlabWireGeometry(double height, double radius, double width)
    : Geometry("Slab"), height_(height), radius_(radius), length_(width) {}

display::Geometries SlabWireGeometry::get_components() const {
  display::Geometries ret;
  const int n = 30;
  for (int i = 1; i <= n; ++i) {
    double f = static_cast<double>(i) / n;
    double f1 = static_cast<double>(i - 1) / n;
    algebra::Vector3D v0(radius_ * sin(2 * IMP::PI * f),
                         radius_ * cos(2 * IMP::PI * f), height_ / 2.0);
    algebra::Vector3D v1(radius_ * sin(2 * IMP::PI * f1),
                         radius_ * cos(2 * IMP::PI * f1), height_ / 2.0);
    ret.push_back(new display::SegmentGeometry(algebra::Segment3D(v0, v1)));
    algebra::Vector3D v0m(radius_ * sin(2 * IMP::PI * f),
                          radius_ * cos(2 * IMP::PI * f), -height_ / 2.0);
    algebra::Vector3D v1m(radius_ * sin(2 * IMP::PI * f1),
                          radius_ * cos(2 * IMP::PI * f1), -height_ / 2.0);
    ret.push_back(new display::SegmentGeometry(algebra::Segment3D(v0m, v1m)));
  }
  return ret;
}

IMPNPCTRANSPORT_END_NAMESPACE
