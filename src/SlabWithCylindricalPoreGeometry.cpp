/**
 *  \file SlabWithCylindricalPoreGeometry.cpp
 *  \brief XXXX.
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 *
 */

#include "IMP/npctransport/SlabWithCylindricalPoreGeometry.h"
#include <IMP/display/geometry.h>
#include <IMP/algebra/Cylinder3D.h>
#include <IMP/constants.h>
#include <IMP/algebra/Plane3D.h>
#ifdef IMP_NPCTRANSPORT_USE_IMP_CGAL
#include <IMP/cgal/internal/polyhedrons.h>
#endif
#include <boost/tuple/tuple.hpp>
#include <algorithm>
#include <set>

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
  // Get a mesh representation of a cylindrical pore in a square slab.
  // Height specifies the slab thickness, radius specifies the cylidner radius,
  // width specifies the edge size of the slab.
  const std::pair<algebra::Vector3Ds, Ints>
  get_mesh
  (double height, double radius, double width)
  {
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

    bool less_vector3D_lexigoraphically
    (algebra::Vector3D v0, algebra::Vector3D v1)
    {
      return std::lexicographical_compare
        (v0.begin(), v0.end(),
         v1.begin(), v1.end());
    }

  // sort the segment vertices such that the first vector
  // is smaller or equal to the second one, lexigoraphically
  algebra::Segment3D sort_segment_lexigraphically
  (algebra::Segment3D s)
  {
    algebra::Vector3D v0=s.get_point(0);
    algebra::Vector3D v1=s.get_point(1);
    if(less_vector3D_lexigoraphically(v0,v1))
      {
        return s;
      }
    else
      {
        return algebra::Segment3D(v1,v0);
      }
  }

  // compares two segmenets lexigographically,
  // after sorting s0 and s1 vertices, such that
  // their directionality does not matter
  double less_segment_lexicographically_undirectional
  (algebra::Segment3D s0, algebra::Segment3D s1)
  {
    s0=sort_segment_lexigraphically(s0);
    s1=sort_segment_lexigraphically(s1);
    if( less_vector3D_lexigoraphically
        (s0.get_point(0), s1.get_point(0)) )
      {
        return true;
      }
    if( less_vector3D_lexigoraphically
        (s1.get_point(0), s0.get_point(0)) )
      {
        return false;
      }
    return less_vector3D_lexigoraphically
      (s0.get_point(1), s1.get_point(1));
  }

  struct CompareSegmentsLexicographicallyUndirectional
  {
    bool operator() (algebra::Segment3D s0, algebra::Segment3D s1) const {
      return less_segment_lexicographically_undirectional(s0, s1);
    }
  };

};

SlabWithCylindricalPoreSurfaceGeometry::SlabWithCylindricalPoreSurfaceGeometry(double height, double radius,
                                         double width)
    : SurfaceMeshGeometry(get_mesh(height, radius, width), "SlabWithCylindricalPore")
{
  boost::tie(vertices_, faces_) =
    get_mesh(height, radius, width);
}

namespace {
};

display::Geometries SlabWithCylindricalPoreSurfaceGeometry::get_components() const {
  display::Geometries ret;
  typedef
    std::set<algebra::Segment3D, CompareSegmentsLexicographicallyUndirectional>
    t_segments_set;
  t_segments_set segments;
  algebra::Vector3Ds cur;
  for (unsigned int i = 0; i < faces_.size(); ++i) {
    if (faces_[i] == -1) {
      for(unsigned int j=0; j<cur.size(); j++){
        for(unsigned int k=j+1; k<cur.size(); k++){
          segments.insert(algebra::Segment3D(cur[j],cur[k]));
        } // k
      } // j
      cur.clear();
    } else {
      IMP_USAGE_CHECK(vertices_.size() > static_cast<unsigned int>(faces_[i]),
                      "Out of range vertex: " << faces_[i]);
      cur.push_back(vertices_[faces_[i]]);
    }
  } // i
  for(t_segments_set::const_iterator it=segments.begin();
      it!=segments.end();
      it++)
    {
      ret.push_back
        ( new display::SegmentGeometry(*it) );
    }

  return ret;
}


#else
SlabWithCylindricalPoreSurfaceGeometry::SlabWithCylindricalPoreSurfaceGeometry(double height, double radius,
                                         double width)
    : SurfaceMeshGeometry(std::pair<algebra::Vector3Ds, Ints>(), "SlabWithCylindricalPore") {
  IMP_UNUSED(height);
  IMP_UNUSED(radius);
  IMP_UNUSED(width);
}
#endif

SlabWithCylindricalPoreWireGeometry::SlabWithCylindricalPoreWireGeometry(double height, double radius, double length)
  : Geometry("SlabWithCylindricalPore"), height_(height), radius_(radius), length_(length) {}

display::Geometries SlabWithCylindricalPoreWireGeometry::get_components() const {
  display::Geometries ret;
  const int n= 30;
  for (int i= 1; i <= n; ++i) {
    double f= static_cast<double>(i) / n;
    double f1= static_cast<double>(i - 1) / n;
    double theta= 2 * IMP::PI * f;
    double theta1= 2 * IMP::PI * f1;
    algebra::Vector3D v0(radius_ * sin(theta),
                         radius_ * cos(theta),
                         height_ / 2.0);
    algebra::Vector3D v1(radius_ * sin(theta1),
                         radius_ * cos(theta1),
                         height_ / 2.0);
    ret.push_back(new display::SegmentGeometry(algebra::Segment3D(v0, v1))); // cylinder top
    algebra::Vector3D v0m(radius_ * sin(theta),
                          radius_ * cos(theta),
                          -height_ / 2.0);
    algebra::Vector3D v1m(radius_ * sin(theta1),
                          radius_ * cos(theta1),
                          -height_ / 2.0);
    ret.push_back(new display::SegmentGeometry(algebra::Segment3D(v0m, v1m))); // cylinder bottom
    ret.push_back(new display::SegmentGeometry(algebra::Segment3D(v0, v0m))); // cylinder curved face
    if(i%2==0) {
      continue;
    }
    double isin= 1.0 / (sin(theta)+0.0000001);
    double icos= 1.0 / (cos(theta)+0.0000001);
    double d= 0.5 * length_ * std::min(std::abs(isin), std::abs(icos));
    algebra::Vector3D v2(d * sin(theta),
                          d * cos(theta),
                          0.5 * height_);
    algebra::Vector3D v2m(d * sin(theta),
                           d * cos(theta),
                           -0.5 * height_);
    ret.push_back(new display::SegmentGeometry(algebra::Segment3D(v0, v2))); // ray over slab top face
    ret.push_back(new display::SegmentGeometry(algebra::Segment3D(v0m, v2m))); // ray over slab bottom face
  }
  // add top and bottom slab rectangles:x
  for(double sign=-0.5; sign<=0.5; sign+=1){
    algebra::Vector3D vNE(0.5*length_, 0.5*length_, sign*height_);
    algebra::Vector3D vNW(-0.5*length_, 0.5*length_, sign*height_);
    algebra::Vector3D vSW(-0.5*length_, -0.5*length_, sign*height_);
    algebra::Vector3D vSE(0.5*length_, -0.5*length_, sign*height_);
    ret.push_back(new display::SegmentGeometry(algebra::Segment3D(vNE, vNW))); // slab bottom face
    ret.push_back(new display::SegmentGeometry(algebra::Segment3D(vNW, vSW))); // slab bottom face
    ret.push_back(new display::SegmentGeometry(algebra::Segment3D(vSW, vSE))); // slab bottom face
    ret.push_back(new display::SegmentGeometry(algebra::Segment3D(vSE, vNE))); // slab bottom face
  }
  return ret;
}

IMPNPCTRANSPORT_END_NAMESPACE
