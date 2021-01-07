/**
 *  \file SlabWithCylindricalPoreGeometry.h
 *  \brief XXXXXXXXXXXXXX
 *
 *  Copyright 2007-2021 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_SLAB_WITH_CYLINDRICAL_PORE_GEOMETRY_H
#define IMPNPCTRANSPORT_SLAB_WITH_CYLINDRICAL_PORE_GEOMETRY_H

#include "npctransport_config.h"
#include <IMP/display/geometry.h>
#include <IMP/display/display_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

//! XXXX
/** XXXXXX.
 */
class IMPNPCTRANSPORTEXPORT SlabWithCylindricalPoreWireGeometry : public display::Geometry {
  double height_, radius_, length_;

 public:
  //! Slab with specified height from top to bottom, length x length area
  //! and a cylindrical pore of specified radius
  SlabWithCylindricalPoreWireGeometry(double height, double radius, double length);

  //! returns the set of geometric components that comprise this geometry
  //! (for e.g. storing in RMF format)
  virtual IMP::display::Geometries get_components() const IMP_OVERRIDE;

  IMP_OBJECT_METHODS(SlabWithCylindricalPoreWireGeometry);
};

class IMPNPCTRANSPORTEXPORT SlabWithCylindricalPoreSurfaceGeometry
    : public display::SurfaceMeshGeometry {
 public:
  //! Slab with specified height from top to bottom, length x length area
  //! and a cylindrical pore of specified radius
  SlabWithCylindricalPoreSurfaceGeometry(double height, double radius, double length);

  //! returns the set of geometric components that comprise this geometry
  //! (for e.g. storing in RMF format)
  virtual IMP::display::Geometries get_components() const IMP_OVERRIDE;

 private:
  IMP::algebra::Vector3Ds vertices_;
  Ints faces_;
};

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SLAB_WITH_CYLINDRICAL_PORE_GEOMETRY_H */
