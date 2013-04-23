/**
 *  \file SlabGeometry.h
 *  \brief XXXXXXXXXXXXXX
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_SLAB_GEOMETRY_H
#define IMPNPCTRANSPORT_SLAB_GEOMETRY_H

#include "npctransport_config.h"
#include <IMP/display/geometry.h>
#include <IMP/display/display_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

//! XXXX
/** XXXXXX.
 */
class IMPNPCTRANSPORTEXPORT SlabWireGeometry: public display::Geometry
{
  double height_, radius_, length_;
public:
  SlabWireGeometry(double height, double radius, double length);
  IMP_GEOMETRY(SlabWireGeometry);
};

class IMPNPCTRANSPORTEXPORT SlabSurfaceGeometry:
  public display::SurfaceMeshGeometry
{
public:
  SlabSurfaceGeometry(double height, double radius, double length);
};


IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_SLAB_GEOMETRY_H */
