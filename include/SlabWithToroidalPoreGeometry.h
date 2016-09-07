/**
 *  \file SlabWithToroidalPoreGeometry.h
 *  \brief XXXXXXXXXXXXXX
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_SLAB_WITH_TOROIDAL_PORE_GEOMETRY_H
#define IMPNPCTRANSPORT_SLAB_WITH_TOROIDAL_PORE_GEOMETRY_H

#include "npctransport_config.h"
#include <IMP/display/geometry.h>
#include <IMP/display/display_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

//! XXXX
/** XXXXXX.
 */
class IMPNPCTRANSPORTEXPORT SlabWithToroidalPoreWireGeometry : public display::Geometry {
  double r_; // minor radius
  double R_; // major radius
  double slab_length_; // length of slab edge

 public:
  //! Slab with specified height from top to bottom, length x length area
  //! and a toroidal pore of specified major radius and slab_height/2.0 minor radius
  SlabWithToroidalPoreWireGeometry(double slab_height, double major_radius, double slab_length);

  //! returns the set of geometric components that comprise this geometry
  //! (for e.g. storing in RMF format)
  virtual IMP::display::Geometries get_components() const IMP_OVERRIDE;

  IMP_OBJECT_METHODS(SlabWithToroidalPoreWireGeometry);
};


IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SLAB_WITH_TOROIDAL_PORE_GEOMETRY_H */
