/**
 *  \file IMP/npctransport/SlabWithCylindricalPore.h
 *  \brief A decorator for a particle that's a slab with a cylindrical pore.
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPNPCTRANSPORT_SLAB_WITH_CYLINDRICAL_PORE_H
#define IMPNPCTRANSPORT_SLAB_WITH_CYLINDRICAL_PORE_H

#include "npctransport_config.h"
#include "SlabWithPore.h"
#include <IMP/Decorator.h>
#include <IMP/decorator_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

//! A decorator for a particle that represents a slab containing
//! a cylindrical pore
/** \ingroup helper
    \ingroup decorators
 */
class IMPNPCTRANSPORTEXPORT SlabWithCylindricalPore
: public SlabWithPore
{
  /** Decorate a particle that represents a slab (e.g. nuclear
      envelope) with specified thickness and a cylindrical pore of
      specified radius. Note that the radius is controlled by set_radius()
      as any other XYZR particle, but the XYZ coordinates are ignored for now
      (assumed to be 0,0,0).

      The slab is parallel to the x,y plain from z=-0.5*thickness to
      z=0.5*thickness, and the central axis of the pore lies on the
      origin.

      @param m the model
      @param pi the particle index
      @param thickness slab thickness
      @param radius pore radius
  */
  static void do_setup_particle(IMP::Model* m,
                                ParticleIndex pi,
				double thickness,
				double radius);


 public:
  IMP_DECORATOR_METHODS(SlabWithCylindricalPore, SlabWithPore);


  /** Decorate a particle that represents a slab (e.g. nuclear
      envelope) with specified thickness and a cylindrical pore of
      specified radius.

      The slab is parallel to the x,y plain from z=-0.5*thickness to
      z=0.5*thickness, and the central axis of the pore lies on the
      origin.

      @param m the model
      @param pi the particle index
      @param thickness slab thickness
      @param radius pore radius
  */
  IMP_DECORATOR_SETUP_2(SlabWithCylindricalPore,
			double, thickness,
			double, radius);

  //! Return true if the particle is an instance of SlabWithCylindricalPore
  static bool get_is_setup(Model *m, ParticleIndex pi) {
    return SlabWithPore::get_is_setup(m, pi) &&
      m->get_has_attribute(get_cylindrical_pore_key(), pi);
  }

  //! Get the decorator key indicating a cylindrical pore
  static IntKey get_cylindrical_pore_key();
};



IMP_DECORATORS(SlabWithCylindricalPore, SlabsWithCylindricalPores, IMP::SlabsWithPores);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SLAB_WITH_CYLINDRICAL_PORE_H */
