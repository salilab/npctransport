/**
 *  \file IMP/npctransport/SlabWithToroidalPore.h
 *  \brief A decorator for a particle that's a slab with a toroidal pore.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPNPCTRANSPORT_SLAB_WITH_TOROIDAL_PORE_H
#define IMPNPCTRANSPORT_SLAB_WITH_TOROIDAL_PORE_H

#include "npctransport_config.h"
#include "SlabWithPore.h"
#include <IMP/Decorator.h>
#include <IMP/decorator_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

//! A decorator for a particle that represents a slab containing
//! a toroidal pore
/** \ingroup helper
    \ingroup decorators
 */
class IMPNPCTRANSPORTEXPORT SlabWithToroidalPore
: public SlabWithPore
{
  /** Decorate a particle that represents a slab (e.g. nuclear
      envelope) with specified thickness and a toroidal pore of
      specified major radius and thickness/2.0 minor radius.
      Note that the radius is controlled by set_pore_radius()
      as any other XYZR particle, but the XYZ coordinates are ignored for now
      (assumed to be 0,0,0).

      The slab is parallel to the x,y plain from z=-0.5*thickness to
      z=0.5*thickness, and the central axis of the pore lies on the
      origin.

      @param m the model
      @param pi the particle index
      @param thickness slab thickness, also twice the minor_radius
      @param major_radius pore major radius
  */
  static void do_setup_particle(IMP::Model* m,
                                ParticleIndex pi,
				double thickness,
				double major_radius);


 public:
  IMP_DECORATOR_METHODS(SlabWithToroidalPore, SlabWithPore);


  /** Decorate a particle that represents a slab (e.g. nuclear
      envelope) with specified thickness and a toroidal pore of
      specified major radius and minor radius of 0.5*thickness.

      The slab is parallel to the x,y plain from z=-0.5*thickness to
      z=0.5*thickness, and the central axis of the pore lies on the
      origin.

      @param m the model
      @param pi the particle index
      @param thickness slab thickness, also twice the minor radius
      @param major_radius pore major radius
  */
  IMP_DECORATOR_SETUP_2(SlabWithToroidalPore,
			double, thickness,
			double, major_radius);

  //! Return true if the particle is an instance of SlabWithToroidalPore
  static bool get_is_setup(Model *m, ParticleIndex pi) {
    return SlabWithPore::get_is_setup(m, pi) &&
      m->get_has_attribute(get_toroidal_pore_key(), pi);
  }

  //! Get the decorator key indicating a toroidal pore
  static IntKey get_toroidal_pore_key();
};



IMP_DECORATORS(SlabWithToroidalPore, SlabsWithToroidalPores, IMP::SlabsWithPores);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SLAB_WITH_TOROIDAL_PORE_H */
