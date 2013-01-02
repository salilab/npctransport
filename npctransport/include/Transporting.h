/**
 *  \file IMP/npctransport/Transporting.h
 *  \brief A decorator for a transporting particle, so as to keep track of its directionality.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPNPCTRANSPORT_TRANSPORTING_H
#define IMPNPCTRANSPORT_TRANSPORTING_H

#include "npctransport_config.h"
#include <IMP/Decorator.h>
#include <IMP/decorator_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

//! A decorator for a particle transporting through a barrier.
/** \ingroup helper
    \ingroup decorators
 */
class IMPNPCTRANSPORTEXPORT Transporting:
  public IMP::Decorator
{

 public:
  IMP_DECORATOR(Transporting, IMP::Decorator);

 public:
  /** Create a decorator with the passed coordinates and D.
  */
  static Transporting setup_particle(Particle *p,
                                  bool is_last_entry_from_top) {
    p->add_attribute(get_is_last_entry_from_top_key(), is_last_entry_from_top);
    return Transporting(p);
  }

  //! Return true if the particle is an instance of an Transporting
  static bool particle_is_instance(Particle *p) {
    return p->has_attribute(get_is_last_entry_from_top_key());
  }

  //! Return true if the particle is an instance of an Transporting
  static bool particle_is_instance(Model *m, ParticleIndex p) {
    return m->get_has_attribute(get_is_last_entry_from_top_key(), p);
  }

  //! sets whether the particle last enetered the transport moiety from its top
  void set_is_last_entry_from_top(bool is_last_entry_from_top) {
    get_particle()->set_value(get_is_last_entry_from_top_key(),
                              is_last_entry_from_top ? 1 : 0);
  }

  //! returns whether the particle last enetered the transport moiety from its top
  bool get_is_last_entry_from_top() const {
    return (get_particle()->get_value(get_is_last_entry_from_top_key()) != 0);
  }
  //! Get the decorator key for is_last_entry_from_top
  static IntKey get_is_last_entry_from_top_key();
};


IMP_DECORATORS(Transporting, Transportings, IMP::Decorators);

IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_TRANSPORTING_H */
