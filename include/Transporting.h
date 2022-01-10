/**
 *  \file IMP/npctransport/Transporting.h
 *  \brief A decorator for a transporting particle, so as to keep track of its
* directionality.
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
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
class IMPNPCTRANSPORTEXPORT Transporting : public IMP::Decorator {
  /** Decorate a transporting particle, mainly for tracking transport
      statistics (number of transports, etc.).  It is assumed the
      transport occurs along a z coordinate, and that the crossed barrier
      is bounded from top and bottom.
      The last tracked z is set to the z coordinate of p, initially,
      and n entries from bottom and top are set to 0
      @see get_last_tracked_z
      @see get_n_entries_bottom
      @see get_n_entries_top

      @param m the model
      @param pi the particle index
      @param is_last_entry_from_top has particle last entered from top of
                                    barrier (rather than bottom or unknown)
  */
  static void do_setup_particle(IMP::Model* m,
                                ParticleIndex pi,
                                bool is_last_entry_from_top = false);

 public:
  IMP_DECORATOR_METHODS(Transporting, Decorator);


  /** Decorate a transporting particle, mainly for tracking transport
      statistics (number of transports, etc.).  It is assumed the
      transport occurs along a z coordinate, and that the crossed barrier
      is bounded from top and bottom.
      The last tracked z is set to the z coordinate of p, initially,
      and n entries from bottom and top are set to 0
      @see get_last_tracked_z
      @see get_n_entries_bottom
      @see get_n_entries_top

      @param m the model
      @param pi particle index
      @param is_last_entry_from_top has particle last entered from top of
     barrier
                                    (rather than bottom or unknown)
  */
  IMP_DECORATOR_SETUP_1(Transporting,  bool, is_last_entry_from_top = false);

  //! Return true if the particle is an instance of an Transporting
  static bool get_is_setup(Model *m, ParticleIndex pi) {
    return m->get_has_attribute(get_is_last_entry_from_top_key(), pi) &&
           m->get_has_attribute(get_last_tracked_z_key(), pi) &&
           m->get_has_attribute(get_n_entries_bottom_key(), pi) &&
           m->get_has_attribute(get_n_entries_top_key(), pi);
  }

  //! sets whether the particle last enetered the transport moiety from its top
  void set_is_last_entry_from_top(bool is_last_entry_from_top) {
    get_particle()->set_value(get_is_last_entry_from_top_key(),
                              is_last_entry_from_top ? 1 : 0);
  }

  //! returns whether the particle last enetered the transport moiety from its
  //top
  bool get_is_last_entry_from_top() const {
    return (get_particle()->get_value(get_is_last_entry_from_top_key()) != 0);
  }

  //! Get the decorator key for is_last_entry_from_top
  static IntKey get_is_last_entry_from_top_key();

  //! set the Z coordinate of the particle, the last
  //! time it was tracked for transport statistics
  void set_last_tracked_z(double last_tracked_z) {
    get_particle()->set_value(get_last_tracked_z_key(), last_tracked_z);
  }

  //! returns the Z coordinate of the particle, the last
  //! time it was tracked for transport statistics
  double get_last_tracked_z() const {
    return get_particle()->get_value(get_last_tracked_z_key());
  }

  //! Get the decorator key last_tracked_z value
  static FloatKey get_last_tracked_z_key();

  //! sets the number of times the particle crossed from the bottom
  //! of the barrier into its interior
  void set_n_entries_bottom(int n) {
    get_particle()->set_value(get_n_entries_bottom_key(), n);
  }

  //! gets the number of times the particle crossed from the bottom
  //! of the barrier into its interior
  int get_n_entries_bottom() const {
    return get_particle()->get_value(get_n_entries_bottom_key());
  }

  //! Get the decorator key n_entries_bottom value
  static IntKey get_n_entries_bottom_key();

  //! sets the number of times the particle crossed from the top
  //! of the barrier into its interior
  void set_n_entries_top(int n) {
    get_particle()->set_value(get_n_entries_top_key(), n);
  }

  //! gets the number of times the particle crossed from the top
  //! of the barrier into its interior
  int get_n_entries_top() const {
    return get_particle()->get_value(get_n_entries_top_key());
  }

  //! Get the decorator key n_entries_top value
  static IntKey get_n_entries_top_key();
};



IMP_DECORATORS(Transporting, Transportings, IMP::Decorators);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_TRANSPORTING_H */
