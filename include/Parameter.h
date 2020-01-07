/**
 *  \file npctransport/Parameter.h
 *  \brief description
 *
 *  Copyright 2007-2020 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_PARAMETER_H
#define IMPNPCTRANSPORT_PARAMETER_H

#include "npctransport_config.h"
#include <IMP/check_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

template <class T>
class Parameter {
 private:
  T t_;
  bool init_;

 public:
 Parameter() : init_(false) {}

 Parameter(T t) : t_(t), init_(true) {}

  T get_value() const {
  IMP_USAGE_CHECK(init_, "npctransort::Parameter Not initialized");
  return t_;
  }

#ifndef SWIG
  operator T() const {
  return get_value();
}

  void operator=(T t) {
  t_ = t;
  init_ = true;
}
#endif

  bool is_init() { return init_; }
};


IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_PARAMETER_H */
