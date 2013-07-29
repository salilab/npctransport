/**
 *  \file npctransport/Parameter.h
 *  \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_PARAMETER_H
#define IMPNPCTRANSPORT_PARAMETER_H

#include "npctransport_config.h"
#include <IMP/base/check_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

#ifndef SWIG
template <class T>
struct Parameter {
  T t_;
  bool init_;

public:
Parameter() : init_(false) {}
  operator T() const {
    IMP_USAGE_CHECK(init_, "Not initialized");
    return t_;
  }
  void operator=(T t) {
    t_ = t;
    init_ = true;
  }
    bool is_init() { return init_; }
};

#endif

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_PARAMETER_H */
