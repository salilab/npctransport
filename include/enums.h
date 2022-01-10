/**
# * \file npctransport/enums.h
# * \brief useful enums and contants

#
# * Copyright 2007-2022 IMP Inventors. All rights reserved.
# */

#ifndef IMPNPCTRANSPORT_ENUMS_H
#define IMPNPCTRANSPORT_ENUMS_H

#include "npctransport_config.h"
#include <limits>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

const double FS_IN_NS = 1.0E+6;

const double MAX_DOUBLE = std::numeric_limits< double >::max();

const double HALF_SQRT_MAX_DOUBLE = 0.5 * std::sqrt( MAX_DOUBLE );

IMPNPCTRANSPORT_END_NAMESPACE


#endif /* IMPNPCTRANSPORT_ENUMS_H */
