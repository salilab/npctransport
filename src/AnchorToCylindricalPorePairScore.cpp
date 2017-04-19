/**
 *  \file AnchorToCylidnricalPorePairScore.cpp
 *  \brief XXXX.
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 *
 */

#include "IMP/npctransport/AnchorToCylindricalPorePairScore.h"
#include <cmath>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

AnchorToCylidnricalPorePairScore
::AnchorToCylidnricalPorePairScore
( Float rot_angle,
  Float radial_d,
  Float z,
  Float k)
  :
  normalized_xy_(cos(rot_angle),
                 sin(rot_angle)),
  pore_radial_d_(radial_d),
  ds_(k)
{
  reference_point_[2]= z; // z is constant
}


AnchorToCylidnricalPorePairScore
::AnchorToCylidnricalPorePairScore
( SlabWithCylindricalPore scp,
  algebra::Vector3D initial_anchor_point,
  Float k )
  :
  ds_(k)
{
  Float x(initial_anchor_point[0]);
  Float y(initial_anchor_point[1]);
  Float r(std::sqrt(x*x+y*y));
  if(r>0.000001){
    normalized_xy_[0]=x/r;
    normalized_xy_[1]=y/r;
  } else {
    // arbitrary values:
    normalized_xy_[0]=1.0;
    normalized_xy_[1]=0.0;
  }
  pore_radial_d_= r - scp.get_pore_radius();
  reference_point_[2]= initial_anchor_point[2]; // same z, always
}

ModelObjectsTemp
AnchorToCylidnricalPorePairScore::do_get_inputs
( Model *m,
  const ParticleIndexes &pis ) const
{
  return IMP::get_particles(m, pis);
}

IMPNPCTRANSPORT_END_NAMESPACE
