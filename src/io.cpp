/**
 *  \file io.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/io.h>
#include <IMP/npctransport/SitesGeometry.h>
#include <IMP/npctransport/SimulationData.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/base/Pointer.h>
#include <IMP/display/restraint_geometry.h>
#include <RMF/decorators.h>
#include <IMP/rmf/atom_io.h>
#include <IMP/rmf/geometry_io.h>
#include <IMP/rmf/restraint_io.h>
#include <IMP/rmf/frames.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

void write_geometry(const ParticlesTemp &kaps,
                    const algebra::Vector3Ds &kap_sites,
                    const ParticlesTemp &chains,
                    const algebra::Vector3Ds &chain_sites,
                    const RestraintsTemp &rs, display::Writer *out) {
  for (unsigned int i = 0; i < kaps.size(); ++i) {
    IMP_NEW(SitesGeometry, g, (kaps[i], kap_sites));
    out->add_geometry(g);
  }
  for (unsigned int i = 0; i < chains.size(); ++i) {
    ParticlesTemp leaves =
        get_as<ParticlesTemp>(atom::get_leaves(atom::Hierarchy(chains[i])));
    for (unsigned int j = 0; j < leaves.size(); ++j) {
      IMP_NEW(SitesGeometry, g, (leaves[j], chain_sites));
      out->add_geometry(g);
    }
  }
  for (unsigned int i = 0; i < rs.size(); ++i) {
    IMP_NEW(display::RestraintGeometry, g, (rs[i]));
    out->add_geometry(g);
  }
}
IMPNPCTRANSPORT_END_NAMESPACE
