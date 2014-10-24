#!/usr/bin/env python

"""Add cylinders of a given radius to the passed file."""
import IMP.em
import sys

scale=4

def prune(mp):
  for xi in range(0, 50):
    for yi in range(0, 100):
      for zi in range(0, 100):
        vi = mp.xyz_ind2voxel(xi, yi, zi)
        mp.set_value(vi, 0)

for m in sys.argv[1:]:
  dm= IMP.em.read_map(m)
  dm.update_voxel_size(scale*dm.get_spacing())
  dm.set_origin(-50*scale, -50*scale, -50*scale)
  nm=m+".new.mrc"
  prune(dm)
  IMP.em.write_map(dm, nm)
  print nm
