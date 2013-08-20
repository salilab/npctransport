import IMP.core
import IMP.em
import IMP.atom
import IMP.display
import IMP.kmeans
import IMP.rmf
import RMF
import os.path
import re
import math
#from optparse import OptionParser
import sys

def translate_particle(p, translation ):
    pXYZ = IMP.core.XYZ( p )
    new_coords = pXYZ.get_coordinates() + translation
    pXYZ.set_coordinates( new_coords )

def copy_and_translate_XYZR_particle(p, translation ):
    new_p = IMP.Particle( p.get_model() )
    new_coords = IMP.core.XYZ( p ).get_coordinates() + translation
    new_radius = IMP.core.XYZR( p ).get_radius()
    new_sphere = IMP.algebra.Sphere3D( new_coords, new_radius )
    IMP.core.XYZR.setup_particle( new_p, new_sphere )
    return new_p

def show_hierarchy( h, prefix="" ):
    print prefix,
    h.show()
    print
    for child in h.get_children():
        show_hierarchy( child, prefix + "> " )

def cluster_MRC_file_with_cache(mrc_filename, k):
    cache = mrc_filename + ".cache"
    centers = []
    mean_loc = None
    ok = False
    try: # get from cache
        CACHE=open(cache,"r")
        begin=False
        mrc=False
        k_ok=False
        should_end=False
        i=0
        for line in CACHE:
            line=line.strip()
            if(line=="BEGIN"):
                begin=True
                continue
            if(line==mrc_filename):
                mrc=True
                continue
            if(mrc):
                mrc=False #reset
                if k <> int(line):
                    print "Wrong k", int(line)
                    break
                k_ok=True
                continue
            if(k_ok and i<k):
                center = [float(x) for x in line.split()]
                if len(center)<>3:
                    print "Invalid center length", center
                    break
                centers.append(center)
                i=i+1
                continue
            if(k_ok and i==k and not should_end):
                mean_loc = [float(x) for x in line.split()]
                if len(mean_loc)<>3:
                    print "Invalid mean loc length", mean_loc
                    break
                should_end = True
                continue
            if(should_end and line=="END"):
                ok=True
                break
    except:
        pass
    if(ok):
        print "Read", k, "centers from cache", cache
    if(not ok):
        print "reading", k, "centers from file", mrc_filename
        kmeans, centers, mean_loc, pos_voxels \
            = cluster_MRC_file(mrc_filename, k)
        CACHE = open(cache, "w")
        print >>CACHE, "BEGIN"
        print >>CACHE, mrc_filename
        print >>CACHE, k
        for pos in centers+[mean_loc]:
            for x in pos:
                print >>CACHE, x,
            print >>CACHE
        print >>CACHE, "END"
        CACHE.close()
    assert(len(centers)==k and len(mean_loc)==3)
    return centers, mean_loc



def cluster_MRC_file(input_fname,
                     K, # number of clusters
                     resolution=8.0,
                     voxel_size=1.5,
                     INTENSE_PIXEL_THRESH=0.1): # the intensity of a voxel to be defined as intense
    print "Reading MRC file ", input_fname
    if not os.path.exists(input_fname):
        raise IOError("File %s not found" % input_fname)
    dmap=IMP.em.read_map(input_fname + "", IMP.em.MRCReaderWriter())
    dmap.get_header_writable().set_resolution(resolution)
    nvoxels = dmap.get_number_of_voxels()
#    print "Bounding box: %s" % IMP.em.get_bounding_box(dmap, 0.001)
#    print "Number of voxels: %d" % nvoxels
    #  dmap_resample = dmap.
    sum_location = IMP.algebra.Vector3D(0,0,0)
    sum_value = 0.0
    max_value = 0.0
    intense_voxel = None # a voxel that is intense enough to be taken as an arbitrary sample
    pos_voxels_locations = []
    for i in range(nvoxels):
        value = dmap.get_value(i)
        location = dmap.get_location_by_voxel(i) * 10 # convert from nanometers to angstroms
        if(value > 0):
            sum_location = sum_location + location * value
            sum_value = sum_value + value
            if value > INTENSE_PIXEL_THRESH:
                pos_voxels_locations.append( (location, value) )
#                print "Voxel %d location %s ; value %f" % (i, location, value)
            # save all intense pixels in a list
            if(value > INTENSE_PIXEL_THRESH and intense_voxel is None):
                intense_voxel = i
#                print "Intense Voxel %d location %s ; value %f" % (i, location, value)
            max_value = max( value, max_value )
    if(len(pos_voxels_locations) == 0):
#        Print "ERROR: No intense voxels in file '" + input_fname + "' with max value", max_value,
#        print "Sum value: ", sum_value
        raise IMP.ValueException("No intense voxels in file %s" % input_fname)
    mean_location = sum_location / sum_value
#    print "Mean location: %s" % mean_location, "Sum value: %.2f" % sum_value
    intense_voxel_location = dmap.get_location_by_voxel(intense_voxel) * 10.0 # from nm to angstroms
#    print "Intense voxel location: %s" % intense_voxel_location
    # generate some clusters of points using kmeans
    ncycles = 5000
#    for i,voxel in enumerate(pos_voxels_locations):
#        print i, voxel
    kmeans = IMP.kmeans.KMeans()
    for voxel in (pos_voxels_locations):
        kmeans.add_data_pt( voxel[0] )
    kmeans.execute(K, IMP.kmeans.KM_LLOYDS, ncycles)
    centers = []
    for i in range(K):
        centers.append( kmeans.get_center(i) )
#        print "%d %s" % (i , centers[i] )
    return (kmeans, centers, mean_location, pos_voxels_locations)



#################### MAIN ################
def main():
    print "DEGENERATED"

if __name__ == "__main__":
    main()
