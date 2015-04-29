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
#import argparse
import sys

def get_options():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-input', metavar='mrcfile', type=str, nargs='+',
                        help='input MRC file(s)')
    parser.add_argument('-output', metavar='rmffile', type=str, nargs=1,
                        help='output RMF file')
    options = parser.parse_args(  )
    print "Parsed options: ", options
    return options

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
    if(not ok):
        print "No cache found - reading", k, "centers from file", mrc_filename
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
    else:
        pass
#        print "Read", k, "centers from cache", cache
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
    options = get_options()
    # Root of hierarchy:
    model = IMP.Model()
    p_root = IMP.Particle( model, "npc" )
    IMP.core.XYZ.setup_particle( p_root, [0,0,0] )
    IMP.core.XYZR.setup_particle( p_root, 50.0 )
    h_root = IMP.atom.Hierarchy.setup_particle( p_root )
    n_inputs = len(options.input)

    # Add all inputs to hierarchy
    if(n_inputs == 1):
        color_gap = 0.5
    else:
        color_gap = 1.0 / (n_inputs - 1)
    color_value = 0.0
    for input_fname in options.input:
        re_search=re.search('([0-9]+)copies',input_fname)
        K=int(re_search.group(1))
#        print "K=%d" % K
        #  for K in [8,16,24,32]: # number of clusters
#        print "COLOR VALUE ", color_value
        try:
            (kmeans,centers,mean_location,pos_voxels_locations)=cluster_MRC_file(input_fname,K)
        except IMP.UsageException:
            print "ERROR: Couldn't open file '" + input_fname + "'"

        # RMF

        # # Intense Samples:
        # p_intense_samples = IMP.Particle( model, "intense samples " + input_fname )
        # IMP.core.XYZ.setup_particle( p_intense_samples, mean_location )
        # IMP.core.XYZR.setup_particle( p_intense_samples, 5.0 )
        # h_intense_samples = IMP.atom.Hierarchy.setup_particle( p_intense_samples )
        # h_root.add_child( h_intense_samples )
        # for p in particles:
        #     hp = IMP.atom.Hierarchy.setup_particle( p )
        #     h_intense_samples.add_child(hp)
        # Centers:
        p_centers = IMP.Particle( model, "Centers " + input_fname )
        IMP.core.XYZ.setup_particle( p_centers, mean_location )
        IMP.core.XYZR.setup_particle( p_centers, 30.0 )
        h_centers = IMP.atom.Hierarchy.setup_particle( p_centers )
        h_root.add_child( h_centers )
        centers_color = IMP.display.get_hot_color( color_value )
        scaling_factor = 1.0
        for (k, center) in enumerate(centers):
            delta = [a-b for a,b in zip(center, mean_location)]
            scaled_pos = [m+d*scaling_factor for m,d in zip(mean_location,delta)]
            scaled_radius = 30.0 * scaling_factor
            label = "center %d" % k
            p = IMP.Particle( model, label )
            IMP.core.XYZR.setup_particle( p, IMP.algebra.Sphere3D( scaled_pos, scaled_radius ) )
            IMP.display.Colored.setup_particle( p, centers_color )
            hp = IMP.atom.Hierarchy.setup_particle( p )
            h_centers.add_child( hp )
        # Voxels:
        assignments = kmeans.get_assignments()
        SSD=0.0 # Sum Square Deviation from centers
        n=0.0
        # create a hierarchy with specific voxels for each center (which can be omitted from RMF)
        for (voxel, k) in zip(pos_voxels_locations, assignments):
            delta = [a-b for a,b in zip(voxel[0], mean_location)]
            scaled_pos = [m+d*scaling_factor for m,d in zip(mean_location,delta)]
            scaled_radius = voxel[1] * 10.0 * scaling_factor # use weight for radius
            p = IMP.Particle( model, "voxel C%d" % k )
            p_sphere = IMP.algebra.Sphere3D( scaled_pos, scaled_radius )
            IMP.core.XYZR.setup_particle( p, p_sphere )
            IMP.display.Colored.setup_particle( p, centers_color )
            hp = IMP.atom.Hierarchy.setup_particle( p )
            #        h_centers.get_child( k ).add_child( hp ) # if uncommented, voxels are added to hierarchy
            SSD = SSD + voxel[1]*sum([(a-b)**2 for a,b in zip(voxel[0],centers[k])])
            n=n+voxel[1]
#        print "fname=%s, k=%d, Root Mean Square Deviation = %.4f" % (input_fname,K,math.sqrt(SSD/n))
        color_value = color_value + color_gap

#        print "Printing hierarchy"
#        show_hierarchy(h_root)

        rmf= RMF.create_rmf_file(options.output[0])
        # IMP.rmf.add_hierarchy(rmf, h_intense_samples)
        IMP.rmf.add_hierarchy(rmf, h_root)
        IMP.rmf.save_frame(rmf)
        del rmf


# # display symmetry and centers
# wchimera = IMP.display.ChimeraWriter("out_chimera.py")
# wpymol = IMP.display.PymolWriter("out_pymol.py")
# for (i,p) in enumerate(particles):
#     print p.get_name(), IMP.core.XYZ(p).get_coordinates()
#     g = IMP.core.XYZRGeometry( p )
#     g.set_name( "p" + str(i) )
#     g.set_color( IMP.display.Colored(p).get_color() )
#     wchimera.add_geometry( g )
#     wpymol.add_geometry( g )
# for (i,center) in enumerate(centers):
#     print "Center " + str(i), center
#     sphere = IMP.algebra.Sphere3D( center, radius )
#     g = IMP.display.SphereGeometry( sphere )
#     g.set_name( "center" + str(i) )
#     g.set_color( IMP.display.Color(0,1,0) )
#     wchimera.add_geometry( g )
#     wpymol.add_geometry( g )
# del wchimera  # make sure that the file is flushed
# del wpymol
# # RMF


if __name__ == "__main__":
    main()
