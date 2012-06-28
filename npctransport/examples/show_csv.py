import IMP.npctransport
import math
fname="/home/drussel/kd_output.sorted.csv"

csv= IMP.npctransport.CSV(fname)
#print dir(csv)
variables=["spring_constant", "attraction_range_factor", "number_of_sites"]
output="fraction_bound"

print csv
outdata= csv.get_maximum(variables, output, ["time_step"], [10000], [30000])
print outdata
w= IMP.display.PivyWriter()
def get_coord(c):
    return IMP.algebra.Vector3D(10*math.log(c[0]+.1, 2), 10*math.log(c[1]+.1, 2), c[2])

minp= [x for x in outdata[0]]
maxp= [x for x in outdata[0]]
for d in outdata:
    center= get_coord(d);
    for i in range(0,3):
        #print "comp", i, d, minp, min(d[i], minp[i])
        minp[i]= min(d[i], minp[i])
        maxp[i]= max(d[i], maxp[i])
        #print "result", d, minp, maxp
    g= IMP.display.SphereGeometry(IMP.algebra.Sphere3D(center, .4));
    g.set_color(IMP.display.get_hot_color(d[-1]));
    w.add_geometry(g)

acolors=[IMP.display.Color(1,0,0),
         IMP.display.Color(0,1,0),
         IMP.display.Color(0,0,1)]
for i in range(0,3):
    start= IMP.algebra.get_zero_vector_3d()
    end= IMP.algebra.get_zero_vector_3d()
    start[i]=minp[i]
    end[i]=maxp[i]
    mp= .5*(start+end)
    for j in range(0,3):
        if i==j: continue
        mp[j]= mp[j]=1
    #print start, end, mp
    seg= IMP.algebra.Segment3D(get_coord(start), get_coord(end))
    #print seg
    g= IMP.display.CylinderGeometry(IMP.algebra.Cylinder3D(seg, .2))
    g.set_color(acolors[i])
    w.add_geometry(g)
    l= IMP.display.LabelGeometry(mp, variables[i])
    w.add_geometry(l)
w.show()
