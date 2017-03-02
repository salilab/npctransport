import glob
import re

def get_xyz(atom_line):
    assert(re.match('^ATOM',atom_line) and len(atom_line)>=54)
    x=float(atom_line[30:38])
    y=float(atom_line[38:46])
    z=float(atom_line[46:54])
    return [x,y,z]


def get_mean_atom_line(lines):
    X=[]
    Y=[]
    Z=[]
    for i,line in enumerate(lines):
        x,y,z= get_xyz(line)
        X.append(x)
        Y.append(y)
        Z.append(z)
    x=sum(X)/len(X)
    y=sum(Y)/len(Y)
    z=sum(Z)/len(Z)
    mean_atom_line= "%s%8.3f%8.3f%8.3f%s" % (line[0:30], x, y, z, line[54:])
    return mean_atom_line


N=3
PDBs=[]
for i,pdb_fname in enumerate(glob.glob('N2movie16_dump5*.pdb')):
    F=open(pdb_fname,'r')
    PDBs.append(F.readlines())
    n_lines=len(PDBs[-1])

print "Smoothing"

for i in range(0,len(PDBs)-N):
    fname= 'smooth_frame%04d.pdb' % i
    F= open(fname,'w')
    for j in range(n_lines):
        lines=[]
        for k in range(i,i+N):
            lines.append(PDBs[k][j])
        if re.match('^ATOM',lines[0]):
            print >>F, get_mean_atom_line(lines),
        else:
            for line in lines[1:]:
                assert(lines[0]==line)
            print >>F, lines[0],
    F.close()
