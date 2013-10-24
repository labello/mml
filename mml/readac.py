
'''
This module reads and returns basic information about ac files (Amber Antechamber). It 
is fragile and not well tested.  It is only intended to provide basic file
reading capabilities for the purpose of developing the more interesting simulation code.

usage:

import readac
Atoms,Bonds,Angles,Torsions = readac.main(acfile)

'''


import sys
from collections import namedtuple
import collections
import numpy as np

Count = namedtuple('Count','NAtoms NBonds') 
Atom  = namedtuple('Atom','index label x y z atomtype charge')
Bond  = namedtuple('Bond','i j order')
Angle = namedtuple('Angle','index i j k')
Torsion = namedtuple('Torsion','index i j k l')

def main(molfile):
    moldata = open(molfile,'r').readlines()
    moldata = filter(moldata)

    counts,Atoms,Bonds = processlines(moldata)
    check_read(counts.NAtoms,Atoms,counts.NBonds,Bonds)
    angles = get_angles(counts.NBonds,Bonds)
    torsions = get_torsions(counts.NBonds,Bonds,angles)
    
    Angles = []
    index = 0
    for angle in angles:
        Angles.append(Angle(index,angle[0],angle[1],angle[2]))
        index += 1
        
    Torsions = []
    index = 0
    for torsion in torsions:
        Torsions.append(Torsion(index,torsion[0],torsion[1],torsion[2],
         torsion[3]))
         
    return (Atoms,Bonds,Angles,Torsions)   

def get_torsions(NBonds,Bonds,Angles):
    torsions = []
    for angle in Angles:
        ai,ak = angle[0],angle[2]
        aibondlist = get_bonds_atom_participates_in(Bonds,ai)
        akbondlist = get_bonds_atom_participates_in(Bonds,ak)

        for bond in aibondlist:
            if (bond.i not in angle):
                torsions.append([bond.i]+angle)
            if (bond.j not in angle):
                torsions.append([bond.j]+angle)

        for bond in akbondlist:
            if (bond.i not in angle):
                torsions.append(angle+[bond.i])
            if (bond.j not in angle):
                torsions.append(angle+[bond.j])
                               
    outlist = []
    for item in torsions:
        revitem = item[:] ; revitem.reverse()
        if (item and revitem) not in outlist:
            outlist.append(item)
    return outlist

def get_bonds_atom_participates_in(Bonds,index):
    bondlist = []
    for bond in Bonds:
        if index in [bond.i,bond.j]:
            bondlist.append(bond)
    return bondlist
    
def get_angles(NBonds,Bonds):
    angles = []
    for i in range(0,NBonds):
        for j in range(0,NBonds):
            if i != j:
                b1 = collections.Counter([Bonds[i].i,Bonds[i].j])
                b2 = collections.Counter([Bonds[j].i,Bonds[j].j])
# if there is an intersection of the indices in bond1 and bond2
# the index that exists in both bonds will be the center of an angle.
# if there is no intersection the list will return empty
                intersection = list((b1 & b2).elements())
                if intersection:
                    center = intersection[0]
                    b1_remainder = list((b1 - b2).elements())[0]
                    b2_remainder = list((b2 - b1).elements())[0]                   
                    angles.append([b1_remainder,center,b2_remainder])
                    
# (1,0,2) and (2,0,1) are the same angle.  Uniqify the list by removing
# the identical reverse angles.
    outlist = []
    for item in angles:
        revitem = item[:]
        revitem.reverse()
        if (item and revitem) not in outlist:
            outlist.append(item)
    return outlist

def check_read(NAtoms,Atoms,NBonds,Bonds):
    if NAtoms != len(Atoms):
        print "Atom count mismatch!"
        sys.exit()
    if NBonds != len(Bonds):
        print "Bond count mismatch!"
        sys.exit()
            
def processlines(moldata):
    atomlines = []
    bondlines = []
    for l in moldata:
        if 'ATOM' in l:
            atomlines.append(l)
        if 'BOND' in l:
            bondlines.append(l)
            
    counts = Count( len(atomlines),len(bondlines) )
    
    Atoms = []
    index = 0
    for l in atomlines:
        notused1,notused2,label,notused3,notused4,x,y,z,charge,atomtype = l.strip().split() 
        Atoms.append(Atom(index,label,float(x),float(y),float(z),atomtype,float(charge)  )) 
        index += 1
        
    Bonds = []
    for l in bondlines:   
        notused1,notused2,atom1index,atom2index,order,notused3,notused4 = l.strip().split()
        Bonds.append(Bond(int(atom1index)-1,int(atom2index)-1,int(order)))

    return (counts,Atoms,Bonds)           

def filter(moldata):
    x = []
    for line in moldata:
        if 'ATOM' or 'BOND' in line:
            x.append(line)
    return x

if __name__ == '__main__':
    main()
