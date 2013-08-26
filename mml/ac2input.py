
import math
import sys

import readac

import datastructures


Atoms,Bonds,Angles,Torsions = readac.main(sys.argv[1])

def main():

    universe = datastructures.Universe()
    
    for atom in Atoms:
        universe.atoms.append(datastructures.Atom(
         symbol=None,atomtype=atom.atomtype,label=atom.label,
         x=atom.x,y=atom.y,z=atom.z,charge=atom.charge))

    for bond in Bonds:
        universe.bonds.append(datastructures.Bond(
         atoms=[universe.atoms[bond.i],universe.atoms[bond.j] ] ) )
           
    for angle in Angles:
        universe.angles.append(datastructures.Angle(
         atoms=[universe.atoms[angle.i],universe.atoms[angle.j],universe.atoms[angle.k] ] ) )

    for tor in Torsions:
        universe.torsions.append(datastructures.Torsion(
         atoms=[universe.atoms[tor.i],universe.atoms[tor.j],
                universe.atoms[tor.k],universe.atoms[tor.l]] ) )

    universe.update_neighbor_assignments()

    return universe 

if __name__ == '__main__':
    main()  
