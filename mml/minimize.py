
import argparse

import numpy as np
import scipy.optimize

import energy
import read_input

def main(molfile,ffield,tol):
  
    universe = read_input.main(molfile)

    energy.energy(ffield,universe,noreport=False)

    coordinates = input2cartesian_array(universe.atoms)  
    flat_coords = np.ndarray.flatten(coordinates)

    res = scipy.optimize.minimize(wrap_energy_for_minimizer,flat_coords,args=(ffield,universe),
     method='L-BFGS-B',tol=tol,bounds=None)
     
    energy.energy(ffield,universe,noreport=False)

    print res

    
def wrap_energy_for_minimizer(flat_coords,ffield,universe):

    c = flat_coords.reshape(-1,3)
    # update the universe with the coords passed in
    for i in xrange(0,len(universe.atoms)):
        universe.atoms[i].x = c[i][0]
        universe.atoms[i].y = c[i][1]
        universe.atoms[i].z = c[i][2]

    # write the structure
    write_cartesian_output(universe.atoms,c)

    e = energy.energy(ffield,universe)
    return e

def write_cartesian_output(Atoms,coordinates):
    minfile = open('MINIMIZE.xyz','a')
    minfile.write( str(len(coordinates)) )
    minfile.write( '\n\n' )
    for i in range(0,len(Atoms)):
        minfile.write( "{0:4s} {1:12.5f} {2:12.5f} {3:12.5f} \n".format(
         Atoms[i].atomtype, coordinates[i][0],coordinates[i][1],coordinates[i][2]) )
    
def input2cartesian_array(Atoms):
    return np.vstack(atom.vec for atom in Atoms )
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='\n\n    Molecular Modeling Lite\n',add_help=True)
    parser.add_argument('--molfile',default='/Users/labello/gitwork/mml/no-git/molecules/acfiles/ethane.mol.ac',
     help='Path to .ac molecule file (Antechamber format)')
    parser.add_argument('--ffield',default='gaff',help='amber94 or gaff')
    parser.add_argument('--tol',default=0.01,help='Energy tolerance for geometry optimization')
    parser.add_argument('--method',default='BFGS',help='Minimization routine (BFGS, Newton-CG, TNC, SLSSQP)')
    args = parser.parse_args()
    main( args.molfile,args.ffield,float(args.tol))

