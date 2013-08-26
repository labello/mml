
import argparse
import importlib

import numpy as np
import scipy.optimize

import energy

def main(molecule,ffield,tol):
  
    importmodule = importlib.import_module( 'molecules.' + molecule )
    system = importmodule.inputdata
    coordinates = input2cartesian_array(system.Atoms)  
    flat_coords = np.ndarray.flatten(coordinates)
    energy.energy(coordinates,ffield,system,noreport=True)

    res = scipy.optimize.minimize(wrap_energy_for_minimizer,flat_coords,args=(ffield,system),
     method='BFGS',tol=tol,bounds=None)
     
    print res

    
def wrap_energy_for_minimizer(flat_coords,ffield,system):
    c = flat_coords.reshape(-1,3)
    write_cartesian_output(system.Atoms,c)
    e = energy.energy(c,ffield,system)
    return e

def write_cartesian_output(Atoms,coordinates):
    minfile = open('MINIMIZE.xyz','a')
    minfile.write( str(len(coordinates)) )
    minfile.write( '\n\n' )
    for i in range(0,len(Atoms)):
        minfile.write( "{0:4s} {1:12.5f} {2:12.5f} {3:12.5f} \n".format(
         Atoms[i].type, coordinates[i][0],coordinates[i][1],coordinates[i][2]) )
    
def input2cartesian_array(Atoms):
    return np.vstack(atom.vec for atom in Atoms )
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='\n\n    Molecular Modeling Lite\n',add_help=True)
    parser.add_argument('--molecule',default='ethane',help='The molecular system from ../mml/molecules')
    parser.add_argument('--ffield',default='gaff',help='amber94 or gaff')
    parser.add_argument('--tol',default=5,help='Energy tolerance for geometry optimization')
    parser.add_argument('--method',default='BFGS',help='Minimization routine (BFGS, Newton-CG, TNC, SLSSQP)')
    args = parser.parse_args()
    main( args.molecule,args.ffield,float(args.tol))

