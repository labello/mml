
import argparse

import numpy as np
import scipy.optimize

import energy
import read_input

def main(molfile,method,tol):
  
    universe = read_input.main(molfile)

    ffield='gaff'
    energy.energy(ffield,universe,noreport=False)

    flat_coords = universe.get_flat_cartesian_atom_array()

    res = scipy.optimize.minimize(wrap_energy_for_minimizer,flat_coords,args=(ffield,universe),
     method=method,tol=tol,bounds=None)
    
    energy.energy(ffield,universe,noreport=False)

    print res

    
def wrap_energy_for_minimizer(flat_coords,ffield,universe):

    c = flat_coords.reshape(-1,3)
    for i in xrange(0,len(universe.atoms)):
        universe.atoms[i].update_coords(c[i][0],c[i][1],c[i][2])

    universe.write_xyz()

    e = energy.energy(ffield,universe)
    return e

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='\n\n    Molecular Modeling Lite\n',add_help=True)
    parser.add_argument('--molfile',default='/Users/labello/gitwork/mml/no-git/molecules/acfiles/ethane.mol.ac',
     help='Path to .ac molecule file (Antechamber format)')
    parser.add_argument('--tol',default=0.01,help='Energy tolerance for geometry optimization')
    parser.add_argument('--method',default='BFGS',help='Minimization routine (BFGS, L-BFGS-B,Anneal,Newton-CG, TNC, SLSSQP)')
    args = parser.parse_args()
    main( args.molfile,args.method,float(args.tol))

