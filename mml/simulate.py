
import argparse

import numpy as np
import scipy.optimize

import energy
import read_input

def main(molfile):
  
    universe = read_input.main(molfile)

    ffield='gaff'
    energy.energy(ffield,universe,noreport=False)
    cartesians = universe.get_cartesian_atom_array()

    universe.assign_random_velocity(0.1)
    
    print universe.get_cartesian_velocity_array()

    # call the simulation function

#def the simulation function

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='\n\n    Molecular Modeling Lite\n',add_help=True)
    parser.add_argument('--molfile',default='../tests/ethane.mol.ac',
     help='Path to .ac molecule file (Antechamber format)')
    parser.add_argument('--method',default='BFGS',help='Not used right now.')
    args = parser.parse_args()
    main(args.molfile)

