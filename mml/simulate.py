
import argparse

import numpy as np
import scipy.optimize

import energy
import read_input

SECONDS_TO_FEMPTOSECONDS = 10**15
FEMPTOSECONDS_TO_SECONDS = 1.0000 / SECONDS_TO_FEMPTOSECONDS 

def main(molfile):
  
    ffield = 'gaff' # make command line option
    ts = 0.1 # make command line option, in femptoseconds
    universe = read_input.main(molfile)
    update_universe_mass(universe,ffield)
    energy.energy(ffield,universe,noreport=False)
    cartesians = universe.get_cartesian_atom_array()
    universe.assign_random_velocity(100)
    
    print "Temperature: ", universe.temperature()

# next step maybe : convert temperature and velocity functions to operate in 
# femptoseconds.  The expand calc_euler_step to actually calcualte the new 
# positions of atoms based on current velocity
    calc_euler_step(universe,deltat)
    # call the simulation function

#def the simulation function

def calc_euler_step(universe,deltat):
    '''
    Document stuff here
    '''
    vels = universe.get_cartesian_velocity_array()
    print vels


    
#    return new_positions


def update_universe_mass(universe,ffield):

    if ffield == 'gaff':
        from amberff import gaffparams as ffparams
        for atom in universe.atoms:
            atom.mass = ffparams.Atoms[atom.atomtype].mass 


    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='\n\n    Molecular Modeling Lite\n',add_help=True)
    parser.add_argument('--molfile',default='../tests/ethane.mol.ac',
     help='Path to .ac molecule file (Antechamber format)')
    parser.add_argument('--method',default='BFGS',help='Not used right now.')
    args = parser.parse_args()
    main(args.molfile)

