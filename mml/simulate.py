
import argparse
import copy

import numpy as np
import scipy.optimize

import constants
import energy
import read_input

ffield = 'gaff' # make command line option

def main(molfile):
  
    ts = 0.015 # make command line option, in picoseconds 
    nsteps = 10**3 # make command line option

    universe = read_input.main(molfile)
    update_universe_mass(universe,ffield)
 #  energy.energy(ffield,universe,noreport=False)
    energy.energy(ffield,universe,noreport=True)

    cartesians = universe.get_cartesian_atom_array()
    universe.assign_random_velocity(100 * 10**-2)  # angstroms/picosecond  (400 * 10**-2)
    
#####################
    simulate(universe,deltat=ts,steps=nsteps)

def simulate(universe,deltat,steps):
    '''
    Document the simulation steps...
    '''

    newuniverse = take_euler_step(universe,deltat) 

    print newuniverse.get_cartesian_atom_array()

    step = 0
    while step < steps:
        calc_forces(newuniverse)
        calc_accelerations(newuniverse)
# verlet definately does not work
        universe,newuniverse = take_verlet_step(universe,newuniverse,deltat)
        print newuniverse.get_cartesian_atom_array()
###### work here ^^^
        newuniverse.write_xyz()
        step += 1
         

def take_verlet_step(prevuniverse,universe,deltat):
    '''
    Verlet algorithm
    '''
    newuniverse = copy.deepcopy(universe)

    for i in xrange(0,newuniverse.natoms()):
        x = universe.atoms[i].x
        y = universe.atoms[i].y
        z = universe.atoms[i].z
        prevx = prevuniverse.atoms[i].x
        prevy = prevuniverse.atoms[i].y
        prevz = prevuniverse.atoms[i].z
        ax = newuniverse.atoms[i].accelerationvec[0] 
        ay = newuniverse.atoms[i].accelerationvec[1] 
        az = newuniverse.atoms[i].accelerationvec[2] 
#        print ax,ay,az 
#        print "AX",ax,ay,az,deltat,deltat**2
        newuniverse.atoms[i].update_coords( 
         ( 2*x - prevx + ax*deltat**2 ),( 2*y - prevy + ay*deltat**2 ),( 2*z - prevz + az*deltat**2 ) ) 

    return (universe,newuniverse)


def calc_forces(universe,stepsize=0.0001):
    '''
     f(x+h) - f(x-h)
     ---------------
          2h
    '''

    for i in xrange(0,universe.natoms()): 

        x,y,z = universe.atoms[i].x,universe.atoms[i].y,universe.atoms[i].z
        p1 = copy.deepcopy(universe)
        
        p1.atoms[i].update_coords( x-stepsize, y,          z)          # minus x universe 
        mxe = energy.energy(ffield,p1)
        p1.atoms[i].update_coords( x+stepsize, y,          z)          # plus  x universe 
        pxe = energy.energy(ffield,p1)
        p1.atoms[i].update_coords( x         , y-stepsize, z)          # minus y universe 
        mye = energy.energy(ffield,p1)
        p1.atoms[i].update_coords( x         , y+stepsize, z)          # plus  y universe 
        pye = energy.energy(ffield,p1)
        p1.atoms[i].update_coords( x         , y,          z-stepsize) # minus z universe 
        mze = energy.energy(ffield,p1)
        p1.atoms[i].update_coords( x         , y,          z+stepsize) # plus  z universe 
        pze = energy.energy(ffield,p1)

        universe.atoms[i].forcevec = -1 * np.array(
         [(pxe-mxe)/(2.0*stepsize), (pye-mye)/(2.0*stepsize), (pze-mze)/(2.0*stepsize)] )
#  forces are in kcal/mol/angstrom
        print "p1: ",mxe,pxe,mye,pye,mze,pze
        print "force vec: ",universe.atoms[i].forcevec 
                

def calc_accelerations(universe):
    '''
    a = F/m
    '''
    for atom in universe.atoms:
        atom.accelerationvec = atom.forcevec / atom.mass
    # acceleration is in kcal/mol/angstrom/amu


def take_euler_step(universe,deltat):
    '''
    positions in angstroms
    velocities in angstroms/picosecond
    deltat in picoseconds 
    '''
    positions = universe.get_cartesian_atom_array()
    velocities = universe.get_cartesian_velocity_array()
    newuniverse = copy.deepcopy(universe)
    for atom in newuniverse.atoms:
        atom.update_coords( (atom.x + atom.vx*deltat), (atom.y + atom.vy*deltat), (atom.z + atom.vz*deltat) )
    return newuniverse

def update_universe_mass(universe,ffield):
    '''mass is a forcefield dependent property that is not assumed or created at the time
       of universe creation.'''

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

