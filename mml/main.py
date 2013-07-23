
import argparse
import logging
import sys
from collections import namedtuple
import math
import numpy as np
import scipy.optimize


import amber
import measure 

# define molecular system input
import molecules
system = molecules.benzene.inputdata

def main():
  
    coordinates = input2cartesian_array(system.Atoms)  
    print_cartesian_output(system.Atoms,coordinates)    
    print amber.ambergaffenergy(coordinates,system)

    flat_coords = np.ndarray.flatten(coordinates)

    res = scipy.optimize.minimize(wrap_energy_for_minimizer,flat_coords,args=(system,),
     method='CG',tol=5.0,bounds=None)
     
    print res

def wrap_energy_for_minimizer(flat_coords,system):
    c = flat_coords.reshape(-1,3)
    print_cartesian_output(system.Atoms,c)
    e = amber.ambergaffenergy(c,system)
#    print e
    return e
    
def print_cartesian_output(Atoms,coordinates):
    print len(coordinates)
    print 
    for i in range(0,len(Atoms)):
        # replace w/ log
        print "{0:4s} {1:12.5f} {2:12.5f} {3:12.5f}".format(
         Atoms[i].type, coordinates[i][0],coordinates[i][1],coordinates[i][2])
    
def input2cartesian_array(Atoms):
    return np.vstack(atom.vec for atom in Atoms )
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='\n\n    Molecular Modeling Lite\n',add_help=True)
    parser.add_argument('--molecule',help='The molecular system from ../mml/molecules')
    parser.add_argument('--routine',default='energy',help='energy, geometry, or dynamics')
    parser.add_argument('--tol',default=5,help='Energy tolerance for gemoetry optimization')
    parser.add_argument('--outputlevel',default=30,help='''The priority of the 
     output that is generated.  Lower numbers correspond to lower priority 
     messages being printed.  Meaningful values are 10 (most output) , 20, 30,
     40, or 50 (least output)''')
    args = parser.parse_args()

    main()
