
import argparse
import importlib
import logging
import sys
import numpy as np

import amber
import measure 


def main(molecule,ffield,report=False,outputlevel=50):
    ''' Calculate the energy for the specified molecular system and forcefield.
        
        * molecule: a module (e.g., ethane.py) from the molecules directory. 
            The molecule data (bond tables, coordinates, angles, torsions, etc.,
            have been assembled into a Python data structure stored in the .py
            file.  
`
        * ffield:   gaff or amber94
    '''

# import the appropriate molecular system
    importmodule = importlib.import_module( 'molecules.' + molecule )
    system = importmodule.inputdata
  
    coordinates = input2cartesian_array(system.Atoms)  

    E = amber.ambergaffenergy(coordinates,system)

    if report:
        print_report(system,E)


    return E

    
def print_report(system,E):
    
    print "\n" 
    print 60*'*'
    print "Molecular Modeling Lite: Energy (Version 0.01)"
    print 60*'*'

    coordinates = input2cartesian_array(system.Atoms)  
    print len(coordinates)
    print 
    for i in range(0,len(system.Atoms)):
        print "{0:4s} {1:12.5f} {2:12.5f} {3:12.5f}".format(
         system.Atoms[i].type, coordinates[i][0],coordinates[i][1],coordinates[i][2])

    print "\n\nEnergy {0:10.3f} kcal/mol\n\n".format(E)

def input2cartesian_array(Atoms):
    return np.vstack(atom.vec for atom in Atoms )

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='\n\n    Molecular Modeling Lite\n',add_help=True)
    parser.add_argument('--molecule',default='ethane',help='The molecular system from ../mml/molecules')
    parser.add_argument('--ffield',default='gaff',help='amber94 or gaff')
    parser.add_argument('--report',default=True ,help='Print a detailed report')
    parser.add_argument('--outputlevel',default=10,help='''The priority of the 
     output that is generated.  Lower numbers correspond to lower priority 
     messages being printed.  Meaningful values are 10 (most output) , 20, 30,
     40, or 50 (least output)''')
    args = parser.parse_args()
    main( args.molecule,args.ffield,args.report,args.outputlevel)

