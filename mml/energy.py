#!/usr/bin/env python

import argparse
import sys
import numpy as np

import amber
import measure 
import read_input


def main(molfile=None,ffield='gaff',noreport=False):
    ''' Calculate the energy for the specified molecular universe and forcefield.
        
        * molecule: a module (e.g., ethane.py) from the molecules directory. 
            The molecule data (bond tables, coordinates, angles, torsions, etc.,
            have been assembled into a Python data structure stored in the .py
            file.  
`
        * ffield:   gaff or amber94
    '''
    universe = read_input.main(molfile)
 
    coordinates = input2cartesian_array(universe.atoms)  

    energy(ffield,universe,noreport)

def energy(ffield,universe,noreport=True):


    if ffield == 'gaff':
        results = amber.ambergaffenergy(universe)
        E = results['E']

    elif ffield == 'amber94':
        results = amber.amber94energy(universe)
        E = results['E']

    if not noreport:
       print_report(universe,results)

    return E

def print_report(universe,results):
    
    ldelim = '\n' + 60*'*' + '\n'
    
    header = ldelim + "Molecular Modeling Lite: Energy (Version 0.01)" + ldelim
    bondheader =  ldelim + "   Contributions from individual bonds"    + ldelim
    angleheader = ldelim + "   Contributions from individual angles"   + ldelim
    torheader =   ldelim + "   Contributions from individual torsions" + ldelim
    vdwheader =   ldelim + "   Contributions from Van der Waals pairs" + ldelim
    elheader =    ldelim + "   Contributions from electrostatic pairs" + ldelim
    
    print header
    
    template = '{0:5d} {1:3s} {2:8.5f} {3:8.5f} {4:8.5f}'
    for atom in universe.atoms:
        print template.format(atom.index,atom.atomtype,atom.x,atom.y,atom.z) 
         
    print bondheader
    
    template = '{0:4d} {1:4d} {2:8s} {3:4.3f} Angs {4:4.3f} kcal/mol'
    for bondresult in results['EBonds']:
        bond = bondresult[0]
        r = bondresult[1]
        e = bondresult[2]
        print template.format(bond.i,bond.j,bond.bondtype,r,e)
        
    print angleheader 
    
    template = '{0:4d} {1:4d} {2:4d} {3:8s} {4:>7.3f} deg {5:>7.3f} kcal/mol'
    for angleresult in results['EAngles']:
        angle = angleresult[0]
        theta = angleresult[1]
        e = angleresult[2]
        print template.format(angle.i,angle.j,angle.k,angle.angletype,theta,e)

    print torheader

    template = '{0:4d} {1:4d} {2:4d} {3:4d} {4:11s} {5:>7.2f} deg {6:>7.3f} kcal/mol'
    for torresult in results['ETorsions']:
        tor = torresult[0]
        phi = torresult[1]
        e = torresult[2]
        print template.format(tor.i,tor.j,tor.k,tor.l,tor.torsiontype,phi,e)

    print vdwheader
 
    template = "{0:4d} {1:4d} {2:3s} {3:3s} {4:>5.2f} Angs {5:>7.3f} kcal/mol"
    for vdwpair in results['EVDWs']:
        ai,aj,r,e = vdwpair
        print template.format(ai.index,aj.index,ai.atomtype,aj.atomtype,r,e)
         
    print elheader
    
    for elpair in results['EEls']:
        ai,aj,r,e = elpair
        print template.format(ai.index,aj.index,ai.atomtype,aj.atomtype,r,e)

    print"\n"
    print "E(BONDS)    %12.6f kcal/mol" % results['Ebond']
    print "E(ANGLES)   %12.6f kcal/mol" % results['EAngle']
    print "E(TORSIONS) %12.6f kcal/mol" % results['Etorsion']
    print "E(VDW)      %12.6f kcal/mol" % results['EVDW']
    print "E(EL)       %12.6f kcal/mol" % results['Eel']
    print ""
    print "E(TOTAL): %f kcal/mol \n" % (results['E'])


def input2cartesian_array(Atoms):
    return np.vstack(atom.vec for atom in Atoms )

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='\n\n    Molecular Modeling Lite\n',add_help=True)
    parser.add_argument('--molfile',default='/Users/labello/gitwork/mml/no-git/molecules/acfiles/ethane.mol.ac',
     help='Path to .ac molecule file (Antechamber format)')
    parser.add_argument('--ffield',default='gaff',help='amber94 or gaff')
    parser.add_argument('--noreport',action='store_true',help='Do not print report.')
    args = parser.parse_args()
    main(args.molfile,args.ffield,args.noreport)

