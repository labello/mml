
import logging
import math

import measure 

# change default level of logger to 10 to print each parmater lookup
logging.basicConfig(filename=__name__,level=30)

def ambergaffenergy(universe):
    ''' Returns a dictionary with floating point values for the energy contributions
    of the five major components of the forcefield (Ebond,Eangle,Etorsion,EVDW,Eel),
    as well as lists with the contributions from each bond, angle, torsion, or atom
    pair (EBonds,EAngles,ETorsions,EVDWs,EEls)
    '''

    from amberff import gaffparams as ffparams
    EBonds,Ebond       = calc_bond_energy(ffparams,universe)
    EAngles,Eangle     = calc_angle_energy(ffparams,universe) 
    ETorsions,Etorsion = calc_torsion_energy(ffparams,universe)
    EVDWs,EEls,EVDW,Eel = calc_nonbonded_energy(ffparams,universe)
    
    E = Ebond + Eangle + Etorsion + EVDW + Eel

    return ({'EBonds':EBonds,'Ebond':Ebond,'EAngles':EAngles,'EAngle':Eangle,
             'ETorsions':ETorsions,'Etorsion':Etorsion,'EVDWs':EVDWs,'EEls':EEls,
             'EVDW':EVDW,'Eel':Eel,'E':E})
    
def amber94energy(universe):
    ''' Returns a dictionary with floating point values for the energy contributions
    of the five major components of the forcefield (Ebond,Eangle,Etorsion,EVDW,Eel),
    as well as lists with the contributions from each bond, angle, torsion, or atom
    pair (EBonds,EAngles,ETorsions,EVDWs,EEls)
    '''
    from amberff import amber94params as ffparams
    EBonds,Ebond       = calc_bond_energy(ffparams,universe)
    EAngles,Eangle     = calc_angle_energy(ffparams,universe) 
    ETorsions,Etorsion = calc_torsion_energy(ffparams,universe)
    EVDWs,EEls,EVDW,Eel = calc_nonbonded_energy(ffparams,universe)
    
    E = Ebond + Eangle + Etorsion + EVDW + Eel

    return ({'EBonds':EBonds,'Ebond':Ebond,'EAngles':EAngles,'EAngle':EAngle,
             'ETorsions':ETorsions,'Etorsion':Etorsion,'EVDWs':EVDWs,'EEls':EEls,
             'EVDW':EVDW,'Eel':Eel,'E':E})
    
def calc_bond_energy(ffparams,universe):
    E = 0.0
    EBonds = []
    for bond in universe.bonds:
        atom1,atom2 = bond.atom1,bond.atom2 
        params = get_bond_params(ffparams,universe,bond)
        r = measure.distance(atom1,atom2)
        e = params.k * (r - params.r0)**2
        EBonds.append((bond,r,e))
        E += e
    return (EBonds,E)
    
    
def calc_angle_energy(ffparams,universe):
    E = 0.0
    EAngles = []
    for angle in universe.angles:
        atom1,atom2,atom3 = angle.atom1,angle.atom2,angle.atom3 
        params = get_angle_params(ffparams,universe,angle)   
        theta = measure.vec_angle(atom1.vec,atom2.vec,atom3.vec)
        # params.k is kcal/mol/radian**2 as read from amber force field 
        # math.radieans is used to convert the angle to radians before evaluation
        e = params.k * (math.radians(theta) - math.radians(params.theta0) )**2
        EAngles.append((angle,theta,e))
        E += e
    return (EAngles,E)

def calc_torsion_energy(ffparams,universe):
    E = 0.0
    ETorsions = []
    counter = 0
    for tor in universe.torsions:
        atom1,atom2,atom3,atom4 = tor.atom1,tor.atom2,tor.atom3,tor.atom4 
        params = get_tor_params(ffparams,universe,tor)
        e = 0
        for term in params:
            Vn2 = term.Vn2
            gamma = term.gamma
            period = term.period
            bondpaths = term.bondpaths
            phi = measure.vec_dihedral(atom1.vec,atom2.vec,atom3.vec,atom4.vec)
            termE = Vn2 * (1 + (math.cos(math.radians( period*phi - gamma ) ) ) )
            e = e + termE
        ETorsions.append((tor,phi,e))
        E += e
    return (ETorsions,E)
    
def calc_nonbonded_energy(ffparams,universe):
    dielectric_constant = 1.0 
    NAtoms = len(universe.atoms)
    EVDWs = []
    EEls  = []
    EVDW = 0
    Eel = 0
    for ai in universe.atoms:
        for aj in universe.atoms:
            if aj not in (ai.fnn + ai.snn):
                iparams = get_vdw_params(ffparams,ai)
                jparams = get_vdw_params(ffparams,aj)
                Rstarij = iparams.R + jparams.R
                epsilonij = math.sqrt(iparams.epsilon * jparams.epsilon)
                Aij = epsilonij * (Rstarij)**12
                Bij = 2 * epsilonij * (Rstarij)**6
                Rij = measure.distance(ai,aj)
                vdw = (Aij / (Rij)**12) - (Bij / (Rij)**6)
                el  = (ai.charge * aj.charge) / (dielectric_constant * Rij)
# apply scale factor for atoms separated by 3 bonds
                if aj in ai.tnn:
                    vdw = vdw * (1.0/2.0) 
                    el  =  el * (1.0/1.2)
            
                EVDWs.append((ai,aj,Rij,vdw))
                EEls.append((ai,aj,Rij,el))
                EVDW += vdw
                Eel += el
    return (EVDWs,EEls,EVDW,Eel)
               

def get_bond_params(ffparams,universe,bond):
    a1 = bond.atom1 
    a2 = bond.atom2 
    dks = ['-'.join([a1.atomtype,a2.atomtype]),
           '-'.join([a2.atomtype,a1.atomtype])]
    params = None
    for dk in dks:
        if dk in ffparams.Bonds:
            return ffparams.Bonds[dk]
    
def get_angle_params(ffparams,universe,angle):
    a1 = angle.atom1 
    a2 = angle.atom2 
    a3 = angle.atom3 
    dks = ['-'.join([a1.atomtype,a2.atomtype,a3.atomtype]),
           '-'.join([a3.atomtype,a2.atomtype,a1.atomtype])]
    for dk in dks:
        if dk in ffparams.Angles:
            return ffparams.Angles[dk]
     
def get_tor_params(ffparams,universe,torsion):
    a1 = torsion.atom1 
    a2 = torsion.atom2 
    a3 = torsion.atom3 
    a4 = torsion.atom4 
    
    dks = ['-'.join([ a1.atomtype,a2.atomtype,a3.atomtype,a4.atomtype]),
    '-'.join([ a4.atomtype,a3.atomtype,a2.atomtype,a1.atomtype]),
    '-'.join([ 'X',a3.atomtype,a2.atomtype,'X']),
    '-'.join([ 'X',a2.atomtype,a3.atomtype,'X']) ]
        
    for dk in dks:
        if dk in ffparams.Torsions:
            return ffparams.Torsions[dk]
    return params

def get_vdw_params(ffparams,atom):
    dk = atom.atomtype
    params = ffparams.VdWs[dk]
    return params
