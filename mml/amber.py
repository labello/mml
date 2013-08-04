
import logging
import math

from molecules import ethane,propane
import measure 

# change default level of logger to 10 to print each parmater lookup
logging.basicConfig(filename=__name__,level=30)

def ambergaffenergy(coordinates,system=ethane.inputdata):
    ''' Returns a dictionary with floating point values for the energy contributions
    of the five major components of the forcefield (Ebond,Eangle,Etorsion,EVDW,Eel),
    as well as lists with the contributions from each bond, angle, torsion, or atom
    pair (EBonds,EAngles,ETorsions,EVDWs,EEls)
    '''

    from amberff import gaffparams as ffparams
    EBonds,Ebond       = calc_bond_energy(ffparams,system,coordinates)
    EAngles,Eangle     = calc_angle_energy(ffparams,system,coordinates) 
    ETorsions,Etorsion = calc_torsion_energy(ffparams,system,coordinates)
    EVDWs,EEls,EVDW,Eel = calc_nonbonded_energy(ffparams,system,coordinates)
    
    E = Ebond + Eangle + Etorsion + EVDW + Eel

    return ({'EBonds':EBonds,'Ebond':Ebond,'EAngles':EAngles,'EAngle':Eangle,
             'ETorsions':ETorsions,'Etorsion':Etorsion,'EVDWs':EVDWs,'EEls':EEls,
             'EVDW':EVDW,'Eel':Eel,'E':E})
    
def amber94energy(coordinates,system=ethane.inputdata):
    ''' Returns a dictionary with floating point values for the energy contributions
    of the five major components of the forcefield (Ebond,Eangle,Etorsion,EVDW,Eel),
    as well as lists with the contributions from each bond, angle, torsion, or atom
    pair (EBonds,EAngles,ETorsions,EVDWs,EEls)
    '''
    from amberff import amber94params as ffparams
    EBonds,Ebond       = calc_bond_energy(ffparams,system,coordinates)
    EAngles,Eangle     = calc_angle_energy(ffparams,system,coordinates) 
    ETorsions,Etorsion = calc_torsion_energy(ffparams,system,coordinates)
    EVDWs,EEls,EVDW,Eel = calc_nonbonded_energy(ffparams,system,coordinates)
    
    E = Ebond + Eangle + Etorsion + EVDW + Eel

    return ({'EBonds':EBonds,'Ebond':Ebond,'EAngles':EAngles,'EAngle':EAngle,
             'ETorsions':ETorsions,'Etorsion':Etorsion,'EVDWs':EVDWs,'EEls':EEls,
             'EVDW':EVDW,'Eel':Eel,'E':E})
    
def calc_bond_energy(ffparams,system,coordinates):
    E = 0.0
    EBonds = []
    for bond in system.Bonds:
        atom1,atom2 = coordinates[bond.i],coordinates[bond.j]
        params = get_bond_params(ffparams,system,bond)
        r = measure.distance(atom1,atom2)
        e = params.k * (r - params.r0)**2
        EBonds.append((bond,r,e))
        E += e
    return (EBonds,E)
    
    
def calc_angle_energy(ffparams,system,coordinates):
    E = 0.0
    EAngles = []
    for angle in system.Angles:
        atom1,atom2,atom3 = coordinates[angle.i],coordinates[angle.j],coordinates[angle.k]
        params = get_angle_params(ffparams,system,angle)   
        k = params.k
        theta0 = params.theta0
        theta = measure.vec_angle(atom1,atom2,atom3)
        e = k * (theta - theta0)**2
        EAngles.append((angle,theta,e))
        E += e
    return (EAngles,E)

def calc_torsion_energy(ffparams,system,coordinates):
    E = 0.0
    ETorsions = []
    for tor in system.Torsions:
        atom1,atom2,atom3,atom4 = coordinates[tor.i],coordinates[tor.j],coordinates[tor.k],coordinates[tor.l]
        params = get_tor_params(ffparams,system,tor)
        Vn2 = params.Vn2
        gamma = params.gamma
        period = params.period
        bondpaths = params.bondpaths
        phi = measure.vec_dihedral(atom1,atom2,atom3,atom4)
        e = Vn2 * (1 + math.cos(period*phi - gamma))
        ETorsions.append((tor,phi,e))
        E += e
    return (ETorsions,E)
    
def calc_nonbonded_energy(ffparams,system,coordinates):
    NAtoms = len(system.Atoms)
    EVDWs = []
    EEls  = []
    EVDW = 0
    Eel = 0
    for i in xrange(0,NAtoms):
        for j in xrange(0,i):
            ai,aj = system.Atoms[i],system.Atoms[j]
            if j not in (ai.fnn + ai.snn):
                iparams = get_vdw_params(ffparams,ai)
                jparams = get_vdw_params(ffparams,aj)
                Rstarij = iparams.R + jparams.R
                epsilonij = math.sqrt(iparams.epsilon * jparams.epsilon)
                Aij = epsilonij * (Rstarij)**12
                Bij = 2 * epsilonij * (Rstarij)**6
                Rij = measure.distance(coordinates[i],coordinates[j])
                vdw = (Aij / (Rij)**12) - (Bij / (Rij)**6)
                el  = (ai.charge * aj.charge) / (epsilonij * Rij)
# apply scale factor for atoms separated by 3 bonds
                if j in ai.tnn:
                    vdw = vdw * (1.0/2.0) 
                    el  =  el * (1.0/1.2)
            
                EVDWs.append((ai,aj,Rij,vdw))
                EEls.append((ai,aj,Rij,el))
                EVDW += vdw
                Eel += el
    return (EVDWs,EEls,EVDW,Eel)
               

def get_bond_params(ffparams,system,bond):
    a1 = system.Atoms[bond.i]
    a2 = system.Atoms[bond.j]
    dks = ['-'.join([a1.type,a2.type]),
           '-'.join([a2.type,a1.type])]
    params = None
    for dk in dks:
        if dk in ffparams.Bonds:
            return ffparams.Bonds[dk]
    
def get_angle_params(ffparams,system,angle):
    a1 = system.Atoms[angle.i]
    a2 = system.Atoms[angle.j]
    a3 = system.Atoms[angle.k]
    dks = ['-'.join([a1.type,a2.type,a3.type]),
           '-'.join([a3.type,a2.type,a1.type])]
    for dk in dks:
        if dk in ffparams.Angles:
            return ffparams.Angles[dk]
     
def get_tor_params(ffparams,system,torsion):
    a1 = system.Atoms[torsion.i]
    a2 = system.Atoms[torsion.j]
    a3 = system.Atoms[torsion.k]
    a4 = system.Atoms[torsion.l] 
    
    dks = ['-'.join([ a1.type,a2.type,a3.type,a4.type]),
    '-'.join([ a4.type,a3.type,a2.type,a1.type]),
    '-'.join([ 'X',a3.type,a2.type,'X']),
    '-'.join([ 'X',a2.type,a3.type,'X']) ]
        
    for dk in dks:
        if dk in ffparams.Torsions:
            return ffparams.Torsions[dk]
    return params

def get_vdw_params(ffparams,atom):
    dk = atom.type
    params = ffparams.VdWs[dk]
    return params
