'''
Define data types and some basic operations on that data.
'''

import numpy as np
import random

import constants
import measure

class Atom:
    '''
    A basic container for atoms.

    =================== ==================================================================
    Attribute           Description
    =================== ==================================================================
    index               Unique integer identifier of the atom
    symbol              The 1-3 letter code for a chemical element
    atomtype            A string representation for forcefield lookup, e.g., 'c3'
    mass                The atomic mass 
    charge              The formal charge of the atom
    label               A string used to identify individual atoms, e.g., C1
    fnn                 List of atoms one bond away (1st nearest neighbors)
    snn                 List of atoms two bonds away (2nd nearest neighbors)
    tnn                 List of atoms three bonds away 3rd nearest neighbors)
    x                   The atom's X Cartesian coordinate
    y                   The atom's Y Cartesian coordinate
    z                   The atom's Z Cartesian coordinate
    vec                 A vector of the atom cooridinates <x,y,z> 
    vx                  The atom's velocity on the X axis
    vy                  The atom's velocity on the Y axis
    vz                  The atom's velocity on the Z axis
    vvec                A vector of the atom velocity
    vmean               The mean velocity in all directions 
    update_coords()
    '''

    atom_id = 0

    def __init__(self,symbol=None,atomtype=None,mass=None,charge=0,label=None,
                fnn=[],snn=[],tnn=[],x=None,y=None,z=None):

        self.index = Atom.atom_id
        self.symbol = symbol
        self.atomtype = atomtype
        self.mass = mass
        self.charge = charge
        self.label = label
        self.fnn = fnn
        self.snn = snn
        self.tnn = tnn
        self.x = x
        self.y = y
        self.z = z
        self.vec = np.array([x,y,z])
        self.vx = 0.0
        self.vy = 0.0
        self.vz = 0.0
        self.vvec = np.array([self.vx,self.vy,self.vz])
        self.vmean = np.linalg.norm(self.vvec)
# Increment the atom_id counter each time a new Atom instance is created
        Atom.atom_id += 1 
         
    def update_coords(self,x,y,z):
        '''Make changes to the coordinates of the atoms using this function to ensure that
           all Cartesian related properties remain in sync.
           
           *Usage
           for atom in universe:
               atom.update_coords(newx,newy,newz) 
        '''
        self.x = x
        self.y = y
        self.z = z
        self.vec = np.array([x,y,z])

    def update_vels(self,vx,vy,vz):
        '''Make changes to the velocities of the atoms using this function to ensure that
           all velocity related properties remain in sync.
           
           *Usage
           for atom in universe:
               atom.update_coords(newvx,newvy,newvz) 

        Velocity should be reported to this function in units of 
        Angstroms/picosecond.
        '''
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.vvec = np.array([vx,vy,vz])
        self.vmean = np.linalg.norm(self.vvec)

class Bond:
    '''
    A bond object is a container for atoms that participate in the bond, plus some 
    additional properties that are unique to the bond. 

    =================== =================================================================
    Attribute           Description
    =================== =================================================================
    index               Unique integer identifier of the bond
    bondtype            A string representation for forcefield lookup, e.g, "c3-c3"
    atoms               A list of Atom types
    atom1               The first atom in the list
    atom2               The second atom in the list
    i                   Convenience attribute to reference atom1.index
    j                   Convenience attribute to reference atom2.index
    length              The distance between atom1 and atom2

    '''

    bond_id = 0

    def __init__(self,bondtype=None,atoms=[]):
        self.index = Bond.bond_id
        self.atoms = atoms
        self.atom1 = atoms[0]
        self.atom2 = atoms[1]
        self.i = self.atom1.index
        self.j = self.atom2.index
        self.bondtype = '-'.join([self.atom1.atomtype,self.atom2.atomtype])

# Increment the bond_id counter each time a new Bond instance is created
        Bond.bond_id += 1

    def length(self):
        return measure.distance(self.atom1,self.atom2)

class Angle:
    '''
    An angle object is a container for atoms that participate in the angle, plus some 
    additional properties that are unique to the angle. 

    =================== =================================================================
    Attribute           Description
    =================== =================================================================
    index               Unique integer identifier of the bond
    angletype           A string representation for forcefield lookup, e.g, "c3-c3-c3"
    atoms               A list of Atom types
    atom1               The first atom in the list
    atom2               The second atom in the list
    atom3               The third atom in the list
    i                   Convenience attribute to reference atom1.index
    j                   Convenience attribute to reference atom2.index
    k                   Convenience attribute to reference atom3.index
    angle               The atom1--atom2--atom3 angle

    '''

    angle_id = 0

    def __init__(self,angletype=None,atoms=[]):
        self.index = Angle.angle_id 
        self.atoms = atoms
        self.atom1 = atoms[0]
        self.atom2 = atoms[1]
        self.atom3 = atoms[2]
        self.i = self.atom1.index
        self.j = self.atom2.index
        self.k = self.atom3.index
        self.angletype = '-'.join([self.atom1.atomtype,self.atom2.atomtype,self.atom3.atomtype])
# Increment the angle_id counter each time a new Angle instance is created
        Angle.angle_id += 1
        
    def angle(self):
        return measure.vec_angle(self.atom1.coords,self.atom2.coords,self.atom3.coords)

class Torsion:
    '''
    A torsion object is a container for atoms that participate in the torsion, plus some 
    additional properties that are unique to the torsion. 

    =================== =================================================================
    Attribute           Description
    =================== =================================================================
    index               Unique integer identifier of the torsion
    torsiontype         A string representation for forcefield lookup, e.g, "c3-c3-c3-c3"
    atoms               A list of Atom types
    atom1               The first atom in the list
    atom2               The second atom in the list
    atom3               The third atom in the list
    atom4               The forth atom in the list
    i                   Convenience attribute to reference atom1.index
    j                   Convenience attribute to reference atom2.index
    k                   Convenience attribute to reference atom3.index
    l                   Convenience attribute to reference atom4.index
    torsion             The atom1-atom2---atom3-atom4 dihedral angle

    '''

    torsion_id = 0

    def __init__(self,torsiontype=None,atoms=[]):
        self.index = Torsion.torsion_id
        self.atoms = atoms
        self.atom1 = atoms[0]
        self.atom2 = atoms[1]
        self.atom3 = atoms[2]
        self.atom4 = atoms[3]
        self.i = self.atom1.index
        self.j = self.atom2.index
        self.k = self.atom3.index
        self.l = self.atom4.index
        self.torsiontype = '-'.join([self.atom1.atomtype,self.atom2.atomtype,self.atom3.atomtype,self.atom4.atomtype]) 

# Increment the angle_id counter each time a new Angle instance is created
        Torsion.torsion_id += 1

    def torsion(self):
        return measure.vec_dihedral(self.atom1.coords,self.atom2.coords,
                                    self.atom3.coords,self.atom4.coords)

class Universe:
    '''
    A universe object is a container for all particles, boundaries, and other entities
    considered in the calculation.

    =============================  =================================================================
    Attribute                      Description
    =============================  =================================================================
    atoms
    bonds
    angles
    torsions
    natoms()                       Calculate the number of atoms in the universe
    nparameters()                  Number of cartesian parameters: 3 * natoms()
    get_cartesian_atom_array()     Return a numpy array with the cartesian coordinates of all atoms
    get_flat_carteisan_array()     Return a vecotr of the cartesian coordinates
    get_cartesian_velocity_array() Description
    assign_random_velocity()       Assign random velocity to all atoms in Universe
    temperature()                  Calculate the temperature of the atoms in the Universe 
    write_xyz()                    Description
    update_neighbor_assignments    Description

    '''

    def __init__(self,atoms=[],bonds=[],angles=[],torsions=[]):
        self.atoms = atoms
        self.bonds = bonds
        self.angles = angles
        self.torsions = torsions

    def natoms(self):
        return len(self.atoms)    

    def nparameters(self):
        return (self.natoms() * 3 )
  

    def get_cartesian_atom_array(self):
        '''Return a single NATOMSx3 Numpy array with the Cartesian coordinates'''
        return np.vstack(atom.vec for atom in self.atoms)

    def get_centered_cartesian_atom_array(self):
        '''Return a single NATOMSx3 Numpy array with the center-of-mass centered Cartesian coordinates.'''
        pass  

    def get_flat_cartesian_atom_array(self):
        '''Return a vector of the cartesian coordinates.  Useful for the built-in SciPy
        routines.    
        '''
        return np.ndarray.flatten(self.get_cartesian_atom_array())

    def get_cartesian_velocity_array(self):
        '''Return a single NATOMSx3 Numpy array with atom velocities.'''
        return np.vstack(atom.vvec for atom in self.atoms)

    def assign_random_velocity(self,absmax=100):
        '''Assign random starting velocity in angstroms/picosecond to all atoms in Universe.  
           absmax: the max absolute value of a velocity in one cartesian direction 
        '''
        for atom in self.atoms:
            vx = random.uniform(-absmax,absmax)
            vy = random.uniform(-absmax,absmax)
            vz = random.uniform(-absmax,absmax)
            atom.update_vels(vx,vy,vz)


    def temperature(self):
        ''' T  =   mv**2
                 [ ----- ]   /  3k
                     N
        N = number of atoms
        k = Boltzman's constant = 1.380658 x 10^-23 kg * m**2 / Kelvin * s**2  

        INCOMING VELOCITY ASSUMED TO BE ANGSTROMS/PICOSECOND
        '''
      
        mvsq = 0
        for atom in self.atoms:
            m = atom.mass 
            v = atom.vmean 
            mvsq = mvsq + (m * v**2)  
        mvsqav = mvsq / self.natoms()
        temperature = mvsqav / (3 * constants.k3)
        return temperature

    def write_xyz(self,filehandle='MINIMIZE.xyz'):
        '''
        Write current coordinates of all atoms in the universe out to a file.  Append to the file if it
        already exists.  Writes in xyz format.  
        '''
        minfile = open(filehandle,'a')
        minfile.write( str(self.natoms() ) )
        minfile.write( '\n\n' )
        for atom in self.atoms: 
            minfile.write( "{0:4s} {1:12.5f} {2:12.5f} {3:12.5f} \n".format(
             atom.atomtype,atom.x,atom.y,atom.z) )

    def update_neighbor_assignments(self):
        '''
        This method updates the fnn, snn, and tnn fields of atoms in the universe. It should
        be run following a change in the universe atoms or bonds. Be sure to run it if the
        universe is instantiated empty, with atoms and bonds added later.  

        universe = Universe()
      # add atoms and bonds to universe.atoms and universe.bonds 
        universe.update_neighbor_assignments()

        '''
        for atom in self.atoms:
            atom.fnn = self.__get_first_nearest_neighbors(atom) 
            atom.snn = self.__get_second_nearest_neighbors(atom)
            atom.tnn = self.__get_third_nearest_neighbors(atom)

    def __get_first_nearest_neighbors(self,atom):
        FNN = []
        for bond in self.bonds:
            if atom in bond.atoms:
                if atom == bond.atom1:
                    FNN.append(bond.atom2)
                else:
                    FNN.append(bond.atom1)
        return FNN

    def __get_second_nearest_neighbors(self,atom):
        FNN = self.__get_first_nearest_neighbors(atom)
        SNN = []
        for fnn in FNN:
            SNN += self.__get_first_nearest_neighbors(fnn)
        SNN = list(set(SNN))
        return SNN

    def __get_third_nearest_neighbors(self,atom):
        SNN = self.__get_second_nearest_neighbors(atom)
        TNN = []
        for snn in SNN:
            TNN += self.__get_first_nearest_neighbors(snn)
        TNN = list(set(TNN))
        return TNN

