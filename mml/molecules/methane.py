

from collections import namedtuple

import numpy as np

Atom = namedtuple('Atom','index type x y z vec fnn snn tnn charge')
Bond = namedtuple('Bond','index label i j length')
Angle = namedtuple('Angle','index label i j k angle')
Torsion = namedtuple('Torsion','index label i j k l angle')

Sysdata = namedtuple('Sysdata','Atoms,Bonds,Angles,Torsions')


Atoms = []
Bonds = []
Angles = []
Torsions = []


Atoms.append(Atom(0,'c3',1.011000,-0.074000,0.066000,np.array([1.011000,-0.074000,0.066000]),[1,2,3,4],[0],[1,2,3,4],-0.108800))
Atoms.append(Atom(1,'hc',2.103000,-0.074000,0.066000,np.array([2.103000,-0.074000,0.066000]),[0],[1,2,3,4],[0],0.026700))
Atoms.append(Atom(2,'hc',0.647000,-0.918000,0.657000,np.array([0.647000,-0.918000,0.657000]),[0],[1,2,3,4],[0],0.026700))
Atoms.append(Atom(3,'hc',0.647000,-0.165000,-0.960000,np.array([0.647000,-0.165000,-0.960000]),[0],[1,2,3,4],[0],0.026700))
Atoms.append(Atom(4,'hc',0.647000,0.859000,0.500000,np.array([0.647000,0.859000,0.500000]),[0],[1,2,3,4],[0],0.026700))
Bonds.append(Bond(0,'c3-hc',0,1,111))
Bonds.append(Bond(1,'c3-hc',0,2,111))
Bonds.append(Bond(2,'c3-hc',0,3,111))
Bonds.append(Bond(3,'c3-hc',0,4,111))
Angles.append(Angle(0,'hc-c3-hc',1,0,2,222))
Angles.append(Angle(1,'hc-c3-hc',1,0,3,222))
Angles.append(Angle(2,'hc-c3-hc',1,0,4,222))
Angles.append(Angle(3,'hc-c3-hc',2,0,3,222))
Angles.append(Angle(4,'hc-c3-hc',2,0,4,222))
Angles.append(Angle(5,'hc-c3-hc',3,0,4,222))
inputdata = Sysdata(Atoms,Bonds,Angles,Torsions)
