

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


Atoms.append(Atom(0,'c3',1.012000,-0.100000,-0.006000,np.array([1.012000,-0.100000,-0.006000]),[1,2,3,4],[0,5,6,7],[1,2,3,4],-0.094100))
Atoms.append(Atom(1,'c3',2.523000,-0.100000,-0.007000,np.array([2.523000,-0.100000,-0.007000]),[0,5,6,7],[1,2,3,4],[0,5,6,7],-0.094100))
Atoms.append(Atom(2,'hc',0.627000,0.792000,-0.511000,np.array([0.627000,0.792000,-0.511000]),[0],[1,2,3,4],[0,5,6,7],0.031700))
Atoms.append(Atom(3,'hc',0.627000,-0.109000,1.018000,np.array([0.627000,-0.109000,1.018000]),[0],[1,2,3,4],[0,5,6,7],0.031700))
Atoms.append(Atom(4,'hc',0.627000,-0.982000,-0.527000,np.array([0.627000,-0.982000,-0.527000]),[0],[1,2,3,4],[0,5,6,7],0.031700))
Atoms.append(Atom(5,'hc',2.908000,0.783000,0.514000,np.array([2.908000,0.783000,0.514000]),[1],[0,5,6,7],[1,2,3,4],0.031700))
Atoms.append(Atom(6,'hc',2.908000,-0.090000,-1.031000,np.array([2.908000,-0.090000,-1.031000]),[1],[0,5,6,7],[1,2,3,4],0.031700))
Atoms.append(Atom(7,'hc',2.908000,-0.991000,0.497000,np.array([2.908000,-0.991000,0.497000]),[1],[0,5,6,7],[1,2,3,4],0.031700))
Bonds.append(Bond(0,'c3-c3',0,1,111))
Bonds.append(Bond(1,'c3-hc',0,2,111))
Bonds.append(Bond(2,'c3-hc',0,3,111))
Bonds.append(Bond(3,'c3-hc',0,4,111))
Bonds.append(Bond(4,'c3-hc',1,5,111))
Bonds.append(Bond(5,'c3-hc',1,6,111))
Bonds.append(Bond(6,'c3-hc',1,7,111))
Angles.append(Angle(0,'c3-c3-hc',1,0,2,222))
Angles.append(Angle(1,'c3-c3-hc',1,0,3,222))
Angles.append(Angle(2,'c3-c3-hc',1,0,4,222))
Angles.append(Angle(3,'c3-c3-hc',0,1,5,222))
Angles.append(Angle(4,'c3-c3-hc',0,1,6,222))
Angles.append(Angle(5,'c3-c3-hc',0,1,7,222))
Angles.append(Angle(6,'hc-c3-hc',2,0,3,222))
Angles.append(Angle(7,'hc-c3-hc',2,0,4,222))
Angles.append(Angle(8,'hc-c3-hc',3,0,4,222))
Angles.append(Angle(9,'hc-c3-hc',5,1,6,222))
Angles.append(Angle(10,'hc-c3-hc',5,1,7,222))
Angles.append(Angle(11,'hc-c3-hc',6,1,7,222))
Torsions.append(Torsion(0,'hc-c3-c3-hc',5,1,0,2,333))
Torsions.append(Torsion(1,'hc-c3-c3-hc',6,1,0,2,333))
Torsions.append(Torsion(2,'hc-c3-c3-hc',7,1,0,2,333))
Torsions.append(Torsion(3,'hc-c3-c3-hc',5,1,0,3,333))
Torsions.append(Torsion(4,'hc-c3-c3-hc',6,1,0,3,333))
Torsions.append(Torsion(5,'hc-c3-c3-hc',7,1,0,3,333))
Torsions.append(Torsion(6,'hc-c3-c3-hc',5,1,0,4,333))
Torsions.append(Torsion(7,'hc-c3-c3-hc',6,1,0,4,333))
Torsions.append(Torsion(8,'hc-c3-c3-hc',7,1,0,4,333))
inputdata = Sysdata(Atoms,Bonds,Angles,Torsions)
