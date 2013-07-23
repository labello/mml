

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


Atoms.append(Atom(0,'o',0.957000,-0.109000,-0.082000,np.array([0.957000,-0.109000,-0.082000]),[1],[0,2,3],[1,4,5,6,7],-0.496000))
Atoms.append(Atom(1,'c',2.168000,-0.033000,0.037000,np.array([2.168000,-0.033000,0.037000]),[0,3,2],[1,4,5,6,7],[0,8,2,3,9],0.624100))
Atoms.append(Atom(2,'oh',2.936000,0.871000,-0.627000,np.array([2.936000,0.871000,-0.627000]),[1,5],[0,2,3],[1,4,5,6,7],-0.575100))
Atoms.append(Atom(3,'c3',3.015000,-0.902000,0.964000,np.array([3.015000,-0.902000,0.964000]),[1,4,6,7],[0,8,2,3,9],[1,4,5,6,7],0.062800))
Atoms.append(Atom(4,'n3',4.452000,-0.734000,0.619000,np.array([4.452000,-0.734000,0.619000]),[3,8,9],[1,4,6,7],[0,8,2,3,9],-0.932800))
Atoms.append(Atom(5,'ho',3.863000,0.709000,-0.322000,np.array([3.863000,0.709000,-0.322000]),[2],[1,5],[0,2,3],0.444000))
Atoms.append(Atom(6,'h1',2.729000,-1.951000,0.849000,np.array([2.729000,-1.951000,0.849000]),[3],[1,4,6,7],[0,8,2,3,9],0.067700))
Atoms.append(Atom(7,'h1',2.848000,-0.570000,1.992000,np.array([2.848000,-0.570000,1.992000]),[3],[1,4,6,7],[0,8,2,3,9],0.067700))
Atoms.append(Atom(8,'hn',4.653000,-1.343000,-0.183000,np.array([4.653000,-1.343000,-0.183000]),[4],[8,9,3],[1,4,6,7],0.368800))
Atoms.append(Atom(9,'hn',5.001000,-1.142000,1.382000,np.array([5.001000,-1.142000,1.382000]),[4],[8,9,3],[1,4,6,7],0.368800))
Bonds.append(Bond(0,'o-c',0,1,111))
Bonds.append(Bond(1,'c-c3',1,3,111))
Bonds.append(Bond(2,'c-oh',1,2,111))
Bonds.append(Bond(3,'oh-ho',2,5,111))
Bonds.append(Bond(4,'c3-n3',3,4,111))
Bonds.append(Bond(5,'c3-h1',3,6,111))
Bonds.append(Bond(6,'c3-h1',3,7,111))
Bonds.append(Bond(7,'n3-hn',4,8,111))
Bonds.append(Bond(8,'n3-hn',4,9,111))
Angles.append(Angle(0,'o-c-c3',0,1,3,222))
Angles.append(Angle(1,'o-c-oh',0,1,2,222))
Angles.append(Angle(2,'c3-c-oh',3,1,2,222))
Angles.append(Angle(3,'c-c3-n3',1,3,4,222))
Angles.append(Angle(4,'c-c3-h1',1,3,6,222))
Angles.append(Angle(5,'c-c3-h1',1,3,7,222))
Angles.append(Angle(6,'c-oh-ho',1,2,5,222))
Angles.append(Angle(7,'n3-c3-h1',4,3,6,222))
Angles.append(Angle(8,'n3-c3-h1',4,3,7,222))
Angles.append(Angle(9,'c3-n3-hn',3,4,8,222))
Angles.append(Angle(10,'c3-n3-hn',3,4,9,222))
Angles.append(Angle(11,'h1-c3-h1',6,3,7,222))
Angles.append(Angle(12,'hn-n3-hn',8,4,9,222))
Torsions.append(Torsion(0,'o-c-c3-n3',0,1,3,4,333))
Torsions.append(Torsion(1,'o-c-c3-h1',0,1,3,6,333))
Torsions.append(Torsion(2,'o-c-c3-h1',0,1,3,7,333))
Torsions.append(Torsion(3,'o-c-oh-ho',0,1,2,5,333))
Torsions.append(Torsion(4,'n3-c3-c-oh',4,3,1,2,333))
Torsions.append(Torsion(5,'h1-c3-c-oh',6,3,1,2,333))
Torsions.append(Torsion(6,'h1-c3-c-oh',7,3,1,2,333))
Torsions.append(Torsion(7,'c3-c-oh-ho',3,1,2,5,333))
Torsions.append(Torsion(8,'o-c-c3-n3',0,1,3,4,333))
Torsions.append(Torsion(9,'c-c3-n3-hn',1,3,4,8,333))
Torsions.append(Torsion(10,'c-c3-n3-hn',1,3,4,9,333))
Torsions.append(Torsion(11,'o-c-c3-h1',0,1,3,6,333))
Torsions.append(Torsion(12,'o-c-c3-h1',0,1,3,7,333))
Torsions.append(Torsion(13,'o-c-oh-ho',0,1,2,5,333))
Torsions.append(Torsion(14,'c3-c-oh-ho',3,1,2,5,333))
Torsions.append(Torsion(15,'hn-n3-c3-h1',8,4,3,6,333))
Torsions.append(Torsion(16,'hn-n3-c3-h1',9,4,3,6,333))
Torsions.append(Torsion(17,'hn-n3-c3-h1',8,4,3,7,333))
Torsions.append(Torsion(18,'hn-n3-c3-h1',9,4,3,7,333))
Torsions.append(Torsion(19,'c-c3-n3-hn',1,3,4,8,333))
Torsions.append(Torsion(20,'c-c3-n3-hn',1,3,4,9,333))
inputdata = Sysdata(Atoms,Bonds,Angles,Torsions)
