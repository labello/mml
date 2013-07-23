

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


Atoms.append(Atom(0,'n3',3.197000,-1.028000,-0.728000,np.array([3.197000,-1.028000,-0.728000]),[9,10,1],[0,2,3,6],[1,4,5,7,8,9,10,12],-0.915800))
Atoms.append(Atom(1,'c3',2.544000,0.029000,0.103000,np.array([2.544000,0.029000,0.103000]),[3,2,6,0],[1,4,5,7,8,9,10,12],[0,2,3,6,11,13,14,15,16,17,18],0.076500))
Atoms.append(Atom(2,'h1',2.901000,1.000000,-0.261000,np.array([2.901000,1.000000,-0.261000]),[1],[0,2,3,6],[1,4,5,7,8,9,10,12],0.082700))
Atoms.append(Atom(3,'c',3.046000,-0.113000,1.553000,np.array([3.046000,-0.113000,1.553000]),[1,5,4],[0,11,2,3,6],[1,4,5,7,8,9,10,12],0.629100))
Atoms.append(Atom(4,'oh',3.686000,-1.299000,1.736000,np.array([3.686000,-1.299000,1.736000]),[3,11],[1,4,5],[0,11,2,3,6],-0.578100))
Atoms.append(Atom(5,'o',2.881000,0.708000,2.441000,np.array([2.881000,0.708000,2.441000]),[3],[1,4,5],[0,11,2,3,6],-0.498000))
Atoms.append(Atom(6,'c3',1.000000,-0.015000,0.065000,np.array([1.000000,-0.015000,0.065000]),[1,8,7,12],[0,2,3,6,13,14,15,16,17,18],[1,4,5,7,8,9,10,12],-0.097700))
Atoms.append(Atom(7,'c3',0.457000,0.406000,-1.305000,np.array([0.457000,0.406000,-1.305000]),[6,13,14,15],[8,1,12,7],[0,2,3,6,13,14,15,16,17,18],-0.096100))
Atoms.append(Atom(8,'c3',0.421000,-1.380000,0.456000,np.array([0.421000,-1.380000,0.456000]),[6,16,17,18],[8,1,12,7],[0,2,3,6,13,14,15,16,17,18],-0.096100))
Atoms.append(Atom(9,'hn',4.144000,-0.703000,-0.953000,np.array([4.144000,-0.703000,-0.953000]),[0],[9,10,1],[0,2,3,6],0.363800))
Atoms.append(Atom(10,'hn',2.726000,-1.062000,-1.635000,np.array([2.726000,-1.062000,-1.635000]),[0],[9,10,1],[0,2,3,6],0.363800))
Atoms.append(Atom(11,'ho',3.644000,-1.761000,0.862000,np.array([3.644000,-1.761000,0.862000]),[4],[11,3],[1,4,5],0.444000))
Atoms.append(Atom(12,'hc',0.626000,0.713000,0.797000,np.array([0.626000,0.713000,0.797000]),[6],[8,1,12,7],[0,2,3,6,13,14,15,16,17,18],0.083700))
Atoms.append(Atom(13,'hc',0.711000,-0.318000,-2.086000,np.array([0.711000,-0.318000,-2.086000]),[7],[14,13,6,15],[8,1,12,7],0.039367))
Atoms.append(Atom(14,'hc',0.856000,1.382000,-1.599000,np.array([0.856000,1.382000,-1.599000]),[7],[14,13,6,15],[8,1,12,7],0.039367))
Atoms.append(Atom(15,'hc',-0.635000,0.489000,-1.275000,np.array([-0.635000,0.489000,-1.275000]),[7],[14,13,6,15],[8,1,12,7],0.039367))
Atoms.append(Atom(16,'hc',0.713000,-2.163000,-0.252000,np.array([0.713000,-2.163000,-0.252000]),[8],[16,17,18,6],[8,1,12,7],0.039367))
Atoms.append(Atom(17,'hc',-0.674000,-1.340000,0.472000,np.array([-0.674000,-1.340000,0.472000]),[8],[16,17,18,6],[8,1,12,7],0.039367))
Atoms.append(Atom(18,'hc',0.752000,-1.680000,1.455000,np.array([0.752000,-1.680000,1.455000]),[8],[16,17,18,6],[8,1,12,7],0.039367))
Bonds.append(Bond(0,'n3-hn',0,9,111))
Bonds.append(Bond(1,'n3-hn',0,10,111))
Bonds.append(Bond(2,'c3-c',1,3,111))
Bonds.append(Bond(3,'c3-h1',1,2,111))
Bonds.append(Bond(4,'c3-c3',1,6,111))
Bonds.append(Bond(5,'c3-n3',1,0,111))
Bonds.append(Bond(6,'c-o',3,5,111))
Bonds.append(Bond(7,'c-oh',3,4,111))
Bonds.append(Bond(8,'oh-ho',4,11,111))
Bonds.append(Bond(9,'c3-c3',6,8,111))
Bonds.append(Bond(10,'c3-c3',6,7,111))
Bonds.append(Bond(11,'c3-hc',6,12,111))
Bonds.append(Bond(12,'c3-hc',7,13,111))
Bonds.append(Bond(13,'c3-hc',7,14,111))
Bonds.append(Bond(14,'c3-hc',7,15,111))
Bonds.append(Bond(15,'c3-hc',8,16,111))
Bonds.append(Bond(16,'c3-hc',8,17,111))
Bonds.append(Bond(17,'c3-hc',8,18,111))
Angles.append(Angle(0,'hn-n3-hn',9,0,10,222))
Angles.append(Angle(1,'hn-n3-c3',9,0,1,222))
Angles.append(Angle(2,'hn-n3-c3',10,0,1,222))
Angles.append(Angle(3,'c-c3-h1',3,1,2,222))
Angles.append(Angle(4,'c-c3-c3',3,1,6,222))
Angles.append(Angle(5,'c-c3-n3',3,1,0,222))
Angles.append(Angle(6,'c3-c-o',1,3,5,222))
Angles.append(Angle(7,'c3-c-oh',1,3,4,222))
Angles.append(Angle(8,'h1-c3-c3',2,1,6,222))
Angles.append(Angle(9,'h1-c3-n3',2,1,0,222))
Angles.append(Angle(10,'c3-c3-n3',6,1,0,222))
Angles.append(Angle(11,'c3-c3-c3',1,6,8,222))
Angles.append(Angle(12,'c3-c3-c3',1,6,7,222))
Angles.append(Angle(13,'c3-c3-hc',1,6,12,222))
Angles.append(Angle(14,'o-c-oh',5,3,4,222))
Angles.append(Angle(15,'c-oh-ho',3,4,11,222))
Angles.append(Angle(16,'c3-c3-c3',8,6,7,222))
Angles.append(Angle(17,'c3-c3-hc',8,6,12,222))
Angles.append(Angle(18,'c3-c3-hc',6,8,16,222))
Angles.append(Angle(19,'c3-c3-hc',6,8,17,222))
Angles.append(Angle(20,'c3-c3-hc',6,8,18,222))
Angles.append(Angle(21,'c3-c3-hc',7,6,12,222))
Angles.append(Angle(22,'c3-c3-hc',6,7,13,222))
Angles.append(Angle(23,'c3-c3-hc',6,7,14,222))
Angles.append(Angle(24,'c3-c3-hc',6,7,15,222))
Angles.append(Angle(25,'hc-c3-hc',13,7,14,222))
Angles.append(Angle(26,'hc-c3-hc',13,7,15,222))
Angles.append(Angle(27,'hc-c3-hc',14,7,15,222))
Angles.append(Angle(28,'hc-c3-hc',16,8,17,222))
Angles.append(Angle(29,'hc-c3-hc',16,8,18,222))
Angles.append(Angle(30,'hc-c3-hc',17,8,18,222))
Torsions.append(Torsion(0,'hn-n3-c3-c',9,0,1,3,333))
Torsions.append(Torsion(1,'hn-n3-c3-h1',9,0,1,2,333))
Torsions.append(Torsion(2,'hn-n3-c3-c3',9,0,1,6,333))
Torsions.append(Torsion(3,'hn-n3-c3-c',10,0,1,3,333))
Torsions.append(Torsion(4,'hn-n3-c3-h1',10,0,1,2,333))
Torsions.append(Torsion(5,'hn-n3-c3-c3',10,0,1,6,333))
Torsions.append(Torsion(6,'o-c-c3-h1',5,3,1,2,333))
Torsions.append(Torsion(7,'oh-c-c3-h1',4,3,1,2,333))
Torsions.append(Torsion(8,'o-c-c3-c3',5,3,1,6,333))
Torsions.append(Torsion(9,'oh-c-c3-c3',4,3,1,6,333))
Torsions.append(Torsion(10,'c-c3-c3-c3',3,1,6,8,333))
Torsions.append(Torsion(11,'c-c3-c3-c3',3,1,6,7,333))
Torsions.append(Torsion(12,'c-c3-c3-hc',3,1,6,12,333))
Torsions.append(Torsion(13,'o-c-c3-n3',5,3,1,0,333))
Torsions.append(Torsion(14,'oh-c-c3-n3',4,3,1,0,333))
Torsions.append(Torsion(15,'c3-c-oh-ho',1,3,4,11,333))
Torsions.append(Torsion(16,'h1-c3-c3-c3',2,1,6,8,333))
Torsions.append(Torsion(17,'h1-c3-c3-c3',2,1,6,7,333))
Torsions.append(Torsion(18,'h1-c3-c3-hc',2,1,6,12,333))
Torsions.append(Torsion(19,'c3-c3-c3-n3',8,6,1,0,333))
Torsions.append(Torsion(20,'c3-c3-c3-n3',7,6,1,0,333))
Torsions.append(Torsion(21,'hc-c3-c3-n3',12,6,1,0,333))
Torsions.append(Torsion(22,'c-c3-c3-c3',3,1,6,8,333))
Torsions.append(Torsion(23,'h1-c3-c3-c3',2,1,6,8,333))
Torsions.append(Torsion(24,'c3-c3-c3-hc',1,6,8,16,333))
Torsions.append(Torsion(25,'c3-c3-c3-hc',1,6,8,17,333))
Torsions.append(Torsion(26,'c3-c3-c3-hc',1,6,8,18,333))
Torsions.append(Torsion(27,'c-c3-c3-c3',3,1,6,7,333))
Torsions.append(Torsion(28,'h1-c3-c3-c3',2,1,6,7,333))
Torsions.append(Torsion(29,'c3-c3-c3-hc',1,6,7,13,333))
Torsions.append(Torsion(30,'c3-c3-c3-hc',1,6,7,14,333))
Torsions.append(Torsion(31,'c3-c3-c3-hc',1,6,7,15,333))
Torsions.append(Torsion(32,'c-c3-c3-hc',3,1,6,12,333))
Torsions.append(Torsion(33,'h1-c3-c3-hc',2,1,6,12,333))
Torsions.append(Torsion(34,'o-c-oh-ho',5,3,4,11,333))
Torsions.append(Torsion(35,'c3-c-oh-ho',1,3,4,11,333))
Torsions.append(Torsion(36,'o-c-oh-ho',5,3,4,11,333))
Torsions.append(Torsion(37,'hc-c3-c3-c3',16,8,6,7,333))
Torsions.append(Torsion(38,'hc-c3-c3-c3',17,8,6,7,333))
Torsions.append(Torsion(39,'hc-c3-c3-c3',18,8,6,7,333))
Torsions.append(Torsion(40,'c3-c3-c3-hc',8,6,7,13,333))
Torsions.append(Torsion(41,'c3-c3-c3-hc',8,6,7,14,333))
Torsions.append(Torsion(42,'c3-c3-c3-hc',8,6,7,15,333))
Torsions.append(Torsion(43,'hc-c3-c3-hc',16,8,6,12,333))
Torsions.append(Torsion(44,'hc-c3-c3-hc',17,8,6,12,333))
Torsions.append(Torsion(45,'hc-c3-c3-hc',18,8,6,12,333))
Torsions.append(Torsion(46,'c3-c3-c3-hc',1,6,8,16,333))
Torsions.append(Torsion(47,'c3-c3-c3-hc',1,6,8,17,333))
Torsions.append(Torsion(48,'c3-c3-c3-hc',1,6,8,18,333))
Torsions.append(Torsion(49,'hc-c3-c3-hc',13,7,6,12,333))
Torsions.append(Torsion(50,'hc-c3-c3-hc',14,7,6,12,333))
Torsions.append(Torsion(51,'hc-c3-c3-hc',15,7,6,12,333))
Torsions.append(Torsion(52,'c3-c3-c3-hc',1,6,7,13,333))
Torsions.append(Torsion(53,'c3-c3-c3-hc',8,6,7,13,333))
Torsions.append(Torsion(54,'c3-c3-c3-hc',1,6,7,14,333))
Torsions.append(Torsion(55,'c3-c3-c3-hc',8,6,7,14,333))
Torsions.append(Torsion(56,'c3-c3-c3-hc',1,6,7,15,333))
Torsions.append(Torsion(57,'c3-c3-c3-hc',8,6,7,15,333))
inputdata = Sysdata(Atoms,Bonds,Angles,Torsions)