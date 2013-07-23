

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


Atoms.append(Atom(0,'oh',1.179000,-0.268000,0.277000,np.array([1.179000,-0.268000,0.277000]),[1,11],[0,10,2],[3,1,11,4,5],-0.599100))
Atoms.append(Atom(1,'c',2.436000,-0.119000,-0.179000,np.array([2.436000,-0.119000,-0.179000]),[0,10,2],[3,1,11,4,5],[0,2,6,10,12,13,14,15],0.584100))
Atoms.append(Atom(2,'c3',3.146000,1.159000,0.301000,np.array([3.146000,1.159000,0.301000]),[1,5,4,3],[0,2,6,10,12,13,14,15],[1,3,4,5,7,11,16,17],0.135500))
Atoms.append(Atom(3,'h1',3.682000,0.894000,1.220000,np.array([3.682000,0.894000,1.220000]),[2],[1,3,4,5],[0,2,6,10,12,13,14,15],0.097700))
Atoms.append(Atom(4,'n3',4.138000,1.569000,-0.716000,np.array([4.138000,1.569000,-0.716000]),[2,12,13],[1,3,4,5],[0,2,6,10,12,13,14,15],-0.896800))
Atoms.append(Atom(5,'c3',2.168000,2.306000,0.596000,np.array([2.168000,2.306000,0.596000]),[2,6,14,15],[1,3,4,5,7,16,17],[0,2,6,8,9,10,12,13,14,15],-0.069400))
Atoms.append(Atom(6,'c3',2.864000,3.577000,1.088000,np.array([2.864000,3.577000,1.088000]),[5,7,16,17],[2,6,8,9,14,15],[1,3,4,5,7,16,17,18,19],-0.150400))
Atoms.append(Atom(7,'c',1.960000,4.720000,1.511000,np.array([1.960000,4.720000,1.511000]),[6,9,8],[5,7,16,17,18,19],[2,6,8,9,14,15],0.655100))
Atoms.append(Atom(8,'n',0.615000,4.560000,1.414000,np.array([0.615000,4.560000,1.414000]),[7,18,19],[8,9,6],[5,7,16,17,18,19],-0.679000))
Atoms.append(Atom(9,'o',2.413000,5.781000,1.936000,np.array([2.413000,5.781000,1.936000]),[7],[8,9,6],[5,7,16,17,18,19],-0.610100))
Atoms.append(Atom(10,'o',2.933000,-0.982000,-0.888000,np.array([2.933000,-0.982000,-0.888000]),[1],[0,10,2],[3,1,11,4,5],-0.557000))
Atoms.append(Atom(11,'ho',0.891000,-1.125000,-0.106000,np.array([0.891000,-1.125000,-0.106000]),[0],[1,11],[0,10,2],0.447000))
Atoms.append(Atom(12,'hn',3.657000,1.764000,-1.595000,np.array([3.657000,1.764000,-1.595000]),[4],[2,12,13],[1,3,4,5],0.366300))
Atoms.append(Atom(13,'hn',4.728000,0.757000,-0.924000,np.array([4.728000,0.757000,-0.924000]),[4],[2,12,13],[1,3,4,5],0.366300))
Atoms.append(Atom(14,'hc',1.460000,1.974000,1.365000,np.array([1.460000,1.974000,1.365000]),[5],[2,14,6,15],[1,3,4,5,7,16,17],0.062700))
Atoms.append(Atom(15,'hc',1.572000,2.529000,-0.298000,np.array([1.572000,2.529000,-0.298000]),[5],[2,14,6,15],[1,3,4,5,7,16,17],0.062700))
Atoms.append(Atom(16,'hc',3.516000,3.974000,0.302000,np.array([3.516000,3.974000,0.302000]),[6],[16,17,5,7],[2,6,8,9,14,15],0.081200))
Atoms.append(Atom(17,'hc',3.490000,3.334000,1.955000,np.array([3.490000,3.334000,1.955000]),[6],[16,17,5,7],[2,6,8,9,14,15],0.081200))
Atoms.append(Atom(18,'hn',0.154000,3.733000,1.071000,np.array([0.154000,3.733000,1.071000]),[8],[18,19,7],[8,9,6],0.311000))
Atoms.append(Atom(19,'hn',0.039000,5.343000,1.703000,np.array([0.039000,5.343000,1.703000]),[8],[18,19,7],[8,9,6],0.311000))
Bonds.append(Bond(0,'oh-c',0,1,111))
Bonds.append(Bond(1,'oh-ho',0,11,111))
Bonds.append(Bond(2,'c-o',1,10,111))
Bonds.append(Bond(3,'c-c3',1,2,111))
Bonds.append(Bond(4,'c3-c3',2,5,111))
Bonds.append(Bond(5,'c3-n3',2,4,111))
Bonds.append(Bond(6,'c3-h1',2,3,111))
Bonds.append(Bond(7,'n3-hn',4,12,111))
Bonds.append(Bond(8,'n3-hn',4,13,111))
Bonds.append(Bond(9,'c3-c3',5,6,111))
Bonds.append(Bond(10,'c3-hc',5,14,111))
Bonds.append(Bond(11,'c3-hc',5,15,111))
Bonds.append(Bond(12,'c3-c',6,7,111))
Bonds.append(Bond(13,'c3-hc',6,16,111))
Bonds.append(Bond(14,'c3-hc',6,17,111))
Bonds.append(Bond(15,'c-o',7,9,111))
Bonds.append(Bond(16,'c-n',7,8,111))
Bonds.append(Bond(17,'n-hn',8,18,111))
Bonds.append(Bond(18,'n-hn',8,19,111))
Angles.append(Angle(0,'c-oh-ho',1,0,11,222))
Angles.append(Angle(1,'oh-c-o',0,1,10,222))
Angles.append(Angle(2,'oh-c-c3',0,1,2,222))
Angles.append(Angle(3,'o-c-c3',10,1,2,222))
Angles.append(Angle(4,'c-c3-c3',1,2,5,222))
Angles.append(Angle(5,'c-c3-n3',1,2,4,222))
Angles.append(Angle(6,'c-c3-h1',1,2,3,222))
Angles.append(Angle(7,'c3-c3-n3',5,2,4,222))
Angles.append(Angle(8,'c3-c3-h1',5,2,3,222))
Angles.append(Angle(9,'c3-c3-c3',2,5,6,222))
Angles.append(Angle(10,'c3-c3-hc',2,5,14,222))
Angles.append(Angle(11,'c3-c3-hc',2,5,15,222))
Angles.append(Angle(12,'n3-c3-h1',4,2,3,222))
Angles.append(Angle(13,'c3-n3-hn',2,4,12,222))
Angles.append(Angle(14,'c3-n3-hn',2,4,13,222))
Angles.append(Angle(15,'hn-n3-hn',12,4,13,222))
Angles.append(Angle(16,'c3-c3-hc',6,5,14,222))
Angles.append(Angle(17,'c3-c3-hc',6,5,15,222))
Angles.append(Angle(18,'c3-c3-c',5,6,7,222))
Angles.append(Angle(19,'c3-c3-hc',5,6,16,222))
Angles.append(Angle(20,'c3-c3-hc',5,6,17,222))
Angles.append(Angle(21,'hc-c3-hc',14,5,15,222))
Angles.append(Angle(22,'c-c3-hc',7,6,16,222))
Angles.append(Angle(23,'c-c3-hc',7,6,17,222))
Angles.append(Angle(24,'c3-c-o',6,7,9,222))
Angles.append(Angle(25,'c3-c-n',6,7,8,222))
Angles.append(Angle(26,'hc-c3-hc',16,6,17,222))
Angles.append(Angle(27,'o-c-n',9,7,8,222))
Angles.append(Angle(28,'c-n-hn',7,8,18,222))
Angles.append(Angle(29,'c-n-hn',7,8,19,222))
Angles.append(Angle(30,'hn-n-hn',18,8,19,222))
Torsions.append(Torsion(0,'o-c-oh-ho',10,1,0,11,333))
Torsions.append(Torsion(1,'c3-c-oh-ho',2,1,0,11,333))
Torsions.append(Torsion(2,'oh-c-c3-c3',0,1,2,5,333))
Torsions.append(Torsion(3,'oh-c-c3-n3',0,1,2,4,333))
Torsions.append(Torsion(4,'oh-c-c3-h1',0,1,2,3,333))
Torsions.append(Torsion(5,'o-c-c3-c3',10,1,2,5,333))
Torsions.append(Torsion(6,'o-c-c3-n3',10,1,2,4,333))
Torsions.append(Torsion(7,'o-c-c3-h1',10,1,2,3,333))
Torsions.append(Torsion(8,'oh-c-c3-c3',0,1,2,5,333))
Torsions.append(Torsion(9,'o-c-c3-c3',10,1,2,5,333))
Torsions.append(Torsion(10,'c-c3-c3-c3',1,2,5,6,333))
Torsions.append(Torsion(11,'c-c3-c3-hc',1,2,5,14,333))
Torsions.append(Torsion(12,'c-c3-c3-hc',1,2,5,15,333))
Torsions.append(Torsion(13,'oh-c-c3-n3',0,1,2,4,333))
Torsions.append(Torsion(14,'o-c-c3-n3',10,1,2,4,333))
Torsions.append(Torsion(15,'c-c3-n3-hn',1,2,4,12,333))
Torsions.append(Torsion(16,'c-c3-n3-hn',1,2,4,13,333))
Torsions.append(Torsion(17,'oh-c-c3-h1',0,1,2,3,333))
Torsions.append(Torsion(18,'o-c-c3-h1',10,1,2,3,333))
Torsions.append(Torsion(19,'c3-c3-c3-n3',6,5,2,4,333))
Torsions.append(Torsion(20,'hc-c3-c3-n3',14,5,2,4,333))
Torsions.append(Torsion(21,'hc-c3-c3-n3',15,5,2,4,333))
Torsions.append(Torsion(22,'c3-c3-n3-hn',5,2,4,12,333))
Torsions.append(Torsion(23,'c3-c3-n3-hn',5,2,4,13,333))
Torsions.append(Torsion(24,'c3-c3-c3-h1',6,5,2,3,333))
Torsions.append(Torsion(25,'hc-c3-c3-h1',14,5,2,3,333))
Torsions.append(Torsion(26,'hc-c3-c3-h1',15,5,2,3,333))
Torsions.append(Torsion(27,'c-c3-c3-c3',1,2,5,6,333))
Torsions.append(Torsion(28,'c3-c3-c3-c',2,5,6,7,333))
Torsions.append(Torsion(29,'c3-c3-c3-hc',2,5,6,16,333))
Torsions.append(Torsion(30,'c3-c3-c3-hc',2,5,6,17,333))
Torsions.append(Torsion(31,'c-c3-c3-hc',1,2,5,14,333))
Torsions.append(Torsion(32,'c-c3-c3-hc',1,2,5,15,333))
Torsions.append(Torsion(33,'hn-n3-c3-h1',12,4,2,3,333))
Torsions.append(Torsion(34,'hn-n3-c3-h1',13,4,2,3,333))
Torsions.append(Torsion(35,'c-c3-n3-hn',1,2,4,12,333))
Torsions.append(Torsion(36,'c3-c3-n3-hn',5,2,4,12,333))
Torsions.append(Torsion(37,'c-c3-n3-hn',1,2,4,13,333))
Torsions.append(Torsion(38,'c3-c3-n3-hn',5,2,4,13,333))
Torsions.append(Torsion(39,'c-c3-c3-hc',7,6,5,14,333))
Torsions.append(Torsion(40,'hc-c3-c3-hc',16,6,5,14,333))
Torsions.append(Torsion(41,'hc-c3-c3-hc',17,6,5,14,333))
Torsions.append(Torsion(42,'c-c3-c3-hc',7,6,5,15,333))
Torsions.append(Torsion(43,'hc-c3-c3-hc',16,6,5,15,333))
Torsions.append(Torsion(44,'hc-c3-c3-hc',17,6,5,15,333))
Torsions.append(Torsion(45,'c3-c3-c3-c',2,5,6,7,333))
Torsions.append(Torsion(46,'c3-c3-c-o',5,6,7,9,333))
Torsions.append(Torsion(47,'c3-c3-c-n',5,6,7,8,333))
Torsions.append(Torsion(48,'c3-c3-c3-hc',2,5,6,16,333))
Torsions.append(Torsion(49,'c3-c3-c3-hc',2,5,6,17,333))
Torsions.append(Torsion(50,'o-c-c3-hc',9,7,6,16,333))
Torsions.append(Torsion(51,'n-c-c3-hc',8,7,6,16,333))
Torsions.append(Torsion(52,'o-c-c3-hc',9,7,6,17,333))
Torsions.append(Torsion(53,'n-c-c3-hc',8,7,6,17,333))
Torsions.append(Torsion(54,'c3-c3-c-o',5,6,7,9,333))
Torsions.append(Torsion(55,'c3-c3-c-n',5,6,7,8,333))
Torsions.append(Torsion(56,'c3-c-n-hn',6,7,8,18,333))
Torsions.append(Torsion(57,'c3-c-n-hn',6,7,8,19,333))
Torsions.append(Torsion(58,'o-c-n-hn',9,7,8,18,333))
Torsions.append(Torsion(59,'o-c-n-hn',9,7,8,19,333))
Torsions.append(Torsion(60,'c3-c-n-hn',6,7,8,18,333))
Torsions.append(Torsion(61,'o-c-n-hn',9,7,8,18,333))
Torsions.append(Torsion(62,'c3-c-n-hn',6,7,8,19,333))
Torsions.append(Torsion(63,'o-c-n-hn',9,7,8,19,333))
inputdata = Sysdata(Atoms,Bonds,Angles,Torsions)