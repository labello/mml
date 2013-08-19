
import math
import numpy as np

def vec_distance(v1,v2):
    ''' 
* Usage 
    v1 = np.array([0,0,0])
    v2 = np.array([1,1,1])
    r = vec_distance(v1,v2) 
* Input
    v1 - 1D Numpy array
    v2 - 1D Numpy array
    Both vectors must have the same number of elements 
* Output
    Floating point number which is the Cartesian distance between the two points 
* Description
    '''
    return np.linalg.norm(v2-v1)

def vec_angle(v1,v2,v3):
    '''
* Usage
    v1 = np.array([0,0,0])
    v2 = np.array([1,1,1])  # v2 is vertex of angle
    v3 = np.array([2,2,2])
    theta = vec_angle(v1,v2,v3)
* Input
    v1,v2,v3 - 1D Numpy array
    Vectors must have the same number of elements
* Output
    Floating point number which is the angle in degrees
* Description
    v1,v2,v3 refers to the vector of Cartesian coordinates describing each point
    vec1, vec2, etc... refers to the vector between 2 points, e.g, vec1 = v1-v2
    '''
    a = v1-v2
    b = v3-v2
    angle = np.arccos((a*b).sum(-1) / ( (a**2).sum(-1) * (b**2).sum(-1) )**0.5)
    return math.degrees(angle)

def vec_dihedral(v1,v2,v3,v4):
    '''
* Usage
    v1 = np.array([0,0,0])
    v2 = np.array([1,1,1])  # 
    v3 = np.array([2,2,2])  # The v2--v3 bond is rotated around
    v4 = np.array([3,3,3])
    phi = vec_dihedral(v1,v2,v3)
* Input
    v1,v2,v3 - 1D Numpy array
    Vectors must have the same number of elements
* Output
    Floating point number which is the angle in degrees
* Description
    v1,v2,v3,v4 refers to the vector of Cartesian coordinates describing each point
    vec1, vec2, etc... refers to the vector between 2 points, e.g, vec1 = v1-v2
    '''
    a1 = v2-v1
    a2 = v3-v2
    a3 = v4-v3

    v1 = np.cross(a1, a2)
    v2 = np.cross(a2, a3)
    
    v1 = v1 / (v1 * v1).sum()**0.5
    v2 = v2 / (v2 * v2).sum()**0.5

    porm = np.sign((v1 * a3).sum())
    rad = np.arccos((v1*v2).sum() / ((v1**2).sum() * (v2**2).sum())**0.5)
    angle = porm * rad
    if np.isnan(angle):
        angle = 1.0
    return math.degrees(angle)
    
def distance(p1,p2):
    '''
* Usage
    distance will first look x,y, and z in class or named tuple notation, e.g,
    p1.x,p1.y,p1.z = 0,0,0
    p2.x,p2.y,p2.z = 1,1,1
    r = distance(p1,p2)

    if that fails it will look for indexes, 0/1/2 assumed to be x/y/z.
 
    p1 = [0,0,0]
    p2 = [1,1,1]
* Input
    p1 - object or named tuple with .x, .y, and .z properties
    p2 - object or named tuple with .x, .y, and .z properties
* Output
    Floating point number which is the Cartesian distance between the two points
* Description
    '''
    try:
        return math.sqrt( (p2.x - p1.x)**2 + (p2.y - p1.y)**2 + (p2.z - p1.z)**2 )
    except AttributeError:
        return math.sqrt( (p2[0] - p1[0])**2 + (p2[1] - p1[1])**2 + (p2[2] - p1[2])**2 )
        

