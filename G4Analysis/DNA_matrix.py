#1/usr/bin/env python
#-*- coding:utf-8 -*-
'''
Created at 2012.02.23
@change:\n
    - 2012.02.23.
        - Copied those functions from parallel_analysis.py to make it simplify.
'''

from Coor import atomlib
import numpy
import Simple_atom
import time as Time
import math
import sys
import os
import usage
from mymath import least_squares_fitting
from math import cos, sin, sqrt
from numpy import matrix, linalg, dot

def Get_rotate_matrix(experim_coor, base_name):
    '''
    Reading in the coor_list data and base name,  return a 3*3 rotation matrix.\n
    The algorithm come form the 3DNA. you can find it from the 3DNA manual. The\
            details are in page 22-23.\n
    note: the unit is A.
    @param experim_coor: The coordinate for the 9 atoms for A G nucleic acid base or \
            6 atoms for C T U nucleic acid base.\n
    @param base_name: Only A,T,C,G,U are allowd.
    @type  experim_coor: numpy.array
    @type  base_name:    char
    @rtype:  list
    @return:  rotation matrix in list[0] and origin coordinate in list[1]
    '''

    if 'A' in base_name :
        standard_coor = atomlib.BASE_A_array 
    elif 'T' in base_name :
        standard_coor = atomlib.BASE_T_array 
    elif 'C' in base_name:
        standard_coor = atomlib.BASE_C_array 
    elif 'G' in base_name :
        standard_coor = atomlib.BASE_G_array 
    elif 'U' in base_name :
        standard_coor = atomlib.BASE_U_array 
    else:
        print "The residue name %s not found in my standard lib." %base_name
        sys.exit()
    
    result=least_squares_fitting.Fitting(standard_coor,experim_coor)
    return result

def Get_group_rotmat(group_rotate_list,group_size):
    '''
    get the rotation_matrix and origin_coor for a group which conatin 1,2 or 4 bases.
    '''
    z_axis=[]
    '''the z_axis for the base group'''
    origin = []
    '''the origin coordinate fot the base group'''
    if group_size  ==1:
        z_axis = group_rotate_list[0][0][:, 2]
        origin = group_rotate_list[0][1]
    elif group_size  ==2:
        origin = [(group_rotate_list[0][1][i]+group_rotate_list[1][1][i])/2 for i in range(3)]
        z_axis = Rotate_2_vector(group_rotate_list[0][0][:, 2], group_rotate_list[1][0][:, 2])         

    elif group_size  == 4:
        origin=[(group_rotate_list[0][1][i] + group_rotate_list[1][1][i] \
                + group_rotate_list[2][1][i] + group_rotate_list[3][1][i])/4 for i \
                in range(3)]

        temp1 = Rotate_2_vector(group_rotate_list[0][0][:, 2], group_rotate_list[1][0][:, 2])         
        temp2 = Rotate_2_vector(group_rotate_list[2][0][:, 2], group_rotate_list[3][0][:, 2])         
        z_axis = Rotate_2_vector(temp1, temp2)

    else:
        pass

    return z_axis,origin

def Get_baseID_list(atom_list,base_serial):

    atom_list = Simple_atom.Get_Atom_in_residue(atom_list, base_serial)
    temp_atom_list = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    for atom in atom_list:
        if ("G" in atom.residue_name ) or ("A" in atom.residue_name):
            if atom.atom_name in atomlib.BASE_AG_LIST:
                temp_atom_list[atomlib.BASE_AG_LIST.index(atom.atom_name)] = atom.atom_serial
        elif ("C" in atom.residue_name ) or ("T" in atom.residue_name) or ("U" in atom.residue_name):
            if atom.atom_name in atomlib.BASE_CTU_LIST:
                temp_atom_list[atomlib.BASE_CTU_LIST.index(atom.atom_name)] = atom.atom_serial
        else:
            print "the base", atom.residue_name, " not supposed."
    if temp_atom_list[6]==0:
        temp_atom_list = temp_atom_list[:6]
    return temp_atom_list

def Rotate_matrix(ga,cv):
    '''
     It's a general rotation matrix describe a rotation of mafnitude gamma about an arbitrary unit vector(cv1,cv2,cv3)
    '''
    rm = numpy.array([
    [ cos(ga)+(1- cos(ga))*cv[0]**2,  (1- cos(ga))*cv[0]*cv[1]-cv[2]* sin(ga),  (1- cos(ga))*cv[0]*cv[2]+cv[1]* sin(ga) ], 
    [ (1- cos(ga))*cv[0]*cv[1]+cv[2]* sin(ga),  cos(ga)+(1- cos(ga))*cv[1]**2,  (1- cos(ga))*cv[1]*cv[2]-cv[0]* sin(ga) ], 
    [ (1- cos(ga))*cv[0]*cv[2]-cv[1]* sin(ga),  (1- cos(ga))*cv[1]*cv[2] + cv[0]* sin(ga),  cos(ga)+(1- cos(ga))*cv[2]**2 ]
    ])
    return rm

def Norm_matrix_in_row(matrix):
    for i in range(3):
        temp_vect=numpy.array([matrix[j,i] for j in range(3)])
        temp_vect=temp_vect/math.sqrt(numpy.dot(temp_vect,temp_vect.T))
#        print numpy.array(temp_vect)[0]
        for j in range(3):
            matrix[j,i]=numpy.array(temp_vect)[j]
         #normalization the row.
    return matrix



def Rotate_2_vector(vector1, vector2):
    '''
    Get the average vector between the two vectors named vector1 and vector2.
    It's a standard method which can be found from the 3DNA manual. 
    '''
    vector1 = numpy.array(vector1)
    vector2 = numpy.array(vector2)
    if numpy.dot(vector1, vector2)<0:
        vector2 = vector2*(-1)
    ga = math.acos(numpy.dot(vector1, vector2))
    '''the half angle between the two vector'''
    cross_vector = numpy.cross(vector1, vector2)
    cv = cross_vector/math.sqrt(numpy.dot(cross_vector, cross_vector))
    '''The unit cross vector.'''
    rm = Rotate_matrix(-ga/2,cv)
    result = numpy.matrix(rm)*(numpy.matrix(vector2).T)
    result = result/math.sqrt(result.T*result)

    return numpy.array(result.T)[0]

def Rotate_2_matrix(mat1,mat2):
    '''
    Get the average vector between the two vectors named vector1 and vector2.
    It's a standard method which can be found from the 3DNA manual. 
    '''
    vector1 = numpy.array(mat1[:,2])
    vector2 = numpy.array(mat2[:,2])
    if numpy.dot(vector1, vector2)<0:
        vector2 = vector2*(-1)
    ga = math.acos(numpy.dot(vector1, vector2))
    '''the half angle between the two vector'''
    cross_vector = numpy.cross(vector1, vector2)
    cv = cross_vector/math.sqrt(numpy.dot(cross_vector, cross_vector))
    '''The unit cross vector.'''
    rm1 = Rotate_matrix(ga/2,cv) 
    '''The rotation matrix for matrix1'''
    rm2 = Rotate_matrix(-ga/2,cv) 
    '''The rotation matrix for matrix2'''

    result1 = numpy.matrix(rm1)*(numpy.matrix(mat1))
    result2 = numpy.matrix(rm2)*(numpy.matrix(mat2))
    return result1,result2


def Get_RMSD(experim_coor,base_name, origin, orient):
    N = 0
    if 'A' in base_name :
        standard_coor = atomlib.BASE_A_array 
        N = 9
    elif 'T' in base_name :
        standard_coor = atomlib.BASE_T_array 
        N = 6
    elif 'C' in base_name:
        standard_coor = atomlib.BASE_C_array 
        N = 6
    elif 'G' in base_name :
        standard_coor = atomlib.BASE_G_array 
        N = 9
    elif 'U' in base_name :
        standard_coor = atomlib.BASE_U_array 
        N = 6
    else:
        pass
    

    Fit_coor=numpy.matrix(standard_coor) * numpy.matrix(orient).T + origin
#    print Fit_coor
#    print experim_coor
    Diff_coor=experim_coor - Fit_coor

    ssum=0.0
    for i in range(N):
        for j in range(3):
            ssum = ssum + Diff_coor[i,j]*Diff_coor[i,j]
            
    RMSD=math.sqrt(ssum/N)
    return RMSD

def Get_Z_RMSD(experim_coor,origin,vect_z):
    '''
    get the rmsd in the z axis for this base.
    '''
    RMSD=0.0
    
    for i_coor in experim_coor:
        temp_vector=[i_coor[i]-origin[i] for i in range(3)]
        '''the vector from origin to base atom i'''
        temp_dist=numpy.dot(numpy.array(temp_vector),numpy.array(vect_z))
        RMSD=RMSD+temp_dist*temp_dist

    return math.sqrt(RMSD/len(experim_coor))


def Get_group_RMSD(base_name_list,expr_coor_list,origin_group,orient_group):
    '''
    get the rmsd for a bases group, if the group contain only one base. Using Get_RMSD to calculate the
    RMSD. or the group contain 2 or 4 bases. the Get_Z_RMSD is used to calculate the RMSD, which only 
    calculate the RMSD in z_axis.
    '''
    RMSD=0.0
    base_num=len(base_name_list)

    if base_num ==1:
        RMSD=Get_Z_RMSD(expr_coor_list[0],origin_group,orient_group)

    elif base_num == 2 or base_num==4:
        N=0
        Sum=0.0

        for i in range(base_num):
            r_i=Get_Z_RMSD(expr_coor_list[i],origin_group,orient_group)
            N_i=len(expr_coor_list[i])
            N=N+N_i
            Sum=Sum+r_i*r_i*N_i
        RMSD=math.sqrt(Sum/N)

    else:
        pass

    return RMSD

