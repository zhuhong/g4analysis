from numpy import array
from numpy import matrix
from numpy import linalg
import math
import numpy


def Fitting(standard_coor,experim_coor):

    sum_standard = array([ sum(array(standard_coor)[:, 0]),  sum(array(standard_coor)[:, 1]),  sum(array(standard_coor)[:, 2]) ])
#   print sum_standard/9
    sum_experim  = array([ sum(array(experim_coor)[:, 0]),   sum(array(experim_coor)[:, 1]),   sum(array(experim_coor)[:, 2]) ])
#    print sum_experim/9
    NUM=len(standard_coor)

    temp_matrix = array([
        [sum_standard[0]*sum_experim[0], sum_standard[0]*sum_experim[1], sum_standard[0]*sum_experim[2]], 
        [sum_standard[1]*sum_experim[0], sum_standard[1]*sum_experim[1], sum_standard[1]*sum_experim[2]], 
        [sum_standard[2]*sum_experim[0], sum_standard[2]*sum_experim[1], sum_standard[2]*sum_experim[2]]
        ])  
#    print temp_matrix
#    print base_name
    co = (array(matrix(standard_coor).T*matrix(experim_coor)) - temp_matrix/NUM ) /(NUM-1) #covar_matrix
#    print co
    symm_matrix = array([
        [co[0, 0]+co[1, 1]+co[2, 2], co[1, 2]-co[2, 1],          co[2, 0]-co[0, 2],           co[0, 1]-co[1, 0]], 
        [co[1, 2]-co[2, 1],          co[0, 0]-co[1, 1]-co[2, 2], co[0, 1]+co[1, 0],           co[2, 0]+co[0, 2]], 
        [co[2, 0]-co[0, 2],          co[0, 1]+co[1, 0],         -co[0, 0]+co[1, 1]-co[2, 2],  co[1, 2]+co[2, 1]], 
        [co[0, 1]-co[1, 0],          co[2, 0]+co[0, 2],          co[1, 2]+co[2, 1],          -co[0, 0]-co[1, 1]+co[2, 2]] 
            ])  
#    print symm_matrix
    temp_list = linalg.eig(symm_matrix)
    '''temp_list contain all the eigenvalues and eigenvectors.temp_list[0] is eigenvalues, and temp_list[1] is eigenvectors.'''
 #   print temp_list
    temp_num = list(temp_list[0]).index(sorted(temp_list[0])[-1])
    '''temp_num is the index of the largest eigenvalue.'''
    q = temp_list[1][:, temp_num]
    '''the eigenvector cooresponding to the largest eigenvalue.'''

    R = array([
        [ q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2,  2*(q[1]*q[2] - q[0]*q[3]),              2*(q[1]*q[3]+q[0]*q[2])],
        [ 2*(q[2]*q[1] + q[0]*q[3]),              q[0]**2 - q[1]**2 + q[2]**2 - q[3]**2,  2*(q[2]*q[3]-q[0]*q[1])],
        [ 2*(q[3]*q[1] - q[0]*q[2]),              2*(q[3]*q[2] + q[0]*q[1]),              q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2]
        ])
    '''the rotation matrix'''

    origin_coor = array(sum_experim/NUM -matrix(sum_standard/NUM)*(matrix(R).T))
    '''the origin coordinate of the base'''

    result = [R, origin_coor[0]]
    return result




