#!/usr/bin/env python
#-*- coding: utf-8 -*-

'''
Craeted at 2011-01-14.
@author: zhuh
@version: 0.1.0
@change:\n
    - 2011-01-14.\n
        - finished the coor data for the five base.
    - 2011-01-18.\n
        - add the base list BASE_AG_LIST and BASE_CTU_LIST.
'''
import numpy as np

BASE_AG_LIST=["N9","C8","N7","C5","C6","N1","C2", "N3","C4"]
'''It's a list of atoms which used to fit the standard structure for the A and G base.'''
BASE_CTU_LIST=["N1","C2","N3","C4","C5","C6"]
'''It's a list of atoms which used to fit the standard structure for the C,T and U base.'''

RESIDUE_NAME_LIST=['A','T','C','G',\
        'DA','DT','DC','DG',\
        'DA5','DT5','DC5','DG5',\
        'DA3','DT3','DC3','DG3',\
        'RA','RT','RC','RG',\
        'RA5','RT5','RC5','RG5',\
        'RA3','RT3','RC3','RG3']
'''
just for nucleic now. pretein residue can be append if need.
'''

BASE_A_DICT={
        "N9":(-1.291, 4.498, 0.000),
        "C8":(0.024, 4.897, 0.000),
        "N7":(0.877, 3.902, 0.000),
        "C5":(0.071, 2.771, 0.000),
        "C6":(0.369, 1.398, 0.000),
        "N1":(-0.668, 0.532, 0.000),
        "C2":(-1.912, 1.023, 0.000),
        "N3":(-2.320, 2.290, 0.000),
        "C4":(-1.267, 3.124, 0.000) 
        }
'''
The cartesian coordinates of non-hydrogen atoms in the standard reference frames \
        of the Adenine base.\
        Note: the unit is A.
        '''


BASE_C_DICT={
        "N1":(-1.285,4.542,0.000),
        "C2":(-1.472,3.158,0.000),
        "N3":(-0.391,2.344,0.000),
        "C4":(0.837,2.868,0.000),
        "C5":(1.056,4.275,0.000),
        "C6":(-0.023,5.068,0.000)
        }
'''The cartesian coordinates of non-hydrogen atoms in the standard reference frames \
        of the Cytosine base.'''

BASE_G_DICT={
        "N9":(-1.289,4.551,0.000),
        "C8":( 0.023,4.962,0.000),
        "N7":( 0.870,3.969,0.000),
        "C5":( 0.071,2.833,0.000),
        "C6":( 0.424,1.460,0.000),
        "N1":(-0.700,0.641,0.000),
        "C2":(-1.999,1.087,0.000),
        "N3":(-2.342,2.364,0.001),
        "C4":(-1.265,3.177,0.000)
        }
'''The cartesian coordinates of non-hydrogen atoms in the standard reference frames \
        of the Guanine base.'''

BASE_T_DICT={
        "N1":(-1.284,4.500,0.000),
        "C2":(-1.462,3.135,0.000),
        "N3":(-0.298,2.407,0.000),
        "C4":(0.994,2.897,0.000),
        "C5":(1.106,4.338,0.000),
        "C6":(-0.024,5.057,0.000)
        }
'''The cartesian coordinates of non-hydrogen atoms in the standard reference frames \
        of the Thymine base.'''

BASE_U_DICT={
        "N1":(-1.284,4.500,0.000),
        "C2":(-1.462,3.131,0.000),
        "N3":(-0.302,2.397,0.000),
        "C4":(0.989,2.884,0.000),
        "C5":(1.089,4.311,0.000),
        "C6":(-0.024,5.053,0.000)
        }
'''The cartesian coordinates of non-hydrogen atoms in the standard reference frames \
        of the Uracil base.'''


BASE_A_array=np.array([
            [-1.291,   4.498,   0.000 ],
            [ 0.024,   4.897,   0.000 ],
            [ 0.877,   3.902,   0.000 ],
            [ 0.071,   2.771,   0.000 ],
            [ 0.369,   1.398,   0.000 ],
            [-0.668,   0.532,   0.000 ],
            [-1.912,   1.023,   0.000 ],
            [-2.320,   2.290,   0.000 ],
            [-1.267,   3.124,   0.000 ]
            ])
'''The cartesian coordinates of non-hydrogen atoms in the standard reference frames \
        of the Adenine base.'''


BASE_C_array=np.array([
           [-1.285,  4.542,  0.000],
           [-1.472,  3.158,  0.000],
           [-0.391,  2.344,  0.000], 
           [ 0.837,  2.868,  0.000], 
           [ 1.056,  4.275,  0.000], 
           [-0.023,  5.068,  0.000]
           ])
'''The cartesian coordinates of non-hydrogen atoms in the standard reference frames \
        of the Cytosine base.'''

BASE_G_array=np.array([
           [-1.289,  4.551,  0.000], 
           [ 0.023,  4.962,  0.000],  
           [ 0.870,  3.969,  0.000], 
           [ 0.071,  2.833,  0.000], 
           [ 0.424,  1.460,  0.000],
           [-0.700,  0.641,  0.000],
           [-1.999,  1.087,  0.000],
           [-2.342,  2.364,  0.001], 
           [-1.265,  3.177,  0.000]
           ])
'''The cartesian coordinates of non-hydrogen atoms in the standard reference frames \
        of the Guanine base.'''

BASE_T_array=np.array([
           [-1.284,  4.500,  0.000],
           [-1.462,  3.135,  0.000], 
           [-0.298,  2.407,  0.000],
           [ 0.994,  2.897,  0.000],
           [ 1.106,  4.338,  0.000], 
           [-0.024,  5.057,  0.000]
           ])
'''The cartesian coordinates of non-hydrogen atoms in the standard reference frames \
        of the Thymine base.'''

BASE_U_array=np.array([
           [-1.284,  4.500,  0.000], 
           [-1.462,  3.131,  0.000], 
           [-0.302,  2.397,  0.000], 
           [ 0.989,  2.884,  0.000], 
           [ 1.089,  4.311,  0.000],
           [-0.024,  5.053,  0.000]
           ])
'''The cartesian coordinates of non-hydrogen atoms in the standard reference frames \
        of the Uracil base.'''

