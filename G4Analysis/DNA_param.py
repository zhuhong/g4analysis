import numpy
import math
import DNA_matrix

def base_pair_parameters(rotation_1,rotation_2,origin_1,origin_2):
    z1_v=numpy.array([rotation_1[i][2] for i in range(3)])
    z2_v=numpy.array([rotation_2[i][2] for i in range(3)])

    if numpy.dot(z1_v,z2_v.T) < 0:
        for i in range(3):
            for j in range(1,3):
                rotation_2[i][j]=-rotation_2[i][j]

    y1_vector=numpy.array([rotation_1[i][1] for i in range(3)])
    y2_vector=numpy.array([rotation_2[i][1] for i in range(3)])

    gamma=math.acos(numpy.dot(y1_vector,y2_vector))  

    #calculate buckleopening angle.
    b0_vector=numpy.cross(y2_vector,y1_vector) 
    #buckle-opening axis
    b0_vector=b0_vector/math.sqrt(numpy.dot(b0_vector,b0_vector.T))
    # normalize the b0 

    rotate_matrix_1=DNA_matrix.Rotate_matrix(-gamma/2,b0_vector)
    rotate_matrix_2=DNA_matrix.Rotate_matrix(gamma/2,b0_vector)
    #create rotate_matrix_2 

    trans_orient_1=numpy.matrix(rotate_matrix_1)*numpy.matrix(rotation_1)
    trans_orient_2=numpy.matrix(rotate_matrix_2)*numpy.matrix(rotation_2)
    #get the orientation matirx of transformed bases.

    MBT_matirx=(trans_orient_1+trans_orient_2)/2 
#//        MBT_origin[i]=(origin_1+origin_2)/2 
   
    MBT_matirx=DNA_matrix.Norm_matrix_in_row(MBT_matirx)
   #normalization the row.

    x_vector_temp_1=numpy.array([trans_orient_1[i,0] for i in range(3)])
    x_vector_temp_2=numpy.array([trans_orient_2[i,0] for i in range(3)])
    x_vector_temp_1=x_vector_temp_1/math.sqrt(numpy.dot(x_vector_temp_1,x_vector_temp_1.T))
    x_vector_temp_2=x_vector_temp_2/math.sqrt(numpy.dot(x_vector_temp_2,x_vector_temp_2.T))
#    print x_vector_temp_1
#    print x_vector_temp_2
    y_vector_MBT=numpy.array([MBT_matirx[i,1] for i in range(3)])

    omega=math.acos(numpy.dot(x_vector_temp_1,x_vector_temp_2.T)) 

    cross_result_temp=numpy.cross(x_vector_temp_2,x_vector_temp_1) 

    if(numpy.dot(cross_result_temp,y_vector_MBT.T)<0):
        omega=-abs(omega)
    else:
        omega=abs(omega)

    x_vector_MBT=numpy.array([MBT_matirx[i,0] for i in range(3)])

    phi=math.acos(numpy.dot(b0_vector,x_vector_MBT.T)) 

    cross_b0_xmat_temp=numpy.cross(b0_vector,x_vector_MBT) 

    if(numpy.dot(cross_b0_xmat_temp,y_vector_MBT)<0):
        phi= - abs(phi)
    else:
        phi=abs(phi)

    kappa=gamma * math.cos(phi) 
    sigma=gamma * math.sin(phi) 

    #get kappa and sigma 

    displacement=numpy.matrix(numpy.array(origin_1)-numpy.array(origin_2))*MBT_matirx

    shear=displacement[0,0] 
    stretch=displacement[0,1] 
    stagger=displacement[0,2] 
    buckle=kappa/math.pi * 180  
    propeller=omega /math.pi * 180 
    opening=sigma/math.pi * 180 

    return shear,stretch,stagger,buckle, propeller,opening

def base_step_parameters(rotation_1,rotation_2,origin_1,origin_2):

    z1_vector=numpy.array([rotation_1[i][2] for i in range(3)])
    z2_vector=numpy.array([rotation_2[i][2] for i in range(3)])

    gamma=math.acos(numpy.dot(z1_vector,z2_vector.T)) 
#    print gamma
    #   //calculate buckleopening angle.
    rt_vector=numpy.cross(z1_vector,z2_vector) 
    #//buckle-opening axis
    
    rt_vector=rt_vector/math.sqrt(numpy.dot(rt_vector,rt_vector.T))
    #   //   normalize the rt_vector b0 

    rotate_matrix_1=DNA_matrix.Rotate_matrix(gamma/2,rt_vector)
    rotate_matrix_2=DNA_matrix.Rotate_matrix(-gamma/2,rt_vector)
    #//create rotate_matrix_1, create rotate_matrix_2 

    trans_orient_1=numpy.matrix(rotate_matrix_1)*numpy.matrix(rotation_1)
    trans_orient_2=numpy.matrix(rotate_matrix_2)*numpy.matrix(rotation_2) 
#    print trans_orient_1
#    print trans_orient_2
    #//get the orientation matirx of transformed bases.

    MST_matirx=(trans_orient_1+trans_orient_2)/2 
    MST_matirx=DNA_matrix.Norm_matrix_in_row(MST_matirx)
    #normalization the row.

    y_vector_temp_1=numpy.array([trans_orient_1[i,1] for i in range(3)])
    y_vector_temp_2=numpy.array([trans_orient_2[i,1] for i in range(3)])
    y_vector_temp_1=y_vector_temp_1/math.sqrt(numpy.dot(y_vector_temp_1,y_vector_temp_1.T))
    y_vector_temp_2=y_vector_temp_2/math.sqrt(numpy.dot(y_vector_temp_2,y_vector_temp_2.T))
    z_vector_MST   =numpy.array([MST_matirx[i,2] for i in range(3)])
    y_vector_MST   =numpy.array([MST_matirx[i,1] for i in range(3)])

    omega=math.acos(numpy.dot(y_vector_temp_1,y_vector_temp_2.T)) 
    cross_result_temp=numpy.cross(y_vector_temp_1,y_vector_temp_2) 

    if(numpy.dot(cross_result_temp,z_vector_MST)<0):
        omega=-abs(omega) 
    else:
        omega=abs(omega)
    # // for omega

    phi=math.acos(numpy.dot(rt_vector,y_vector_MST.T)) 

    cross_rt_ymat_temp=numpy.cross(rt_vector,y_vector_MST)

    if(numpy.dot(cross_rt_ymat_temp,z_vector_MST)<0):
        phi= - abs(phi) 
    else:
        phi=abs(phi)
    #//for phi 

    kappa=gamma * math.cos(phi) 
    sigma=gamma * math.sin(phi) 

    #//get kappa and sigma 

    displacement=numpy.matrix(origin_2 - origin_1)*MST_matirx 

    shift=displacement[0,0] 
    slide=displacement[0,1] 
    rise=displacement[0,2] 
    roll=kappa/math.pi * 180  
    twist=omega /math.pi * 180 
    tilt=sigma/math.pi * 180 

    return shift,slide,rise,tilt,roll,twist


def middle_frame(rotation_1,rotation_2,origin_1,origin_2):

    z1_v=numpy.array([rotation_1[i][2] for i in range(3)])
    z2_v=numpy.array([rotation_2[i][2] for i in range(3)])

    if numpy.dot(z1_v,z2_v.T) < 0:
        for i in range(3):
            for j in range(1,3):
                rotation_2[i][j]=-rotation_2[i][j]

    y1_vector=[rotation_1[i][1] for i in range(3)]
    y2_vector=[rotation_2[i][1] for i in range(3)]

    gamma=math.acos(numpy.dot(y1_vector,y2_vector)) 
    #  //calculate buckleopening angle.
    b0_vector=numpy.cross(y2_vector,y1_vector)
    # //buckle-opening axis
    
    b0_vector=b0_vector/math.sqrt(numpy.dot(b0_vector,b0_vector.T)) 
    #   //   normalize the b0 

    rotate_matrix_1=DNA_matrix.Rotate_matrix(-gamma/2,b0_vector)
    rotate_matrix_2=DNA_matrix.Rotate_matrix(gamma/2,b0_vector)
    #//create rotate_matrix_1 
    #//create rotate_matrix_2 

    trans_orient_1=numpy.matrix(rotate_matrix_1)*numpy.matrix(rotation_1) 
    trans_orient_1=DNA_matrix.Norm_matrix_in_row(trans_orient_1)
    
    trans_orient_2=numpy.matrix(rotate_matrix_2)*numpy.matrix(rotation_2) 
    trans_orient_2=DNA_matrix.Norm_matrix_in_row(trans_orient_2)
#//get the orientation matirx of transformed bases.

    middle_matirx=(trans_orient_1+trans_orient_2)/2 

    middle_matirx=DNA_matrix.Norm_matrix_in_row(middle_matirx)
   #normalization the row.

    middle_origin=(numpy.array(origin_1)+numpy.array(origin_2))/2 

    return numpy.array(middle_matirx),middle_origin



def Calc_dihedral(atom1, atom2, atom3, atom4):
    b1=numpy.array([atom2.atom_coor_x-atom1.atom_coor_x,\
            atom2.atom_coor_y-atom1.atom_coor_y,\
            atom2.atom_coor_z-atom1.atom_coor_z])

    b2=numpy.array([atom3.atom_coor_x-atom2.atom_coor_x,\
            atom3.atom_coor_y-atom2.atom_coor_y,\
            atom3.atom_coor_z-atom2.atom_coor_z])

    b3=numpy.array([atom4.atom_coor_x-atom3.atom_coor_x,\
            atom4.atom_coor_y-atom3.atom_coor_y,\
            atom4.atom_coor_z-atom3.atom_coor_z])

    phi=math.atan2( numpy.dot( math.sqrt(numpy.dot(b2,b2))*b1 , numpy.cross(b2,b3) ),\
            numpy.dot( numpy.cross(b1,b2) , numpy.cross(b2,b3) )) *180 /math.pi
    if phi < -120:
        phi_c="t"
    elif phi < 0:
        phi_c="g-"
    elif phi < 120:
        phi_c="g+"
    else:
        phi_c="t"


    return phi,phi_c



def Get_Dihedral(Atom_list, base_serial):
    '''
    calculate the dihedral from the top file.

    '''

    for i in Atom_list:
  #      print Atom_list[i].atom_name
        if Atom_list[i].atom_name   in ["O3'","O3*"] and Atom_list[i].residue_serial==base_serial-1: 
            index_O3_0=i

        elif Atom_list[i].atom_name in ["P","H5T"]   and Atom_list[i].residue_serial==base_serial:
            index_P=i

        elif Atom_list[i].atom_name in ["O5'","O5*"] and Atom_list[i].residue_serial==base_serial:
            index_O5=i
        elif Atom_list[i].atom_name in ["C5'","C5*"] and Atom_list[i].residue_serial==base_serial:
            index_C5=i
        elif Atom_list[i].atom_name in ["C4'","C4*"] and Atom_list[i].residue_serial==base_serial:
            index_C4=i
        elif Atom_list[i].atom_name in ["C3'","C3*"] and Atom_list[i].residue_serial==base_serial:
            index_C3=i
        elif Atom_list[i].atom_name in ["O3'","O3*"] and Atom_list[i].residue_serial==base_serial:
            index_O3=i

        elif Atom_list[i].atom_name in ["P","H3T"]   and Atom_list[i].residue_serial==base_serial+1:
            index_P_2=i
        elif Atom_list[i].atom_name in ["O5'","O5*"] and Atom_list[i].residue_serial==base_serial+1:
            index_O5_2=i

        elif Atom_list[i].atom_name in ["O4'","O4*"] and Atom_list[i].residue_serial==base_serial:
            index_O4=i

        elif Atom_list[i].atom_name in ["C1'","C1*"] and Atom_list[i].residue_serial==base_serial:
            index_C1=i

        elif Atom_list[i].atom_name in ["N9"] and Atom_list[i].residue_serial==base_serial and \
                ("A" in Atom_list[i].residue_name or "G" in Atom_list[i].residue_name):
            index_N_base=i
        elif Atom_list[i].atom_name=="N1" and Atom_list[i].residue_serial==base_serial and \
                ("C" in Atom_list[i].residue_name or "T" in Atom_list[i].residue_name):
            index_N_base=i
        elif Atom_list[i].atom_name=="C4" and Atom_list[i].residue_serial==base_serial and \
                ("A" in Atom_list[i].residue_name or "G" in Atom_list[i].residue_name):
            index_C_base=i
        elif Atom_list[i].atom_name=="C2" and Atom_list[i].residue_serial==base_serial and \
                ("C" in Atom_list[i].residue_name or "T" in Atom_list[i].residue_name):
            index_C_base=i

    dihedral=dict()

    try:
        alpha,alpha_c =Calc_dihedral(Atom_list[index_O3_0],Atom_list[index_P ],Atom_list[index_O5 ],Atom_list[index_C5  ]) 
    except:
        alpha = "-"
        alpha_c="-"
    try:
        beta,beta_c  =Calc_dihedral(Atom_list[index_P   ],Atom_list[index_O5],Atom_list[index_C5 ],Atom_list[index_C4  ]) 
    except: 
        beta = "-"
        beta_c = "-"
    gamma,gamma_c =Calc_dihedral(Atom_list[index_O5  ],Atom_list[index_C5],Atom_list[index_C4 ],Atom_list[index_C3  ]) 
    delta,delta_c =Calc_dihedral(Atom_list[index_C5  ],Atom_list[index_C4],Atom_list[index_C3 ],Atom_list[index_O3  ]) 
    try:
        epslon,epslon_c=Calc_dihedral(Atom_list[index_C4  ],Atom_list[index_C3],Atom_list[index_O3 ],Atom_list[index_P_2 ]) 
    except:
        epslon = "-"
        epslon_c="-"

    try:
        zeta,zeta_c  =Calc_dihedral(Atom_list[index_C3  ],Atom_list[index_O3],Atom_list[index_P_2],Atom_list[index_O5_2]) 
    except:
        zeta = "-"
        zeta_c = "-"

    chi ,chi_c  =Calc_dihedral(Atom_list[index_O4],Atom_list[index_C1],Atom_list[index_N_base],Atom_list[index_C_base]) 
#    print index_O4,index_C1,index_N_base,index_C_base
#    print Atom_list[index_O4].atom_coor_x
    if chi < 0 or chi > 90:
        chi_c="anti"
    elif chi < 90 and chi > 0:
        chi_c="syn"


    dihedral[alpha ] = alpha
    dihedral[beta  ] = beta
    dihedral[gamma ] = gamma
    dihedral[delta ] = delta
    dihedral[zeta  ] = zeta
    dihedral[epslon] = epslon
    dihedral[chi   ] = chi
    dihedral[alpha_c ] = alpha_c
    dihedral[beta_c  ] = beta_c
    dihedral[gamma_c ] = gamma_c
    dihedral[delta_c ] = delta_c
    dihedral[epslon_c] = epslon_c
    dihedral[zeta_c  ] = zeta_c
    dihedral[chi_c   ] = chi_c
    return [alpha,beta,gamma,delta,epslon,zeta,chi,alpha_c,beta_c,gamma_c,delta_c,epslon_c,zeta_c,chi_c]

