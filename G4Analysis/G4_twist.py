#*- coding : utf-8 -*-

'''
Created on 2011-01-18\n
The twist angle for G-DNA is defined from a JCTC article. which DOI is 10.1021/ct100253m.
The twist angle is defined using the angle between the line of C1' atoms in a G-quartet layer. 
@version: 0.1.0
@author: zhuh
@change:
    - 2011-01-18\n
        - Create this file.
    - 2011-01-19\n
        - finish the function B{Get_twist_in_GDNA()} and test it
    - 2011-01-25\n
        - modified the function B{Get_twist_in_GDNA()}, \
                so both gro and pdb can be used for coor_file.
        - add the version to B{0.1.0}
'''

import MDAnalysis
import numpy
import math
import Simple_atom
import usage
import DNA_matrix
import os

def Get_twist_in_GDNA(traj_file,coor_file,base_list_1,base_list_2,output_file):
    '''
    Input the layer 1 (G11,G12,G13,G14) and layer 2 (G21,G22,G23,G24),Calculate the angle 
    between G1i-G1(i+1) and G2i-G2(i+1). write the result to output_file.\n
    B{traj_file:} the GMX trajectory file, in trr or xtc format.\n
    B{coor_file:} the GMX coordinate file, in pdb or gro format.\n
    B{base_list_1:} the frist group contain four guanine bases.\n
    B{base_list_2:} the second group contain four guanine bases.\n
    B{output_file:} the output file. 
    '''
    C1_list_1=[]
    C1_list_2=[]
    print " init......"
    fp=open(output_file,"w")
    fp.write("#Group 1: ")
    for i in base_list_1:
        fp.write("%d\t " %i)
    fp.write("\n")
    fp.write("#Group 2: ")
    for i in base_list_2:
        fp.write("%d\t " %i)
    fp.write("\n")
    fp.write("#time\t angle_1 \t angle _2 \t angle _3 \t angle_4\n")
    Atom_list=Simple_atom.Get_Simple_atom_list(coor_file)

    for base in base_list_1:
        atom_list=Simple_atom.Get_Atom_in_residue(Atom_list,base)
        for atom in atom_list:
            if atom.atom_name =="C1'":
                C1_list_1.append(atom.atom_serial)

    for base in base_list_2:
        atom_list=Simple_atom.Get_Atom_in_residue(Atom_list,base)
        for atom in atom_list:
            if atom.atom_name =="C1'":
                C1_list_2.append(atom.atom_serial)

#    print C1_list_1
#    print C1_list_2

    u=MDAnalysis.Universe(coor_file,traj_file)    
    for ts in u.trajectory:
        angle=[]
        for i in range(4):
            vector1=[ts._x[C1_list_1[(i+1)%4]-1]-ts._x[C1_list_1[i]-1],\
                    ts._y[C1_list_1[(i+1)%4]-1]-ts._y[C1_list_1[i]-1],\
                    ts._z[C1_list_1[(i+1)%4]-1]-ts._z[C1_list_1[i]-1]]
            vector2=[ts._x[C1_list_2[(i+1)%4]-1]-ts._x[C1_list_2[i]-1],\
                    ts._y[C1_list_2[(i+1)%4]-1]-ts._y[C1_list_2[i]-1],\
                    ts._z[C1_list_2[(i+1)%4]-1]-ts._z[C1_list_2[i]-1]]
            vector1=numpy.array(vector1)
            vector2=numpy.array(vector2)
#            print vector1,vector2
            gamma=numpy.dot(vector1,vector2)/(math.sqrt(numpy.dot(vector1,vector1)*numpy.dot(vector2,vector2)))
            angle.append(math.acos(gamma)/3.1416*180)

        fp.write("%6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f\n" \
                %(ts.time/1000,angle[0],angle[1],angle[2],angle[3]))             
        #if ts.frame % 100 ==0:
        #    print " analysis frame %6d......" %ts.frame
        usage.echo(" analysis frame %6d......\r" %ts.frame)
    fp.close()
    print "The result are in the file: ",output_file


def Get_twist_in_GDNA2(traj_file,coor_file,base_list_1,base_list_2,output_name,skip=1,dt=1,begin=0,end=-1):
    '''
    Input the layer 1 (G11,G12,G13,G14) and layer 2 (G21,G22,G23,G24),Calculate the angle 
    between G1i-G1(i+1) and G2i-G2(i+1). write the result to output_file.\n
    B{traj_file:} the GMX trajectory file, in trr or xtc format.\n
    B{coor_file:} the GMX coordinate file, in pdb or gro format.\n
    B{base_list_1:} the frist group contain four guanine bases.\n
    B{base_list_2:} the second group contain four guanine bases.\n
    B{output_file:} the output file. 
    '''

    LIST_NUM=len(base_list_1)

    Atom_list=Simple_atom.Get_Simple_atom_list(coor_file)
    residue_list=Simple_atom.Get_Residue_list(Atom_list)
#    print residue_list

    base_name_list_1=list()
    base_name_list_2=list()
    base_atom_list_1=list()
    base_atom_list_2=list()

    C1_list_1=list()
    C1_list_2=list()

    print " init......"

    for i in range(LIST_NUM):

        if os.path.isfile(output_name[i]):
            print "backup %s to %s" %(output_name[i],"#"+output_name[i]+"#")
            try:
                os.rename(output_name[i],"#"+output_name[i]+"#")
            except OSError,e: 
                print e
                print "the file %s will be overwrited!" %output_name[i]

        fp=open(output_name[i],"w")
        fp.write("#Group 1: ")
        for li in base_list_1[i]:
            fp.write("%d\t " %li)
        fp.write("\n")
        fp.write("#Group 2: ")
        for li in base_list_2[i]:
            fp.write("%d\t " %li)
        fp.write("\n")
        fp.write("#time\t angle_1 \t angle _2 \t angle _3 \t angle_4\n")
        fp.close()

        base_name_list_1.append( [residue_list[j-1] for j in base_list_1[i]])
#        print base_name_list_1
        base_name_list_2.append( [residue_list[j-1] for j in base_list_2[i]])
#        print base_name_list_2
        base_atom_list_1.append([DNA_matrix.Get_baseID_list(Atom_list,j) for j in base_list_1[i]])
#        print base_atom_list_1
        base_atom_list_2.append([DNA_matrix.Get_baseID_list(Atom_list,j) for j in base_list_2[i]])
#        print base_atom_list_2


        for base in base_list_1[i]:
            atom_list=Simple_atom.Get_Atom_in_residue(Atom_list,base)
            for atom in atom_list:
                if atom.atom_name =="C1'":
                    C1_list_1.append(atom.atom_serial)
 #       print C1_list_1
     
        for base in base_list_2[i]:
            atom_list=Simple_atom.Get_Atom_in_residue(Atom_list,base)
            for atom in atom_list:
                if atom.atom_name =="C1'":
                    C1_list_2.append(atom.atom_serial)
 #       print C1_list_2
     
    u=MDAnalysis.Universe(coor_file,traj_file)    

    if traj_file.endswith("mdcrd") or traj_file.endswith("dcd"):
        pass
    else:
        dt=u.trajectory.dt

    for ts in u.trajectory:
        time=float((ts.frame-1)*dt)
        if time < begin:
            continue
        if time > end and end !=-1:
            break
        if ts.frame % skip == 0 :
            for i in range(LIST_NUM):
                r1=[]
                '''the group 1 rotate list'''
                r2=[]
                '''the group 2 rotate list'''
                c1=[]
                '''the group 1 coordinate list'''
                c2=[]
                '''the group 2 coordinate list'''
                for m in range(len(base_name_list_1[i])):
                    temp_list = [ [ts._x[x-1], ts._y[x-1], ts._z[x-1]] for x in base_atom_list_1[i][m] ]
                    result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_1[i][m][0])
                    #base_name_list_1[index of the groups][index of the base of group 1][base_name,base_serial]
                    c1.append(numpy.array(temp_list))
                    r1.append(result)

                for m in range(len(base_name_list_2[i])):
                    temp_list = [ [ts._x[x-1], ts._y[x-1], ts._z[x-1]] for x in base_atom_list_2[i][m] ]
                    result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_2[i][m][0])
                    c2.append(numpy.array(temp_list))
                    r2.append(result)

                orient_group_1,origin_group_1 = DNA_matrix.Get_group_rotmat(r1,len(base_name_list_1[i]))
                orient_group_2,origin_group_2 = DNA_matrix.Get_group_rotmat(r2,len(base_name_list_2[i]))
                RMSD1=DNA_matrix.Get_group_RMSD(base_name_list_1[i],c1,origin_group_1,orient_group_1)
                RMSD2=DNA_matrix.Get_group_RMSD(base_name_list_2[i],c2,origin_group_2,orient_group_2)

                angle=[]
                for k in range(4):
                    vector1=[ts._x[C1_list_1[(k+1)%4]-1]-ts._x[C1_list_1[k]-1],\
                            ts._y[C1_list_1[(k+1)%4]-1]-ts._y[C1_list_1[k]-1],\
                            ts._z[C1_list_1[(k+1)%4]-1]-ts._z[C1_list_1[k]-1]]

                    vector2=[ts._x[C1_list_2[(k+1)%4]-1]-ts._x[C1_list_2[k]-1],\
                            ts._y[C1_list_2[(k+1)%4]-1]-ts._y[C1_list_2[k]-1],\
                            ts._z[C1_list_2[(k+1)%4]-1]-ts._z[C1_list_2[k]-1]]

                    vector1=numpy.array(vector1)

                    vector1_1=numpy.cross(numpy.cross(orient_group_1,vector1),orient_group_1)

                    vector2=numpy.array(vector2)

                    vector2_2=numpy.cross(numpy.cross(orient_group_2,vector2),orient_group_2)

 #                   print vector1_1,vector2_2
                    gamma=numpy.dot(vector1_1,vector2_2)/(math.sqrt(numpy.dot(vector1_1,vector1_1)*numpy.dot(vector2_2,vector2_2)))
                    angle.append(math.acos(abs(gamma))/3.1416*180)
         
                fp=open(output_name[i],'a')
                fp.write("%6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f\n" \
                        %(ts.time/1000,angle[0],angle[1],angle[2],angle[3]))             
                if ts.frame % 100 ==0 and i ==0:
                    usage.echo(" analysis frame %6d\r" %ts.frame)
                fp.close()
    print "The result are in the file: ",output_name
    
    
   
            
