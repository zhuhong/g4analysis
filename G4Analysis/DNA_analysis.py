#1/usr/bin/env python
#-*- coding:utf-8 -*-
'''
Created at 2012.03.06
'''

import Simple_atom
import time as Time
import sys
import os
import usage
import MDAnalysis
import DNA_param
from math import cos, sin, sqrt
import math
import numpy
from numpy import matrix
from numpy import dot

import DNA_matrix

def Get_parm_fromTRJ(traj_file, coor_file, base_list_1, base_list_2, output_name,CALCU="step",skip=1, dt=1,begin=0,end=-1):
    '''
    Note: this function not finish yet.
    base_list_1=[[base or base pair],[base or base pair],...]
    base_list_2=[[base or base pair],[base or base pair],...]
    '''
    print "  init......"
#    if coor_file.endswith('.top'):
#        Atom_list=amber_top
    START_TIME=Time.time()
    if len(base_list_1)==len(base_list_2):
        LIST_NUM=len(base_list_1)
    else:
        print "ERROR: The length of the base list not match."
        return -1
    if len(output_name)!=LIST_NUM:
        print "ERROR: The number of the output file not match the size of the base list."

    Atom_list=Simple_atom.Get_Simple_atom_list(coor_file)
    residue_list=Simple_atom.Get_Residue_list(Atom_list)

    base_name_list_1=list()
    base_name_list_2=list()
    base_atom_list_1=list()
    base_atom_list_2=list()

    for i in range(LIST_NUM):

        if os.path.isfile(output_name[i]):
            print "backup %s to %s" %(output_name[i],"#"+output_name[i]+"#")
            try:
                os.rename(output_name[i],"#"+output_name[i]+"#")
            except OSError,e: 
                print e
                print "the file %s will be overwrited!" %output_name[i]

        fp = open(output_name[i], 'w')
        fp.write("#Group 1: ")
        for j in base_list_1[i]:
            fp.write("%d\t " %j)
        fp.write("\n")
        fp.write("#Group 2: ")
        for j in base_list_2[i]:
            fp.write("%d\t " %j)
        fp.write("\n")
        fp.write("#skip:%d\n" %skip)
        if CALCU=="step":
            fp.write("#  shift       slide        rise        tilt        roll       twist\n")
        else:
            fp.write("#  shear     stretch     stagger      buckle   propeller     opening\n")
        fp.close()

#        base_name_list_1.append( [residue_list[j-1] for j in base_list_1[i]])
#        base_name_list_2.append( [residue_list[j-1] for j in base_list_2[i]])

        temp_list=list()
        for m in base_list_1[i]:
            for n in residue_list:
                if n[1]==m:
                    temp_list.append(n)
                    break
                else:
                    pass
        base_name_list_1.append(temp_list)

        temp_list=list()
        for m in base_list_2[i]:
            for n in residue_list:
                if n[1]==m:
                    temp_list.append(n)
                    break
                else:
                    pass
        base_name_list_2.append(temp_list)


        base_atom_list_1.append([DNA_matrix.Get_baseID_list(Atom_list,j) for j in base_list_1[i]])
        base_atom_list_2.append([DNA_matrix.Get_baseID_list(Atom_list,j) for j in base_list_2[i]])


    u=MDAnalysis.Universe(coor_file,traj_file)

    if traj_file.endswith("mdcrd") or traj_file.endswith("dcd"):
        pass
    else:
        try:
            dt=u.trajectory.dt
        except:
            dt=0.0

    for ts in u.trajectory:
        time=float((ts.frame-1)*dt)
        if dt > 0.0:
            if time < float(begin):
                continue
            if time > float(end) and end !=-1:
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
                    r1.append(result[0])
                    c1.append(result[1])

                for m in range(len(base_name_list_2[i])):
                    temp_list = [ [ts._x[x-1], ts._y[x-1], ts._z[x-1]] for x in base_atom_list_2[i][m] ]
                    result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_2[i][m][0])
                    r2.append(result[0])
                    c2.append(result[1])

                fp = open(output_name[i], 'a')

                if CALCU=="pair":
                    a=DNA_param.base_pair_parameters(r1[0],r2[0],c1[0],c2[0])
                    fp.write("%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f\n" %(a[0],a[1],a[2],a[3],a[4],a[5]))
                else:
                    middle_r1,middle_c1=DNA_param.middle_frame(r1[0],r1[1],c1[0],c1[1])
                    middle_r2,middle_c2=DNA_param.middle_frame(r2[0],r2[1],c2[0],c2[1])
                    a=DNA_param.base_step_parameters(middle_r1,middle_r2,middle_c1,middle_c2)
                    fp.write("%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f\n" %(a[0],a[1],a[2],a[3],a[4],a[5]))

                if ts.frame % 10 ==0 and i==0:
                    NOW_TIME=Time.time()
                    usage.echo("  analysis frame %6d, time %8.1f ps, time used %8.2f s\r" %(ts.frame, time,NOW_TIME-START_TIME))

                fp.close()

    print "The DNA helical analysis finished"
    print "The result are in file: %s" %output_name

def Get_para_fromTOP( coor_file, base_list_1, base_list_2,CALCU="step",PRINT=True):
    '''
    get pair parameters from coor_file
    '''

    Atom_list=Simple_atom.Get_Simple_atom_list(coor_file)
    residue_list=Simple_atom.Get_Residue_list(Atom_list)
    atom_list=Simple_atom.Get_atom_list(coor_file)

    base_name_list_1=list()
    base_name_list_2=list()
    base_atom_list_1=list()
    base_atom_list_2=list()

    for m in base_list_1:
        for n in residue_list:
            if n[1]==m:
                base_name_list_1.append(n)
                break
            else:
                pass

    for m in base_list_2:
        for n in residue_list:
            if n[1]==m:
                base_name_list_2.append(n)
                break
            else:
                pass

    base_atom_list_1=[DNA_matrix.Get_baseID_list(Atom_list,j) for j in base_list_1]
    base_atom_list_2=[DNA_matrix.Get_baseID_list(Atom_list,j) for j in base_list_2]

    r1=[]
    '''the group 1 rotate list'''
    r2=[]
    '''the group 2 rotate list'''
    c1=[]
    '''the group 1 coordinate list'''
    c2=[]
    '''the group 2 coordinate list'''
    for m in range(len(base_name_list_1)):
        temp_list = [ [atom_list[x].atom_coor_x*10, atom_list[x].atom_coor_y*10,atom_list[x].atom_coor_z*10] \
                for x in base_atom_list_1[m] ]
#        print temp_list
        result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_1[m][0])
        r1.append(result[0])
        c1.append(result[1])

    for m in range(len(base_name_list_2)):
        temp_list = [ [atom_list[x].atom_coor_x*10, atom_list[x].atom_coor_y*10,atom_list[x].atom_coor_z*10] \
                for x in base_atom_list_2[m] ]
#        print temp_list
        result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_2[m][0])
        r2.append(result[0])
        c2.append(result[1])

    if CALCU=="pair":
        a=DNA_param.base_pair_parameters(r1[0],r2[0],c1[0],c2[0])
        if PRINT:
            print "   shear     stretch     stagger      buckle   propeller     opening"
            print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])
        return a
    else:
        middle_r1,middle_c1=DNA_param.middle_frame(r1[0],r1[1],c1[0],c1[1])
        middle_r2,middle_c2=DNA_param.middle_frame(r2[0],r2[1],c2[0],c2[1])
        a=DNA_param.base_step_parameters(middle_r1,middle_r2,middle_c1,middle_c2)
        if PRINT:
            print "   shift       slide        rise        tilt        roll       twist"
            print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])
        return a


def Get_Dihedral_fromTOP( coor_file, base_list_1,PRINT=False ):
    '''
    calculate the dihedral from the top file.

    '''
    Atom_list=Simple_atom.Get_atom_list(coor_file)

    if PRINT == True:
        print "%8s%10s%10s%10s%10s%10s%10s%10s%8s%8s%8s%8s%8s%8s%8s"\
                %("baseID","alpha","beta","gamma","delta","epslon","zeta","chi",\
                "alpha","beta","gamma","delta","epslon","zeta","chi")

    for m in base_list_1:
        resu=DNA_param.Get_Dihedral(Atom_list,m)   

        if PRINT == True:
            print "%8d" %m,
            for i in range(7):
                if resu[i]=="-":
                    print "%9s" %("-"*4), 
                else:
                    print "%9.2f" %(resu[i]),
            for i in range(7):
                if resu[i]=="-":
                    print "%7s" %("-"*4), 
                else:
                    print "%7s" %(resu[i+7]),

            print "" 

#a        return [alpha,beta,gamma,delta,epslon,zeta,chi]


def Get_Dihedral_fromTRJ(traj_file, coor_file, base_list, output_name,skip=1, dt=1,begin=0,end=-1):
    '''
    note: if the atom index in coor file not start with 1, there will be an error.
    '''
    START_TIME=Time.time()

    Atom_list=Simple_atom.Get_atom_list(coor_file)


    for file_name in output_name:

        if os.path.isfile(file_name):
            print "backup %s to %s" %(file_name,"#"+file_name+"#")
            try:
                os.rename(file_name,"#"+file_name+"#")
            except OSError,e: 
                print e
                print "the file %s will be overwrited!" %file_name

        fp = open(file_name, 'w')
        fp.write("#Residue Index: ")
        print base_list[output_name.index(file_name)]
        fp.write("%d\n" %base_list[output_name.index(file_name)])
        fp.write("#skip:%d\n" %skip)

        fp.write("%12s%10s%10s%10s%10s%10s%10s%10s%8s%8s%8s%8s%8s%8s%8s\n"\
                %("Time(ns)","alpha","beta","gamma","delta","epslon","zeta","chi",\
                "alpha","beta","gamma","delta","epslon","zeta","chi"))
        fp.close()

    u=MDAnalysis.Universe(coor_file,traj_file)

    if traj_file.endswith("mdcrd") or traj_file.endswith("dcd"):
        pass
    else:
        try:
            dt=u.trajectory.dt
        except:
            dt=0.0

    for ts in u.trajectory:
        time=float((ts.frame-1)*dt)
        if dt > 0.0:
            if time < float(begin):
                continue
            if time > float(end) and end !=-1:
                break

        if ts.frame % skip == 0 :
            for i in Atom_list:
                Atom_list[i].atom_coor_x=ts._x[i-1]
                Atom_list[i].atom_coor_y=ts._y[i-1]
                Atom_list[i].atom_coor_z=ts._z[i-1]


            for i in range(len(base_list)):
                fp = open(output_name[i], 'a')
                resu=DNA_param.Get_Dihedral(Atom_list,base_list[i])
                fp.write("%12.4f" %(time/1000))
                for j in range(7):
                    if resu[j]=="-":
                        fp.write("%10s" %("-"*4)) 
                    else:
                        fp.write("%10.2f" %(resu[j]))

                for j in range(7):
                    if resu[j]=="-":
                        fp.write("%8s" %("-"*4)) 
                    else:
                        fp.write("%8s" %(resu[j+7]))

                fp.write("\n") 


                if ts.frame % 10 ==0 and i==0:
                    NOW_TIME=Time.time()
                    if time < 1000:
                        usage.echo("  analysis frame %6d, time %8.1f ps, time used %8.2f s\r" %(ts.frame, time,NOW_TIME-START_TIME))
                    elif time > 1000 and ts.frame %100 ==0 :
                        usage.echo("  analysis frame %6d, time %8.2f ns, time used %8.2f s\r" %(ts.frame, time/1000,NOW_TIME-START_TIME))

                fp.close()

    print "\nThe DNA helical analysis finished"
    print "The result are in file: %s" %output_name


