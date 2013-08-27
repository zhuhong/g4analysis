#1/usr/bin/env python
#-*- coding:utf-8 -*-
'''
Created at 2010.01.14
@version: 0.0.3
@change:\n
    - 2010.01.14.
        - Create this file.
    - 2010.01.18\n
        - Test the function B{Get_rotate_matrix()}
        - Add the function B{Write_rotate_matrix()}
        - Add the function B{Rotate_2_vector()}
        - Finished the function B{Get_parallel_result()},  but not test it.
    - 2010.01.19\n
        - Modified some bugs. 
            - Calculate the angle between the two vectors.
            - Close the file after writting.
            - modified a bug on generate rotation matirx
            - modified a bug on calculate the origin vector for 2 and 4 base in a group.
    - 2011.02.01\n
        - Modified this module, so both *.pdb and *.gro are allowd for coordinate file.
    - 2011.02.18\n
        - Modified a bug. change all the reside to residue.
    - 2011.03.18.\n
        - Rewrite some functions.
    - 2011.03.20.\n
        - write functions for calculate RMSD.
    - 2011.04.03.\n
        - modified the output, change from print to usage.echo
'''

import Simple_atom
import time as Time
import sys
import os
import usage
import MDAnalysis
from MDAnalysis import coordinates
from math import cos, sin, sqrt
import math
import numpy
from numpy import matrix

import DNA_matrix


def Get_rise_fromTRJ(traj_file, coor_file, base_list_1, base_list_2, output_name,skip=1,dt=1,begin=0,end=-1):
    '''
    Reading the traj file and the coordinate file like *.pdb or *.gro. With the base serial choosed,  get the
    rotation matrix for this base. and write it's to a output file with some syntax.
    base_list_1 format list(list())
    base_list_2 format list(list())
    output_name format list()
    '''
    print "  init......"
#    if coor_file.endswith('.top'):
#        Atom_list=amber_top
    START_TIME=Time.time()
    LIST_NUM=len(base_list_1)

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
        fp.write("#time(ns)   distance(A)\t   angle(degree)    RMSD1(A)\t    RMSD2(A)\n")
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


    try:
        u=MDAnalysis.Universe(coor_file,traj_file).trajectory
    except:
        if traj_file.endswith("mdcrd"):
            u=coordinates.TRJ.TRJReader(traj_file,len(Atom_list),delta=dt)
        else:
            print "ERROR: The trajectory file %s not support new." %traj_file
            sys.exit()

    if traj_file.endswith("mdcrd") or traj_file.endswith("dcd"):
        pass
    else:
        try:
            dt=u.dt
        except:
            dt=0.0

    for ts in u:
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
    
                if numpy.dot(orient_group_1, orient_group_2)<0:
                    orient_group_2 = orient_group_2*(-1)

                orient_total = DNA_matrix.Rotate_2_vector(orient_group_1, orient_group_2)
                dist_vector = numpy.array(origin_group_2)-numpy.array(origin_group_1)

                dist = abs(numpy.dot(orient_total, dist_vector))
                gamma = numpy.dot(orient_group_1, orient_group_2) 
                gamma = math.acos(gamma)*180/3.1416


                if ts.frame % 10 ==0 and i==0:
                    NOW_TIME=Time.time()
                    if time < 1000:
                        usage.echo("  analysis frame %6d, time %8.1f ps, time used %8.2f s\r" %(ts.frame, time,NOW_TIME-START_TIME))
                    elif time > 1000 and ts.frame %200 == 0:
                        usage.echo("  analysis frame %6d, time %8.2f ns, time used %8.2f s\r" %(ts.frame, time/1000,NOW_TIME-START_TIME))

                fp = open(output_name[i], 'a')
                fp.write( " %7.4f\t  %6.3f\t   %6.3f\t    %6.4f\t   %6.4f\n" %(time/1000, dist, gamma,RMSD1,RMSD2))
                fp.close()

    print "The parallel analysis finished"
    print "The result are in file: %s" %output_name

def Get_rise_fromTOP( coor_file, base_list_1, base_list_2):
    '''
    Get the parallel parameters from coor_file.
    '''
#    if coor_file.endswith('.top'):
#        Atom_list=amber_top

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

#    print base_atom_list_1
#    print base_atom_list_2


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
        c1.append(numpy.array(temp_list))
        r1.append(result)

    for m in range(len(base_name_list_2)):
        temp_list = [ [atom_list[x].atom_coor_x*10, atom_list[x].atom_coor_y*10,atom_list[x].atom_coor_z*10] \
                for x in base_atom_list_2[m] ]
#        print temp_list
        result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_2[m][0])
        c2.append(numpy.array(temp_list))
        r2.append(result)

    orient_group_1,origin_group_1 = DNA_matrix.Get_group_rotmat(r1,len(base_name_list_1))
    orient_group_2,origin_group_2 = DNA_matrix.Get_group_rotmat(r2,len(base_name_list_2))
    RMSD1=DNA_matrix.Get_group_RMSD(base_name_list_1,c1,origin_group_1,orient_group_1)
    RMSD2=DNA_matrix.Get_group_RMSD(base_name_list_2,c2,origin_group_2,orient_group_2)

    if numpy.dot(orient_group_1, orient_group_2)<0:
        orient_group_2 = orient_group_2*(-1)

    orient_total = DNA_matrix.Rotate_2_vector(orient_group_1, orient_group_2)
    dist_vector = numpy.array(origin_group_2)-numpy.array(origin_group_1)

    dist = abs(numpy.dot(orient_total, dist_vector))
    gamma = numpy.dot(orient_group_1, orient_group_2) 
    gamma = math.acos(gamma)*180/3.1416

    print  "rise:  %6.3f\tangle:   %6.3f\tRMSD_1:    %6.4f\tRMSD_2:   %6.4f\n" %(dist, gamma,RMSD1,RMSD2)


def Get_RMSD_fromTRJ(traj_file, coor_file, base_list, output_name,skip=1, dt=1,begin=0,end=-1):
    '''
    Reading the traj file and the coordinate file like *.pdb or *.gro. With the base serial choosed,  get the
    rotation matrix for this base. and write it's to a output file with some syntax.
    '''
    print "  init......"
    START_TIME=Time.time()
    LIST_NUM=len(base_list)

    
    Atom_list=Simple_atom.Get_Simple_atom_list(coor_file)
    residue_list=Simple_atom.Get_Residue_list(Atom_list)
    
    base_name_list=list()
    base_atom_list=list()

    for i in range(LIST_NUM):
        if os.path.isfile(output_name[i]):
            print "backup %s to %s" %(output_name[i],"#"+output_name[i]+"#")
            try:
                os.rename(output_name[i],"#"+output_name[i]+"#")
            except OSError,e: 
                print e
                print "the file %s will be overwrited!" %output_name[i]

        fp = open(output_name[i], 'w')
        fp.write("#Group list: ")
        for j in base_list[i]:
            fp.write("%d\t " %j)
        fp.write("\n")
        fp.write("#skip:%d\n" %skip)
        fp.write("#time(ns)   RMSD(A)\n")

#        base_name_list.append( [residue_list[j-1] for j in base_list[i]])
        temp_list=list()
        for m in base_list[i]:
            for n in residue_list:
                if n[1]==m:
                    temp_list.append(n)
                    break
                else:
                    pass
        base_name_list.append(temp_list)
#        print base_name_list

        base_atom_list.append( [DNA_matrix.Get_baseID_list(Atom_list,j) for j in base_list[i]])

    #u=MDAnalysis.Universe(coor_file,traj_file)
    try:
        u=MDAnalysis.Universe(coor_file,traj_file).trajectory
    except:
        if traj_file.endswith("mdcrd"):
            u=coordinates.TRJ.TRJReader(traj_file,len(Atom_list),delta=dt)
        else:
            print "ERROR: The trajectory file %s not support new." %traj_file
            sys.exit()

    if traj_file.endswith("mdcrd") or traj_file.endswith("dcd"):
        pass
    else:
        dt=u.dt


    for ts in u:
        time=(ts.frame-1)*dt
   #     print begin,end
        if time < float(begin):
   #         print time
            continue
        if end > 0 and time > float(end):
            break

        if ts.frame % skip == 0 :
            for i in range(LIST_NUM):
                r1=[]
                '''the group 1 rotate list'''
                c1=[]
                '''the group 1 coordinate list'''
                for m in range(len(base_name_list[i])):
                    temp_list = [ [ts._x[x-1], ts._y[x-1], ts._z[x-1]] for x in base_atom_list[i][m] ]
                    result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list[i][m][0])
                    c1.append(numpy.array(temp_list))
                    r1.append(result)

                orient_group,origin_group = DNA_matrix.Get_group_rotmat(r1,len(base_name_list[i]))
                RMSD=DNA_matrix.Get_group_RMSD(base_name_list[i],c1,origin_group,orient_group)

                time=(ts.frame-1)*dt

                if ts.frame % 100 ==0:
                    NOW_TIME=Time.time()
                    usage.echo("  analysis frame %6d, time %8.1f ps, time used %8.2f s\r" %(ts.frame, time,NOW_TIME-START_TIME))
                fp=open(output_name[i], 'a')
                fp.write( " %7.4f\t  %6.4f\n" %(time/1000,RMSD))
                fp.close()

    print "Finished calculating the RMSD."
    print "The result are in file: %s" %output_name


def Get_RMSD_fromTOP(coor_file, base_list):
    '''
    Calculate the RMSD from coor_file.
    '''
    LIST_NUM=len(base_list)

    atom_list=Simple_atom.Get_Simple_atom_list(traj_file)
    residue_list=Simple_atom.Get_Segment_list(atom_list)

    base_name_list=list()
    base_atom_list=list()
    for i in range(LIST_NUM):
        base_name_list.append( [residue_list[j-1] for j in base_list[i]])
        base_atom_list.append( [DNA_matrix.Get_baseID_list(atom_list,j) for j in base_list[i]])

    for i in range(LIST_NUM):
        r1=[]
        '''the group 1 rotate list'''
        c1=[]
        '''the group 1 coordinate list'''
        for m in range(len(base_name_list[i])):
            temp_list = [ [atom_list[x-1].atom_coor_x, atom_list[x-1].atom_coor_y, atom_list[x-1].atom_coor_z] \
                    for x in base_atom_list[i][m] ]
            result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list[i][m][0])
            c1.append(numpy.array(temp_list))
            r1.append(result)

        orient_group,origin_group = DNA_matrix.Get_group_rotmat(r1,len(base_name_list[i]))
        RMSD=DNA_matrix.Get_group_RMSD(base_name_list[i],c1,origin_group,orient_group)

        print "RMSD=%f" %RMSD
    return 


