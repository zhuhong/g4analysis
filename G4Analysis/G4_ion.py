#1/usr/bin/env python
#-*- coding:utf-8 -*-
'''
'''

import Simple_atom
import time as Time
import sys
import os
import usage
import MDAnalysis
from math import cos, sin, sqrt
import math
import numpy
from numpy import matrix

import DNA_matrix

def Get_dist(traj_file, coor_file, G4_bases, ion_residue, output_name,skip=1, dt=1,begin=0,end=-1):
    print "  init......"
#    if coor_file.endswith('.top'):
#        Atom_list=amber_top
    START_TIME=Time.time()

    Atom_list=Simple_atom.Get_Simple_atom_list(coor_file)
    residue_list=Simple_atom.Get_Residue_list(Atom_list)

    base_name_list_1=list()
    base_atom_list_1=list()


    if os.path.isfile(output_name):
        print "backup %s to %s" %(output_name,"#"+output_name+"#")
        try:
            os.rename(output_name,"#"+output_name+"#")
        except OSError,e: 
            print e
            print "the file %s will be overwrited!" %output_name

    fp = open(output_name, 'w')
    fp.write("#Group 1: ")
    for j in G4_bases:
        fp.write("%d\t " %j)
    fp.write("\n")
    fp.write("#Group 2: ")
    fp.write("%d\t " %ion_residue)
    fp.write("\n")
    fp.write("#skip:%d\n" %skip)
    fp.write("#time(ns)   distance(A)\n")
    fp.close()

    for m in G4_bases:
        for n in residue_list:
            if n[1]==m:
                base_name_list_1.append(n)
                break
            else:
                pass

    base_atom_list_1=[DNA_matrix.Get_baseID_list(Atom_list,j) for j in G4_bases]
    ion_serial=Simple_atom.Get_Atom_in_residue(Atom_list,ion_residue)[0].atom_serial
    print ion_serial

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
            r1=[]
            '''the group 1 rotate list'''
            c1=[]
            '''the group 1 coordinate list'''
            for m in range(len(base_name_list_1)):
                temp_list = [ [ts._x[x-1], ts._y[x-1], ts._z[x-1]] for x in base_atom_list_1[m] ]
                result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_1[m][0])
                r1.append(result)

            c2=  [ts._x[ion_serial-1], ts._y[ion_serial-1], ts._z[ion_serial-1]] 

            orient_group_1,origin_group_1 = DNA_matrix.Get_group_rotmat(r1,len(base_name_list_1))
    

            dist_vector = numpy.array(c2)-numpy.array(origin_group_1)

            dist = numpy.dot(orient_group_1, dist_vector)

            if ts.frame % 10 ==0:
                NOW_TIME=Time.time()
                usage.echo("  analysis frame %6d, time %8.1f ps, time used %8.2f s\r" %(ts.frame, time,NOW_TIME-START_TIME))

            fp = open(output_name, 'a')
            fp.write( " %7.4f\t  %6.3f\n" %(time/1000, dist))
            fp.close()

    print "The parallel analysis finished"
    print "The result are in file: %s" %output_name
