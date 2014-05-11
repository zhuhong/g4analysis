import Simple_atom
import numpy
import os
import sys
import MDAnalysis
import time as Time

'''
Should note that unit is A and A^2.
'''

def Get_Area_fromTOP(coor_file, base_list):
    '''
    Calculate the RMSD from coor_file.
    '''
    LIST_NUM=len(base_list)

    atom_list=Simple_atom.Get_atom_list(coor_file)
    residue_list=Simple_atom.Get_Residue_list(atom_list)

    base_name_list=list()
    base_atom_list=list()
    for i in range(LIST_NUM):
        base_name_list.append( [residue_list[j-1] for j in base_list[i]])
        base_atom_list.append( [_Get_O6(atom_list,j) for j in base_list[i]])
    # print base_atom_list

    for i in range(LIST_NUM):
        r1=[]
        '''the group 1 rotate list'''
        c1=[]
        '''the group 1 coordinate list'''
        temp_list = [ [atom_list[x-1].atom_coor_x*10, atom_list[x-1].atom_coor_y*10, atom_list[x-1].atom_coor_z*10] \
                for x in base_atom_list[i] ]
        # print temp_list
        # result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list[i][m][0])
        vect = list()
        for x in range(4):
            for y in range(x):
                if x > y:
                    vect.append(numpy.array([temp_list[x][z] - temp_list[y][z] for z in range(3)]))
        # print vect
        vect2 = [numpy.linalg.norm(vv) for vv in vect]
        # print vect2
        vect2_sort = sorted(vect2,reverse = True) 
        # print vect2_sort

        ind_1 = vect2.index(vect2_sort[0])
        ind_2 = vect2.index(vect2_sort[1])
        # print ind_1,ind_2

        # print numpy.cross(vect[ind_1],vect[ind_2])
        area = 0.5*numpy.linalg.norm(numpy.cross(vect[ind_1],vect[ind_2]))
        # print area

        print "area=%f" %area
    return 


def _Get_O6(atom_list,residue_id):
	for atom in atom_list:
		if atom.atom_name == "O6" and atom.residue_id == residue_id:
			return atom.atom_id
		else:
			continue


def Get_Area_fromTRJ(traj_file, coor_file, base_list, output_name,skip=1,dt=1,begin=0,end=-1):
    '''
    Reading the traj file and the coordinate file like *.pdb or *.gro. With the base serial choosed,  get the
    rotation matrix for this base. and write it's to a output file with some syntax.
    base_list format list(list())
    output_name format list()
    '''
    print "  init......"
#    if coor_file.endswith('.top'):
#        Atom_list=amber_top
    START_TIME=Time.time()
    LIST_NUM=len(base_list)

    atom_list=Simple_atom.Get_atom_list(coor_file)
    residue_list=Simple_atom.Get_Residue_list(atom_list)

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
        fp.write("#Group 1: ")
        for j in base_list[i]:
            fp.write("%d\t " %j)
        fp.write("\n")
        fp.write("#skip:%d\n" %skip)
        fp.write("#time(ns)   Area(A^2)\n")
        fp.close()

#        base_name_list_1.append( [residue_list[j-1] for j in base_list_1[i]])
#        base_name_list_2.append( [residue_list[j-1] for j in base_list_2[i]])

        base_name_list.append([residue_list[j-1]    for j in base_list[i]])
        base_atom_list.append([_Get_O6(atom_list,j) for j in base_list[i]])


    try:
        u=MDAnalysis.Universe(coor_file,traj_file).trajectory
    except:
        if traj_file.endswith("mdcrd"):
            u=coordinates.TRJ.TRJReader(traj_file,len(atom_list),delta=dt)
        else:
            print "ERROR: The trajectory file %s not support now." %traj_file
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
                temp_list = [ [ts._x[x-1], ts._y[x-1], ts._z[x-1]] for x in base_atom_list[i] ]

                vect = list()
                for x in range(4):
                    for y in range(x):
                        if x > y:
                            vect.append(numpy.array([temp_list[x][z] - temp_list[y][z] for z in range(3)]))
                # print vect
                vect2 = [numpy.linalg.norm(vv) for vv in vect]
                # print vect2
                vect2_sort = sorted(vect2,reverse = True) 
                # print vect2_sort
            
                ind_1 = vect2.index(vect2_sort[0])
                ind_2 = vect2.index(vect2_sort[1])
                # print ind_1,ind_2
            
                # print numpy.cross(vect[ind_1],vect[ind_2])
                area = 0.5*numpy.linalg.norm(numpy.cross(vect[ind_1],vect[ind_2]))
                # print area
            
                # print "area=%f" %area

                fp = open(output_name[i], 'a')
                fp.write( " %7.4f\t  %6.3f\n" %(time/1000, area))
                fp.close()
