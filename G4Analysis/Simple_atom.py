#-*- coding : utf-8 -*-
'''
Create on 2011-01-25
@change: \n
    - 2011-01-25\n
        - created and not finished.
    - 2011-01-27\n
        - It may be finished.
    - 2011-02-17\n
        - Add B{atom_2_PDBformat()} and B{atom_2_GROformat()} to the Simple_atom class.
    - 2011-02-21\n
        - Modified some bugs.
    - 2011-02-25.\n
        - Add function B{Check_list()}.
    - 2011.08.23.\n
        - Move the class Simple_atom to unit_atom and change it's name to unit_atom.
        - Add the version to 0.02.
    - 2011.12.26.\n
        - Add  functions Get_residue and Get_atom in this module.
        - Change Get_Segment_list to Get_Residue_list and modified it.
'''

from Coor import PDB
from Coor import GRO
from Coor import unit_atom
from Coor import amber_top
from Coor import atomlib

import string
import re
import os

def Get_Atom_in_residue(atom_list,resid):
    '''
    Return a atom class list which contain the atoms in a residue.
    '''
    ls = []
    for atom in atom_list:
        if atom.residue_serial == resid:
            ls.append(atom)
    return ls

def Get_Simple_atom_list(filename,crd_file=""):
    '''
    Read in a structure file like pdb,gro. return a Simple_atom class list.
    '''
    atom_list=[]
    if filename.endswith(".pdb"):
        atom_list=PDB.Get_Atom_list(filename)
    elif filename.endswith(".gro"):
        atom_list=GRO.Get_Atom_list(filename)
    elif filename.endswith(".top"):
        if crd_file == "":
            atom_list=amber_top.Read_top(filename)
        else:
            atom_list=amber_top.Read_crd(filename,crd_file)
    else:
        print "file name %s in invalid." %filename
        return []

    simple_list=[]
    for atom in atom_list:
        simple=unit_atom.unit_atom(atom_name=atom.atom_name,\
                atom_serial=atom.atom_serial,\
                residue_name=atom.residue_name,\
                residue_serial=atom.residue_serial,\
                atom_coor_x=atom.atom_coor_x,\
                atom_coor_y=atom.atom_coor_y,\
                atom_coor_z=atom.atom_coor_z)
        simple_list.append(simple)

    return simple_list

def Get_atom_list(filename,crd_file=""):
    '''
    Read in a structure file like pdb,gro. return a Simple_atom class dict.
    rewrite from Get_Simple_atom_list, use dict but not list. for the bug
    which may happen when the atom serial in pdb file not start from 1.
    '''
    atom_list=list()
    if filename.endswith(".pdb"):
        atom_list=PDB.Get_Atom_list(filename)
    elif filename.endswith(".gro"):
        atom_list=GRO.Get_Atom_list(filename)
    elif filename.endswith(".top"):
        if crd_file == "":
            atom_list=amber_top.Read_top(filename)
        else:
            atom_list=amber_top.Read_crd(filename,crd_file)
    else:
        print "file name %s in invalid." %filename
        return []

    simple_list=dict()
    for atom in atom_list:
        simple=unit_atom.unit_atom(atom_name=atom.atom_name,\
                atom_serial=atom.atom_serial,\
                residue_name=atom.residue_name,\
                residue_serial=atom.residue_serial,\
                atom_coor_x=atom.atom_coor_x,\
                atom_coor_y=atom.atom_coor_y,\
                atom_coor_z=atom.atom_coor_z)
        simple_list[simple.atom_serial]=simple

    return simple_list



def Check_list(atom_list):
    '''
    Check the atom list, the list which may be deleted some groups. It's used
    in Modify_Coor.py.
    '''
    i=1
    MAX=len(atom_list)
    temp_residue=0
    while i < MAX:
        if atom_list[i].atom_serial - atom_list[i-1].atom_serial == 1:
            pass
        else:
            if atom_list[i].residue_serial != temp_residue:
                atom_list[i].atom_serial = atom_list[i-1].atom_serial + 1
                temp_residue = atom_list[i].residue_serial
                atom_list[i].residue_serial = atom_list[i-1].residue_serial + 1
            else:
                atom_list[i].atom_serial = atom_list[i-1].atom_serial + 1
                temp_residue = atom_list[i].residue_serial
                atom_list[i].residue_serial = atom_list[i-1].residue_serial 
        i+=1

    return atom_list


def Get_Residue_list(atom_list):
    '''
    return a residue list from atom list
    '''
    seg_list=list()
    seg_list.append([atom_list[0].residue_name,atom_list[0].residue_serial])
    '''
    seg_list=[residue_name,residue_serial]
    '''
    for atom in atom_list:
        if atom.residue_serial > seg_list[-1][1]:
            seg_list.append([atom.residue_name,atom.residue_serial])
    return seg_list

def Get_residue(coor_file,show=True):
    '''
    Get_residue is rewrite from Get_list. 
    '''
    atom_list=Get_Simple_atom_list(coor_file)
    reside_list=Get_Residue_list(atom_list)
    chain=[]
    for reside in reside_list:
        if reside[0] in atomlib.RESIDUE_NAME_LIST:
            chain.append(reside)  
        else:
            pass
    if show== True:
        for ca in chain:
            print "%4d\t (%8s )" % (ca[1],ca[0])

    while True:
        list1=raw_input("Choose the reside numbers for the group (like 1 2 3 4):")
        list1=re.split("\s+",string.strip(list1))
        numlist1=[]
        for aa in list1:
            numlist1.append(int(aa))
        if max(numlist1) > reside_list[-1][1]:
            print "base number %d is out of range." %max(numlist1) 
            continue
        else:
            print "you choose base : ",numlist1
            break

    return numlist1

def Get_atom(atom_list,residue_number):
    '''
    choose atoms from a residue.
    '''
    for atom in atom_list:
        if atom.residue_serial==residue_number:
            print "%8d\t %4s" %(atom.atom_serial,atom.atom_name)

    while True:
        list1=raw_input("Choose the atom number for fitting:")
        list1=re.split("\s+",string.strip(list1))
        numlist1=[]
        for aa in list1:
            numlist1.append(int(aa))  
        print "you choose base : ",numlist1
        break

    return numlist1


def Save_file(filename,atom_list):
    ''' 
    Save the result to a structure file like pdb or gro.
    '''
    if os.path.isfile(filename):
        try:
            os.rename(filename,"#"+filename+"#")
            print "backup file %s to %s" %(filename,"#"+filename+"#")
        except OSError,e:
            print e
    try:
        fp=open(filename,'w+')
    except:
        print "except in opening %s" %filename
        return False

    if filename.endswith(".pdb"):
        fp.write("TITLE    PDB file created by Modify_Coor.py\n")
        for a in atom_list:
            fp.write(a.atom_2_PDBformat()+"\n")
        fp.write("END\n")
    elif filename.endswith(".gro"):
        fp.write("%d\n" %len(atom_list))
        for a in atom_list:
            fp.write(a.atom_2_GROformat()+"\n")
    elif filename.endswith(".pqr"):
        if isinstance(atom_list[0],PQR.PQR):
            for a in atom_list:
                fp.write(a.atom_2_PQRformat()+"\n")
        else:
            print "Error:"
    else:
        print "The coordinate file : %s is invilad, both *.gro and *.pdb are allowd." %filename
        return False

    print "Save to the file %s" %filename
    return True


