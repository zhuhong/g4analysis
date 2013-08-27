#-*- coding : utf-8 -*-

'''
Created on 2011-01-22\n
@version: 0.0.0
@change:\n
    - 2011-01-22.\n
        - The function B{Get_Atom_list()} not finished.\n
    - 2011-01-23.\n
        - Finished the function B{Get_Atom_list()}.\n
    - 2011-02-01.\n
        - Finished the function B{Get_Segment_list()}.\n
    - 2011-02-16.\n
        - Add functions B{atom_2_GROformat()} and B{atom_2_PDBformat()}
    
'''
import unit_atom
import string
import re

class GRO():
    '''
    A class design to read gro file.
    '''
    def __init__(self,residue_name=" ",residue_serial = 0 , atom_name=" ",atom_serial = 0,
            atom_coor_x = 0.0 , atom_coor_y = 0.0 , atom_coor_z = 0.0,
            atom_velo_x = 0.0 , atom_velo_y = 0.0 , atom_velo_z = 0.0):
        self.residue_name=residue_name
        self.residue_serial=residue_serial
        self.atom_name=atom_name
        self.atom_serial=atom_serial
        self.atom_coor_x=atom_coor_x
        self.atom_coor_y=atom_coor_y
        self.atom_coor_z=atom_coor_z
        self.atom_velo_x=atom_velo_x
        self.atom_velo_y=atom_velo_y
        self.atom_velo_z=atom_velo_z

    def atom_2_GROformat(self):
        '''
        return a string in gro fromat.
        '''
        s="%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f" \
                %(self.residue_serial , string.ljust(self.residue_name,5) , self.atom_name , self.atom_serial , \
                self.atom_coor_x , self.atom_coor_y , self.atom_coor_z , self.atom_velo_x , \
                self.atom_velo_y , self.atom_velo_z)
        return s

    def atom_2_PDBformat(self):
        temp_atom=""
        if re.match('^\d',self.atom_name) != None:
            temp_atom=self.atom_name.ljust(4)
        else:
            if len(self.atom_name)<4:
                temp_atom=(" "+self.atom_name).ljust(4)
            else:
                temp_atom=self.atom_name

        s = "%s%5d %s %3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f "  \
                  % ("ATOM".ljust(6) , self.atom_serial , temp_atom , self.residue_name.rjust(3) , \
                  self.residue_serial , self.atom_coor_x*10 , self.atom_coor_y*10 , self.atom_coor_z*10,\
                  1.00,0.00)
        return s

def Get_Atom_list(filename):
    '''
    Read in a gro file and return a GRO class list. 
    '''
    atom_list=list()
    try:
        fp=open(filename,'r')
    except:
        print "Error: No such file: "+filename 
        return atom_list
    allLines = fp.readlines() 
    gro_name=allLines[0]
    gro_atom_number=int(allLines[1])
    gro_box_size=re.split(string.strip(allLines[-1]),"\s")
    for i in range(len(allLines)-3):
        atom_temp=GRO()
        try:
            atom_temp.residue_serial=int(string.strip(allLines[i+2][0:5]))
        except:
            pass
        atom_temp.residue_name=string.strip(allLines[i+2][5:10])
        atom_temp.atom_name=string.strip(allLines[i+2][10:15])
        try:
            atom_temp.atom_serial=int(string.strip(allLines[i+2][15:20]))
        except:
            pass
        try:
            atom_temp.atom_coor_x=float(string.strip(allLines[i+2][20:28]))
            atom_temp.atom_coor_y=float(string.strip(allLines[i+2][28:36]))
            atom_temp.atom_coor_z=float(string.strip(allLines[i+2][36:44]))

            if len(allLines[i+2]) > 50:
                atom_temp.atom_velo_x=float(string.strip(allLines[i+2][44:52]))
                atom_temp.atom_velo_y=float(string.strip(allLines[i+2][52:60]))
                atom_temp.atom_velo_z=float(string.strip(allLines[i+2][60:68]))
        except:
            pass
        atom_list.append(atom_temp)

    if len(atom_list) != gro_atom_number :
        print "Warning: the atom number in gro file %d not match the gro file title %d." %(len(atom_list),gro_atom_number)

    return atom_list

def Read_GRO_2_SimpleAtom_list(filename):
    '''
    Read in a gro file and return a GRO class list. 
    '''
    try:
        fp=open(filename,'r')
    except:
        print "Error: No such file: "+filename 
        exit(1) 
    allLines = fp.readlines() 
    atom_list=[]
    gro_name=allLines[0]
    gro_atom_number=int(allLines[1])
    gro_box_size=re.split(string.strip(allLines[-1]),"\s")
    for i in range(len(allLines)-3):
        atom_temp=unit_atom.unit_atom()
        try:
            atom_temp.residue_serial=int(string.strip(allLines[i+2][0:5]))
        except:
            pass
        atom_temp.residue_name=string.strip(allLines[i+2][5:10])
        atom_temp.atom_name=string.strip(allLines[i+2][10:15])
        try:
            atom_temp.atom_serial=int(string.strip(allLines[i+2][15:20]))
        except:
            pass
        try:
            atom_temp.atom_coor_x=float(string.strip(allLines[i+2][20:28]))
            atom_temp.atom_coor_y=float(string.strip(allLines[i+2][28:36]))
            atom_temp.atom_coor_z=float(string.strip(allLines[i+2][36:44]))
        except:
            pass
        atom_list.append(atom_temp)

    if len(atom_list) != gro_atom_number :
        print "Warning: the atom number in gro file %d not match the gro file title %d." %(len(atom_list),gro_atom_number)

    return atom_list


