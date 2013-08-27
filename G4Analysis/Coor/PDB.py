#-*- coding: utf-8 -*-

'''
Created on 2011-1-1\n
@author: zhuh
@version: 0.1.0
@change: \n
    - 2011-1-10 \n
        - finished the B{get_atom_list()} function.
    - 2011-01-18\n
        - Add the function B{Get_Atom_in_residue()}.
    - 2011-02-16\n
        - Add the function B{atom_2_PDBformat()} and B{atom_2_GROformat()}
    - 2011-04-21\n
        - Modified a bug.
'''
import unit_atom
import string
import re

class Atom_class:
    '''
    The core class, which contain all the atom pramaters in PDB file.
    It's can be initialized with all the paramaters or part.
    '''
    
    def __init__(self , character = "" , atom_serial = 0 , atom_name = "" ,
            alter_local_indicater = "" , residue_name = "" , chain_indicater = "" ,
            residue_serial = 0 , code_for_insertions_of_residues = "" ,
            atom_coor_x = 0.0 , atom_coor_y = 0.0 , atom_coor_z = 0.0 , occupancy = 0.0 , 
            temp_factor = 0 , segment_indent = "" , element_symbol = "" , 
            charge = ""):
        self.character = character
        self.atom_serial = atom_serial
        self.atom_name = atom_name
        self.alter_local_indicater = alter_local_indicater
        self.residue_name = residue_name
        self.chain_indicater = chain_indicater
        self.residue_serial = residue_serial
        self.code_for_insertions_of_residues = code_for_insertions_of_residues
        self.atom_coor_x = atom_coor_x
        self.atom_coor_y = atom_coor_y
        self.atom_coor_z = atom_coor_z
        self.occupancy = occupancy
        self.temp_factor = temp_factor
        self.segment_indent = segment_indent
        self.element_symbol = element_symbol
        self.charge = charge
              
    def atom_2_PDBformat(self):
        '''得到一个格式化的Atom sting'''
        temp_atom=""
        if re.match('^\d',self.atom_name) != None:
            temp_atom=self.atom_name.ljust(4)
        else:
            if len(self.atom_name) <4:
                temp_atom=(" "+self.atom_name).ljust(4)
            else:
                temp_atom=self.atom_name

        s = "%s%5d %s %3s %1s%4d%s    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s" \
                % (self.character.ljust(6) , self.atom_serial , temp_atom,  self.residue_name.rjust(3) , \
                self.chain_indicater , self.residue_serial , self.code_for_insertions_of_residues , \
                self.atom_coor_x*10 , self.atom_coor_y*10 , self.atom_coor_z*10 , self.occupancy ,\
                self.temp_factor , self.segment_indent.ljust(4) , \
                self.element_symbol.rjust(2) , self.charge)
        return s

    def atom_2_GROformat(self):
        s="%5d%5s%5s%5d%8.3f%8.3f%8.3f" \
                %(self.residue_serial , string.ljust(self.residue_name,5) , self.atom_name , self.atom_serial , \
                self.atom_coor_x , self.atom_coor_y , self.atom_coor_z )
        return s

def Get_Atom_list(filename):
    '''把每个atom的信息写到一个Atom_class的类中，然后输出一个Atom_class的list , 
    部分异常处理未完成。'''
    alist = list()
    try:
        fp = open(filename , 'r')
    except:
        print "Error: No such file: "+filename
        return alist
    allLines = fp.readlines()
    for eachline in allLines:
        if ( eachline.startswith("ATOM") or eachline.startswith("HETATM") ):
            temp_Atom_class = Atom_class()
            temp_Atom_class.character = string.strip(eachline[0:6])
            try:
                temp_Atom_class.atom_serial = int(string.strip(eachline[6:11]))
            except:
                pass
            temp_Atom_class.atom_name = string.strip(eachline[12:16])
            temp_Atom_class.alter_local_indicater = string.strip(eachline[16])
            temp_Atom_class.residue_name = string.strip(eachline[17:20])
            temp_Atom_class.chain_indicater = string.strip(eachline[21])
            try:
                temp_Atom_class.residue_serial = int(string.strip(eachline[22:26]))
            except:
                pass
            temp_Atom_class.code_for_insertions_of_residues = string.strip(eachline[27])
            try:
                temp_Atom_class.atom_coor_x = float(string.strip(eachline[30:38]))/10
                temp_Atom_class.atom_coor_y = float(string.strip(eachline[38:46]))/10
                temp_Atom_class.atom_coor_z = float(string.strip(eachline[46:54]))/10
                if len(eachline) > 59:
                    temp_Atom_class.occupancy = float(string.strip(eachline[54:60]))
                if len(eachline) > 65:
                    temp_Atom_class.temp_factor = float(string.strip(eachline[60:66]))
            except:
                pass
            if len(eachline) > 75:
                temp_Atom_class.segment_indent = string.strip(eachline[72:76])
            if len(eachline) > 77:
                temp_Atom_class.element_symbol = string.strip(eachline[76:78])
            if len(eachline) > 79:
                temp_Atom_class.charge = string.strip(eachline[78:80])
            
            alist.append(temp_Atom_class)
            
    return alist

def Read_PDB_2_SimpleAtom_list(filename):
    '''
    把每个atom的信息写到一个Simple_atom的类中，然后输出一个Atom_class的list , 
    部分异常处理未完成。
    '''
    try:
        fp = open(filename , 'r')
    except:
        print "Error: No such file: "+filename
        exit(1)
    allLines = fp.readlines()
    list = []
    for eachline in allLines:
        if ("ATOM" in eachline) or ("HETEAM" in eachline):
            temp_Atom = unit_atom.unit_atom()
            try:
                temp_Atom.atom_serial = int(string.strip(eachline[6:11]))
            except:
                pass
            temp_Atom.atom_name = string.strip(eachline[12:16])
            temp_Atom.residue_name = string.strip(eachline[17:20])
            try:
                temp_Atom.residue_serial = int(string.strip(eachline[22:26]))
            except:
                pass
            try:
                temp_Atom.coor_x = float(string.strip(eachline[30:38]))
                temp_Atom.coor_y = float(string.strip(eachline[38:46]))
                temp_Atom.coor_z = float(string.strip(eachline[46:54]))
            except:
                pass
            list.append(temp_Atom)
            
    return list
