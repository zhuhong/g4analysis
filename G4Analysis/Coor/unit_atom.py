#-*- coding : utf-8 -*-
'''
Create on 2011-08-23
@version: 0.0.1
@change: \n
    - 2011-01-25\n
        - copy from Simple_atom.
'''

import string
import re

class unit_atom():
    '''
    It's a simple atom class. 
    '''
    def __init__(self,atom_name="",atom_serial=0,residue_name="",\
            residue_serial=0,atom_coor_x=0.0,atom_coor_y=0.0,atom_coor_z=0.0):
        self.atom_name=atom_name
        self.atom_serial=atom_serial
        self.residue_name=residue_name
        self.residue_serial=residue_serial
        self.atom_coor_x=atom_coor_x
        self.atom_coor_y=atom_coor_y
        self.atom_coor_z=atom_coor_z

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

    def atom_2_GROformat(self):
        s="%5d%5s%5s%5d%8.3f%8.3f%8.3f" \
                %(self.residue_serial , string.ljust(self.residue_name,5) , self.atom_name , self.atom_serial , \
                self.atom_coor_x , self.atom_coor_y , self.atom_coor_z )
        return s


