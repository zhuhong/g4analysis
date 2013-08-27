#!/usr/bin/env python
'''
Create on 2011-01-19\n
@change:
    - 2011-01-20.\n
        - Now it's runing.
    - 2011-01-22.\n
        - Add code for checking the input. Make sure the user input 
        the right size and number of the base for a group.
    2011-02-01.\n
        - Modified this script, So both *.gro and *.pdb are allowd for \
                coordinare input file.
        - Changed the version to 0.2.1\n
    2011-04-26.\n 
        - Modified the usage part.
    2012.02.24.\n
        - Combining the G4_twist part to this script.
'''

import sys
import string
import re
import getopt
import os
from G4Analysis import Simple_atom
from G4Analysis import G4_rise
from G4Analysis import G4_twist
from G4Analysis import usage
from G4Analysis import DNA_analysis
from G4Analysis.Coor import atomlib

def Usage(coor_file="coor_file",traj_file="traj_file",output_file="output_file",\
        parm_file="para_analysis.in",skip=1,show_help="yes",\
        calcu_helical="False", calcu_dihedral="False",\
        calcu_rise="False",calcu_twist="False",calcu_rmsd="False",begin=0,end=-1):
    '''
    print the usage information.
    '''
    print ""
    usage.File_input()
    usage.Print_line()
    usage.Coor_file(coor_file,"Input")
    usage.Traj_file(traj_file,"Input")
    usage.Xvgr_file(output_file,"Input")
    usage.Parm_file(parm_file,"Input")
    usage.Print_line()
    print ""
    usage.Type_input()
    usage.Print_line()
    usage.Show("--helical","bool",calcu_helical,"Calculate the helical parameters of nucleic acids.")
    usage.Show("--dihedral","bool",calcu_dihedral,"Calculate the backbone dihedral parameters of nucleic acids.")
    usage.Show("--rise","bool",calcu_rise,"Calculate the distance of DNA bases groups.")
    usage.Show("--twist","bool",calcu_twist,"Calculate the twist of DNA bases groups.")
    usage.Show("--rmsd","bool",calcu_rmsd,"Calculate the RMSD of DNA bases groups.")
    usage.Show("--begin","int",begin,"First frame (ps) to read from trajectory.")
    usage.Show("--end","int",end,"Last frame (ps) to read from trajectory.")
    usage.Show_skip(skip)
    usage.Show_help(show_help)
    print ""

def Check_argv(argv):
    '''
    Check the sys.argv list. and return a hash contain the input info.
    '''
    coor_file=""
    traj_file=""
    output_file=""
    parm_file=""
    skip=1
    resu_hash={}
    calcu_rmsd=False
    calcu_rise=False
    calcu_twist=False
    calcu_helical=False
    calcu_dihedral=False
    begin=1
    end=-1

    if len(argv)==1:
        Usage()
        sys.exit()

    try:
        opts,args=getopt.getopt(sys.argv[1:],"p:f:o:i:h",["skip=","begin=","end=","rmsd","rise","twist","dihedral","helical"])
    except getopt.GetoptError,e:
        print e
        sys.exit()

    if len(opts)==0:
        print "Use 'G4Analysis.py -h' to see the help information "

    for a,b in opts:
        if a == "-p" :
            coor_file=b
        elif a == "-f" :
            traj_file=b
        elif a == "-o":
            output_file=b
        elif a == "--skip" :
            try:
                skip = int(b)
            except:
                pass

        elif a == "--rmsd":
            calcu_rmsd=True

        elif a=="--rise":
            calcu_rise=True

        elif a=="--twist":
            calcu_twist=True

        elif a=="--dihedral":
            calcu_dihedral=True

        elif a=="--helical":
            calcu_helical=True

        elif a=="-i":
            parm_file=b

        elif a=="-h":
            Usage()
            sys.exit()
        elif a=="--begin":
            try:
                begin=int(b)
            except:
                pass
        elif a=="--end":
            try:
                end=int(b)
            except:
                pass

    if os.path.isfile(coor_file):
        resu_hash["coor_file"]=coor_file
        resu_hash["traj_file"]=traj_file
        resu_hash["parm_file"]=parm_file
        resu_hash["output_file"]=output_file

#    if output_file=="" and parm_file=="":
#        print "Error: No file for output."
#        sys.exit()

    resu_hash["skip"]=skip
    resu_hash["calcu_rmsd"]=calcu_rmsd
    resu_hash["calcu_rise"]=calcu_rise
    resu_hash["calcu_twist"]=calcu_twist
    resu_hash["calcu_dihedral"]=calcu_dihedral
    resu_hash["calcu_helical"]=calcu_helical
    resu_hash["begin"]=begin
    resu_hash["end"]=end

    Usage(coor_file,traj_file,output_file,parm_file,skip,"no",\
            calcu_helical,calcu_dihedral,\
            calcu_rise,calcu_twist,calcu_rmsd,begin,end)

    return resu_hash


def Print_ProgInfo():
    '''
    Print program information. the version, some notice and tips.
    '''
    print " "
    print "  "*13,"-)","G4_analysis","(-"
    print " "
    print "  "*12,"-)"," Version: %s " %usage.version ,"(-" 
    print " "
    Print_Description()
    print " Note:You should choose 1, 2 or 4 bases for a group."
    print " "

def Print_Description():
    '''
    Print the description of this script.
    '''
    print "DESCRIPTION"
    print "-----------"
    a='''G4_analysis.py is used for calculate the rise and twist of two G-quartet layers.
Usage:
    interactive format:
    using -o result.xvg
    script format:
    using -i parm.in
para.in format:
    group_1(ID1:ID2:...) group_2(ID1:ID2:...) result.xvg
    '''
    print a
    print ""

if __name__=="__main__":
    Print_ProgInfo()
    argc=len(sys.argv)
    resu=Check_argv(sys.argv)
    is_get_dt=False
    have_parm_file=os.path.isfile(resu["parm_file"])

    list_group_1=list()
    list_group_2=list()
    list_output=list()


# step 1, get the dt information of the trajectory.
    if resu["traj_file"].endswith("mdcrd") or resu["traj_file"].endswith("dcd"):
        dt=raw_input("Input the time step between frames for Amber trajectory file (ps).")
        is_get_dt=True

# step 2, calculatin the rise paramters.
    if resu["calcu_rise"]==True:
        if resu["traj_file"]=="":
            l1=Simple_atom.Get_residue(resu["coor_file"],True)
            l2=Simple_atom.Get_residue(resu["coor_file"],False)

            G4_rise.Get_parallel_fromTOP(resu["coor_file"],l1,l2)
            sys.exit()
            
        if have_parm_file:
            fp=open(resu["parm_file"])
            lines=fp.readlines()
            for line in lines:
                try:
                    [group1,group2,outputname]=line.split()
                except ValueError,e:
                    print e
                    print "Warning: The line [%s] was not standard, and this line will be ignored." %line
                    continue
                g1=group1.split(":")
                g2=group2.split(":")
                gg1=list()
                gg2=list()
                gg1=([int(i) for i in g1])
                gg2=([int(i) for i in g2])
                list_group_1.append(gg1)
                list_group_2.append(gg2)
                list_output.append(outputname)

        else:
            l1=Simple_atom.Get_residue(resu["coor_file"],True)
            l2=Simple_atom.Get_residue(resu["coor_file"],False)

            list_group_1.append(l1)
            list_group_2.append(l2)
            list_output.append(resu["output_file"])

        if is_get_dt:
            G4_rise.Get_rise_fromTRJ(resu["traj_file"],resu["coor_file"],\
                    list_group_1,list_group_2,list_output,resu["skip"],float(dt),\
                    resu["begin"],resu["end"])
        else:
            G4_rise.Get_rise_fromTRJ(resu["traj_file"],resu["coor_file"],\
                    list_group_1,list_group_2,list_output,resu["skip"],\
                    begin=resu["begin"],end=resu["end"])

# step 3, calculating the twise paramters.
    if resu["calcu_twist"]==True:
        if have_parm_file:
            fp=open(resu["parm_file"])
            lines=fp.readlines()
            for line in lines:
                try:
                    [group1,group2,outputname]=line.split()
                except ValueError,e:
                    print e
                    print "I guess you forget the --rmsd for calculating the rmsd."
                g1=group1.split(":")
                g2=group2.split(":")
                gg1=list()
                gg2=list()
                gg1=([int(i) for i in g1])
                gg2=([int(i) for i in g2])
                list_group_1.append(gg1)
                list_group_2.append(gg2)
                list_output.append(outputname)

        else:
            l1=Simple_atom.Get_residue(resu["coor_file"],True)
            l2=Simple_atom.Get_residue(resu["coor_file"],False)

            list_group_1.append(l1)
            list_group_2.append(l2)
            list_output.append(resu["output_file"])

        if is_get_dt:
            G4_twist.Get_twist_in_GDNA2(resu["traj_file"],resu["coor_file"],\
                    list_group_1,list_group_2,list_output,resu["skip"],float(dt),\
                    resu["begin"],resu["end"])
        else:
            G4_twist.Get_twist_in_GDNA2(resu["traj_file"],resu["coor_file"],\
                    list_group_1,list_group_2,list_output,resu["skip"],\
                    begin=resu["begin"],end=resu["end"])

#step 4, calculating the rmsd parameters.
    if resu["calcu_rmsd"]==True:
        if have_parm_file:
            fp=open(resu["parm_file"])
            lines=fp.readlines()
            for line in lines:
                [group1,outputname]=line.split()
                g1=group1.split(":")
                list_group_1.append([int(i) for i in g1])
                list_output.append(outputname)

        else:
            l1=Simple_atom.Get_residue(resu["coor_file"],True)
            list_group_1.append(l1)
            list_output.append(resu["output_file"])

        if is_get_dt:
            G4_rise.Get_RMSD_fromTRJ(resu["traj_file"],resu["coor_file"],\
                    list_group_1,list_output,resu["skip"],float(dt),\
                    resu["begin"],resu["end"])
        else:
            print list_output
            G4_rise.Get_RMSD_fromTRJ(resu["traj_file"],resu["coor_file"],\
                    list_group_1,list_output,resu["skip"],\
                    begin=resu["begin"],end=resu["end"])

# step 5, calculating the helical parameters.
    if resu["calcu_helical"]==True:
        if have_parm_file:
            pass

# step 6, calculating the dihedral parameters.
    if resu["calcu_dihedral"]==True:
        if resu["traj_file"]=="":
            atom_list=Simple_atom.Get_Simple_atom_list(resu["coor_file"])
            reside_list=Simple_atom.Get_Residue_list(atom_list)
            chain=list()
            for residue in reside_list:
                if residue[0] in atomlib.RESIDUE_NAME_LIST:
                    chain.append(residue[1])
                else:
                    pass

            DNA_analysis.Get_Dihedral_fromTOP(resu["coor_file"],chain,True)
            sys.exit(0)

        if have_parm_file:
            fp=open(resu["parm_file"])
            lines=fp.readlines()
            for line in lines:
                [group1,outputname]=line.split()
                list_group_1.append(int(group1))
                list_output.append(outputname)

        else:
            l1=Simple_atom.Get_residue(resu["coor_file"],True)
            list_group_1=l1
            list_output.append(resu["output_file"])

        if is_get_dt:
            DNA_analysis.Get_Dihedral_fromTRJ(resu["traj_file"],resu["coor_file"],\
                    list_group_1, list_output,resu["skip"],float(dt),resu["begin"],resu["end"])

        else:
            DNA_analysis.Get_Dihedral_fromTRJ(resu["traj_file"],resu["coor_file"],\
                    list_group_1, list_output,resu["skip"],begin=resu["begin"],end=resu["end"])


