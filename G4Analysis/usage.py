'''
Created on 2011-03-05
@change:\n
    - 2011-04-02.\n
        - add function echo
'''
import sys

version="0.1.0"

def File_input():
    print "Option \t %6s  %20s    %10s " %("Type","Filename","Description")

def Coor_file(filename,IO="Input",Option="-p"):
    print "%6s   %6s  %20s \t Structure file: gro pdb etc." %(Option,IO,filename)

def Traj_file(filename,IO="Input",Option="-f"):
    print "%6s   %6s  %20s \t Trajectory: xtc trr. " %(Option,IO,filename)

def Ind_file(filename,IO="Input",Option="-n"):
    print "%6s   %6s  %20s \t Index file." %(Option,IO,filename)
    
def Xvgr_file(filename,IO="Output",Option="-o"):
    print "%6s   %6s  %20s \t xvgr/xmgr file. " %(Option,IO,filename)

def Parm_file(filename,IO="Input",Option="-i"):
    print "%6s   %6s  %20s \t input parmarter file. " %(Option,IO,filename)

def Print_line():
    print "-"*65

def Type_input():
    print "Option         %10s   %6s \t %10s " %("Type","Value","Description")

def Show_help(value,Option="-h"):
    print "%12s   %10s   %6s \t Print help info and quit " %(Option,"bool",value)

def Show_skip(value,Option="--skip"):
    print "%12s   %10s   %6d \t Get frames when frame MOD skip = 0 " %(Option,"int",value)

def Show(Option,Type,value,description):
    print "%12s   %10s   %6s \t %s" %(Option,Type,value,description)

def echo(s=''):
    """Simple string output that immediately prints to the console."""
    sys.stderr.write(s)
    sys.stderr.flush()


