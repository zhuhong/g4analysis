g4analysis
==========


---------------
Introduction
---------------


A python package for DNA structural analysis


------------
Usage
------------

**Files**


Option|Type|Filename|Description
------|----|--------|-----------
-p       | Input |  coor_file   | Structure fie: gro pdb etc.
-f       | Input |  traj_file   | Trajectory: xtc trr.
-o       | Input |  output_file | xvgr/xmgr file.
-i       | Input |  para_an.in  | input parmarter file.
------|----|--------|-----------

**Other options**

==========    ======    ===========  ==============================================================
Option          Type    Value        Description
==========    ======    ===========  ==============================================================
--helical       bool    False 	     Calculate the helical parameters of nucleic acids.
--dihedral      bool    False 	     Calculate the backbone dihedral parameters of nucleic acids.
--rise          bool    False        Calculate the distance of DNA bases groups.
--twist         bool    False        Calculate the twist of DNA bases groups.
--rmsd          bool    False        skip Calculate the RMSD of DNA bases groups.
--begin         int     0            First frame (ps) to read from trajectory.
--end           int     -1           Last frame (ps) to read from trajectory.
--skip          int     1            Get frames when frame MOD skip = 0
-h              bool    yes          Print help info and quit
==========    ======    ===========  ==============================================================

Two modes to use this program.

Interactive Mode.
-----------------

	G4Analysis -p coor_file -f traj_file -o output_file --rise/twist/rmsd [--begin/end/skip]

Input File Mode.
-----------------

	G4Analysis -p coor_file -f traj_file -i input_file [--begin/end/skip]

para.in format:

    group_1(ID1:ID2:...) group_2(ID1:ID2:...) result.xvg




version 0.2.0
---------------
* Rewrited Simple_atom and modified some bugs.
* Add G4_area.py to calculate the area of four O6 atoms in the G-quartet.



version 0.1.0
----------------
* For B-DNA. helical parameters and dihedral parameters.
* For G-quadruplex, rise, twist and RMSD_z for G-quartets.
