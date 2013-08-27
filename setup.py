
#from setuptools import setup
from distutils.core import setup
import sys

try:
    import numpy
except ImportError:
    print "*** package 'numpy' not found ***"
    print "Please get it from http://numpy.scipy.org/ or install it through your package manager."
    sys.exit(-1)



try:
    import MDAnalysis
except ImportError:
    print "*** package 'MDAnalysis' not found ***"
    sys.exit(-1)


setup(
        name='G4Analysis',
        version='0.1.0',
        author="zhuh",
        packages=['G4Analysis','G4Analysis.Coor','G4Analysis.mymath'],
        scripts=['scripts/G4analysis.py']
        )

