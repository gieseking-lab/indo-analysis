#!/usr/bin/python

#-------------------------------------------------------------------------
# Compute absorption cross-section from Reimers INDO/CIS oscillator strengths
#
# Version 2.0
#
# Rebecca Gieseking, 12/2/2016
# 
#-------------------------------------------------------------------------

import sys
from utils.molecule import Molecule, read_types
import getopt

# Constants for later
gamma = 0.1088
max_e = 8.0
e_step = 0.02
tfilename = 'types'
types = []
prog = ' '

helpfile = """
indo_sigma.py -i <inputfile> -o <outputfile> -e <energy> -s <step> -t <type_file> 

    -i    Base file name for input and output (must be included)
    -o    Base file name for all output files, if different from input

    -e    Maximum energy of the absorption spectrum (eV)        Default = 8.0
    -s    Step size for the absorption spectrum (eV)            Default = 0.02

    -t    File that includes a list of the types of absorption spectra to compute
          Format:
             Atom    1 20             From atoms 1-20 to 21-max
             Attype  Ag               From Ag atoms to all other atoms
             Orbtype Ag D             From Ag D to all others
             Orb     D                From all D orbitals to all others
"""

# Parse input options
try:
    options, remainder = getopt.getopt(sys.argv[1:],"hi:o:e:s:t:",['--help','--input=','--output=','--energy=','--step=','--type='])
except getopt.GetoptError as err:
    print(str(err))
    print(helpfile)
    sys.exit(2)

if remainder != []:
    print("Error: options not read - %r" % remainder)
    print(helpfile)
    sys.exit(2)


for opt, arg in options:
    if opt in ('-h','--help'):
        print(helpfile)
        sys.exit(2)
    elif opt in ('-i','--input'):
        if '.' in arg:
            arg = arg[:arg.rfind('.')]
        infilename = arg
        outfilename = arg
    elif opt in ('-o','--output'):
        if '.' in arg:
            arg = arg[:arg.rfind('.')]
        outfilename = arg
    elif opt in ('-e','--energy'):
        max_e = float(arg)
    elif opt in ('-s','--step'):
        e_step = float(arg)
    elif opt in ('-t','--type'):
        tfilename = arg

types = read_types(tfilename)

out = Molecule(infilename)

# Get info and write preliminary output files
out.read_configs(types)
out.write_exc(outfilename, types)

out.file.close()

