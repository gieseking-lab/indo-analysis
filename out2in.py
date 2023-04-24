#-------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------

import string,sys,math,re
from mol_assign import *
import getopt

prog = ' '
template = 'def_xyz'

helpfile = """
out2in.py -i <inputfile> -o <outputfile> -z <program>

    -i    Base file name for input and output (must be included)
    -o    Base file name for all output files, if different from input

    -t    Template for output                                       Default = def_xyz
             def_xyz     xyz file
             def_mopac   Basic MOPAC file
             filename    Use file as template (use NAT, CHRG, ELEM, XXX, YYY, ZZZ to insert info)

    -z    Program generating output (Mopac, ADF, Qchem, etc.)       Default = Determine automatically
"""

# Parse input options
try:
    options, remainder = getopt.getopt(sys.argv[1:],"hi:o:t:z:",['--help','--input=','--output=','--template=','--program='])
except getopt.GetoptError as err:
    print str(err)
    print helpfile
    sys.exit(2)

if remainder != []:
    print "Error: options not read - %r" % remainder
    print helpfile
    sys.exit(2)


for opt, arg in options:
    if opt in ('-h','--help'):
        print helpfile
        sys.exit(2)
    elif opt in ('-i','--input'):
        infilename = arg
        if '.' in arg:
            arg = arg[:arg.rfind('.')]
        outfilename = arg
    elif opt in ('-o','--output'):
        if '.' in arg:
            arg = arg[:arg.rfind('.')]
        outfilename = arg
    elif opt in ('-t','--template'):
        template = arg
    elif opt in ('-z','--program'):
        prog = arg.lower()

out, prog = open_mol(infilename, prog)

out.write_input(outfilename, template)

out.file.close()

