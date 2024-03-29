#-------------------------------------------------------------------------
# Compute absorption cross-section from Reimers INDO/CI oscillator strengths
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
absorption.py -i <inputfile> -o <outputfile> -e <energy> -s <step> -g <gamma> 

Required:
    -i    Base file name for input file
          NOTE: The script automatically attempts to remove file extensions, 
                starting with the final "." in the file name. If your file name 
                includes periods besides the one before the extension, use the
                full file name (including the extension) here. If not, either 
                the full file name or the base without the extension works. 
                The script will try to append a .out extension for the
                MOPAC output file.

Optional:
    -o    Base file name for output file                   Default = matches -i
    -e    Maximum energy of absorption spectrum (eV)       Default = 8.0
    -s    Step size for absorption spectrum (eV)           Default = 0.02
    -g    Gamma (Lorentzian broadening factor)             Default = 0.1088 eV (= 0.004 au)
    
"""

# Parse input options
try:
    options, remainder = getopt.getopt(sys.argv[1:],"hi:o:e:s:g:t:",['--help','--input=','--output=','--energy=','--step=','--gamma=','--type='])
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
    elif opt in ('-g','--gamma'):
        gamma = float(arg)
    elif opt in ('-t','--type'):
        tfilename = arg

types = read_types(tfilename)

out = Molecule(infilename)

# Get info and write preliminary output files
out.read_configs(types)
out.write_orbs(outfilename, types)
out.write_osc(outfilename, types)
#out.write_exc(outfilename, types)
out.write_sigma(outfilename, types, max_e, e_step, gamma)

out.file.close()

