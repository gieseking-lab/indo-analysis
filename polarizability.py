#-------------------------------------------------------------------------
# Compute INDO/CI sum-over-states polarizability
#-------------------------------------------------------------------------

import sys
from utils.molecule import Molecule
import getopt

# Initialize defaults
omega = 0.0
gamma = 0.1088j
states = 0

helpfile = """
polarizability.py -i <inputfile> -o <outputfile> -w <omega> -g <gamma> -s <# states>

Required:
    -i    Base file name for input file
          NOTE: The script automatically attempts to remove file extensions, 
                starting with the final "." in the file name. If your file name 
                includes periods besides the one before the extension, use the
                full file name (including the extension) here. If not, either 
                the full file name or the base without the extension works. 
                The script will try to append a .log or .out extension for the
                MOPAC output file.

Optional:
    -o    Base file name for output file                        Default = matches -i
    -w    Energy at which polarizabilities are computed (eV)    Default = 0.0
    -g    Lifetime (broadening) of excited states (eV)          Default = 0.1088i (= 0.004 a.u.)
    -s    Number of states to include in SOS expression         Default = All states
"""

# Parse input options
try:
    options, remainder = getopt.getopt(sys.argv[1:],"hi:o:w:g:s:",['--help','--input=','--output=','--omega=','--gamma=','--states='])
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
        infilename = arg
        if '.' in arg:
            arg = arg[:arg.rfind('.')]
        #infilename = arg
        outfilename = arg
    elif opt in ('-o','--output'):
        if '.' in arg:
            arg = arg[:arg.rfind('.')]
        outfilename = arg
    elif opt in ('-w','--omega'):
        omega = float(arg)
    elif opt in ('-g','--gamma'):
        gamma = complex(0.0,float(arg))
    elif opt in ('-s','--states'):
        states = int(arg)

out = Molecule(infilename)

out.read_states()
out.write_alpha(outfilename, omega, gamma, states)

out.file.close()



