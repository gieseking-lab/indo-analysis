#-------------------------------------------------------------------------
# Calculate and print the electron density
# Options for orbital density, excitation/state transition density, or
# excitation/state change in density
#-------------------------------------------------------------------------

import sys
from utils.molecule import Molecule
import getopt
import utils.cubeutils as cubeutils
import utils.constants as cnst
import utils.param as param

#########################
# Start of main program #
#########################

# Set up constants and initial parameters
gap = 0.25/cnst.bohr2ang
extra = 3.0/cnst.bohr2ang
dtype = 0
amin = 1
amax = 1
pfile = None

helpfile = """
cubegen.py -i <inputfile> -o <outputfile> -t <type> -m <min> -x <max> -p <paramfile>
           -v <voxelsize> -e <extra>

Required:
    -i    Base file name for input file
          NOTE: The script automatically attempts to remove file extensions, 
                starting with the final "." in the file name. If your file name 
                includes periods besides the one before the extension, use the
                full file name (including the extension) here. If not, either 
                the full file name or the base without the extension works. 
                The script will try to append a .log or .out extension for the
                MOPAC output file.
    -t  Type of cube file to calculate
        Options:
            orbital (or 0)       Wavefunction of a molecular orbital
            ctrans  (or 1)       Transition density for a CI configuration
            config  (or 2)       Change in electron density upon excitation from 
                                  the SCF ground state to a CI configuration
            trans   (or 3)       Transition density for an excited state 
            exc     (or 4)       Change in electron density upon excitation from 
                                  the SCF ground state to an excited state
        NOTES: 
        The 'orbital' option is implemented for all methods in MOPAC. All other
          options are implemented only for INDO/CI.
        Min/max orbital numbers may be selected either purely by number or in
          reference to HOMO and LUMO (syntax: H-1, L+1, etc.)
        INDO CI configurations are listed in the MOPAC output immediately after 
          the keyword "CI excitations". The first column is the configuration 
          number, and the last few columns indicate which MOs are involved.
        INDO excited states are listed immediately after the keyword "CI trans."
          State 1 is the ground state. 
                                  
Optional:
    -o    Base file name for output file                   Default = matches -i
    -m    Minimum state/orbital to compute                 Default = 1
    -x    Maximum state/orbital to compute                 Default = matches -m
    -p    Parameter file for non-default parameters        Default = None
    -v    Size of voxels, in Angstroms                     Default = 0.25 A
    -e    Size of extra space surrounding the molecule     Default = 3.00 A
"""

# Parse input options
try:
    options, remainder = getopt.getopt(sys.argv[1:],"hi:o:t:m:x:p:v:e:",
                ['--help','--input=','--output=','--type=',"--min=","--max=","--param=","--voxelsize=","--extra="])
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
    elif opt in ('-t','--type'):
        if 'orb' in arg or '0' in arg:
            dtype = 0
        elif 'ct' in arg or '1' in arg:
            dtype = 1
        elif 'con' in arg or '2' in arg:
            dtype = 2
        elif 'tr' in arg or '3' in arg:
            dtype = 3
        elif 'ex' in arg or '4' in arg:
            dtype = 4
    elif opt in ('-m','--min'):
        amin = arg
    elif opt in ('-x','--max'):
        amax = arg
    elif opt in ('-p','--param'):
        pfile = arg
    elif opt in ('-v','--voxelsize'):
        gap = float(arg)/cnst.bohr2ang
    elif opt in ('-e','--extra'):
        extra = float(arg)/cnst.bohr2ang

# Open input and output files
out = Molecule(infilename)
out.read_method()
out.read_atoms()
par = param.Param(pfile)
par.get_orb_params(out.at_types, out.method)
out.read_norbs()
out.read_orbs()

# Convert MO input to MO numbers
if dtype == 0:
    try:
        mo_min = int(amin)
    except ValueError:
        mo_min_inp = amin.upper()
        if 'H' in mo_min_inp:
            mo_min = int(amin[1:]) + out.homo - 1
        else:
            mo_min = int(amin[1:]) + out.homo
    try:
        mo_max = int(amax)
    except ValueError:
        mo_max_inp = amax.upper()
        if 'H' in mo_max_inp:
            mo_max = int(amax[1:]) + out.homo - 1
        else:
            mo_max = int(amax[1:]) + out.homo
    if mo_max < mo_min:
        mo_max = mo_min
elif dtype < 3:
    out.read_configs()
    config_min = int(amin)
    config_max = int(amax)
    if config_max < config_min:
        config_max = config_min
else:
    out.read_configs()
    out.read_states()
    ex_min = int(amin)
    ex_max = int(amax)
    if ex_max < ex_min:
        ex_max = ex_min


# Find the number of voxels in each dimension and offset the molecule appropriately
nvox, minxyz = cubeutils.get_vnum(out.atoms, gap, extra)

print('Number of voxels: ',nvox)

# Build and print densities
#   dtype=0    orbital   Compute density of an orbital
#   dtype=1    ctrans    Compute a configuration transition density
#   dtype=2    config    Compute the change in density for a particular CI configuration
#   dtype=3    trans     Compute an excited state transition density
#   dtype=4    exc       Compute the change in density upon excitation to an excited state
if dtype == 0:
    for mo_curr in range(mo_min-1,mo_max):
        vox = cubeutils.gen_cube_orb(mo_curr, out.at_types, out.atoms, out.mos, 
                                     par.params, extra, nvox, minxyz, gap)
        cubeutils.write_cub(outfilename, dtype, mo_curr+1, vox, out.atoms, out.at_types, 
                            par.params, nvox, minxyz, gap)
elif dtype == 1:
    for config_curr in range(config_min,config_max+1):
        vox = cubeutils.gen_cube_ctrans(config_curr, out.at_types, out.atoms, out.mos, out.configs, 
                                        par.params, extra, nvox, minxyz, gap)
        cubeutils.write_cub(outfilename, dtype, config_curr, vox, out.atoms, out.at_types, 
                            par.params, nvox, minxyz, gap)

elif dtype == 2:
    for config_curr in range(config_min,config_max+1):
        vox = cubeutils.gen_cube_config(config_curr, out.at_types, out.atoms, out.mos, out.configs, 
                                        par.params, extra, nvox, minxyz, gap)
        cubeutils.write_cub(outfilename, dtype, config_curr, vox, out.atoms, out.at_types, 
                            par.params, nvox, minxyz, gap)
elif dtype == 3:
    for ex_curr in range(ex_min,ex_max+1):
        if ex_curr == 1:
            print("Error: Transition density from the ground state (state 1) to state 1 is not defined")
            print("For the first excited state, use state 2")
        else:
            vox = cubeutils.gen_cube_trans(ex_curr, out.at_types, out.atoms, out.mos, out.configs, 
                                           out.states, par.params, extra, nvox, minxyz, gap)
            cubeutils.write_cub(outfilename, dtype, ex_curr, vox, out.atoms, out.at_types, 
                                par.params, nvox, minxyz, gap)
else:
    for ex_curr in range(ex_min,ex_max+1):
        if ex_curr == 1:
            print("Error: Change in density from the ground state (state 1) to state 1 is zero")
            print("For the first excited state, use state 2")
        else:
            vox = cubeutils.gen_cube_exc(ex_curr, out.at_types, out.atoms, out.mos, out.configs, 
                                           out.states, par.params, extra, nvox, minxyz, gap)
            cubeutils.write_cub(outfilename, dtype, ex_curr, vox, out.atoms, out.at_types, 
                                par.params, nvox, minxyz, gap)



