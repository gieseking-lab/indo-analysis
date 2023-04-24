#-------------------------------------------------------------------------
# Compute Raman cross-sections from changes to the INDO/CI static 
# polarizabilities with geometric displacement along normal modes 
# from another level of theory
#
# By default, the Raman intensities are computed using all excited
# states in the INDO/CI output.
#
# Alternatively, the Raman intensities may be
# 
#-------------------------------------------------------------------------

import string,sys,math
from utils.molecule import Molecule
from utils.vibrations import VibAll
import getopt

# Constants for later
types = []
omega = 0.0
gamma = 0.1088j
disp_str = '0.01'
disp = float(disp_str)
coord = ''
nstates = 0
nstatesmin, nstatesmax, nstatesstep = None, None, None
avg = False
width = 20.0
scale = 1E+34

vmin = 350
vmax = 2500

helpfile = """
raman_avg.py -i <inputfile> -o <outputfile> -e <energy> -g <gamma> -d <displacement> -c <coordinate> -n <nrsfile> -z <program>

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
    -o    Base file name for output file                          Default = matches -i
    -e    Energy at which Raman intensities are computed (eV)     Default = 0.0
    -g    Lifetime (broadening) for polarizabilities (eV)         Default = 0.1088 (= 0.004 a.u.)
    -d    Displacement of geometries along vibrational modes (A)  Default = 0.01
    -c    Coordinate along which to compute Raman intensities     Default = isotropic (alternatives are x, y, z)
    -n    Number of states to include in SOS expression           Default = All states
    -f    Minimum vibrational frequency of modes to use (cm-1)    Default = 350
    -v    Maximum vibrational frequency of modes to use (cm-1)    Default = 2500

To compute the Raman intensities averaged over a series of numbers of excited states
instead of one number, use the following options. If one is used, all must be used.
This option is sometimes useful for large systems where the Raman intensities are 
not fully converged with respect to the number of states in the SOS expression.
    -m    Minimum number of excited states to include in SOS expression
    -x    Maximum number of excited states to include in SOS expression
    -s    Step size in number of excited states to include in SOS expression
"""

# Parse input options
try:
    options, remainder = getopt.getopt(sys.argv[1:],"hi:o:e:g:d:c:f:v:m:x:s:",[
        '--help','--input=','--output=','--energy=','--gamma=',
        '--disp=','--coord=','--nstates=','--freqmin=', '--freqmax=',
        '--nstatesmin=','--nstatesmax=','--nstatesstep='])
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
        omega = float(arg)
    elif opt in ('-g','--gamma'):
        gamma = complex(0.0,float(arg))
    elif opt in ('-d','--disp'):
        disp = float(arg)
    elif opt in ('-c','--coord'):
        coord = arg
    elif opt in ('-n','--nstates'):
        nstates = int(arg)
    elif opt in ('-f','--freqmin'):
        vmin = float(arg)
    elif opt in ('-v','--freqmax'):
        vmax = float(arg)
    elif opt in ('-m','--nstatesmin'):
        nstatesmin = int(arg)
    elif opt in ('-x','--nstatesmax'):
        nstatesmax = int(arg)
    elif opt in ('-s','--nstatesstep'):
        nstatesstep = int(arg)

out = Molecule(infilename)

if nstatesmin and nstatesmax and nstatesstep:
    avg = True
    states = range(nstatesmin, nstatesmax, nstatesstep)
elif nstatesmin or nstatesmax or nstatesstep:
    print("Error: To average over a range of states, the min, max, and step size must all be specified")
    print("Continuing by using one value for number of states, ", nstates)
    states = [nstates]
else:
    states = [nstates]

# Get basic info
afac = width/(2.0*math.pi)
bfac = width/2.0

crs_avg = []
out.read_atoms()
vibs = VibAll()
for v in vibs.modes:
    if v.freq > vmin and v.freq < vmax:
        crs_list = []
        ls = ''
        for s in states:
            if s == states[0]: v.norm_mode(out.atoms)
            v.alpha_slope(outfilename,disp_str,omega,gamma,s)
            v.raman_scat(coord)
            v.raman_cross(omega)
            crs_list.append(v.crs_tot[0])
            ls += string.rjust(str(v.crs_tot[0]*scale*afac/bfac**2),17) + ' '

        # Find avg and stdev
        if avg and len(states) > 1:
            mean = 0.0
            stdev = 0.0
            for i in crs_list:
                mean += i
            mean /= len(crs_list)
            for i in crs_list:
                stdev += (i - mean)**2
            stdev /= len(crs_list)-1
            stdev = math.sqrt(stdev)
            #print ls
            if mean > 0.0 and stdev/mean > 0.33:
                print(v.index, v.freq, mean*scale*afac/bfac**2, stdev*scale*afac/bfac**2, stdev/mean)
                print(ls)
            else:
                print(v.index, v.freq, mean*scale*afac/bfac**2, stdev*scale*afac/bfac**2)
            crs_avg.append([v.freq,mean,stdev]) 
        else:
            print(v.index, v.freq, crs_list[0]*scale*afac/bfac**2)
            crs_avg.append([v.freq]) 

# Write stick spectrum
stick = open(outfilename+'.raman_stick_'+str(omega)+'_avg','w')
for k in crs_avg:
    stick.write(string.rjust('%.3f'%k[0],10))
    if len(k) > 1:
        stick.write(string.rjust(str(k[1]*scale*afac/bfac**2),20) + string.rjust(str(k[2]*scale*afac/bfac**2),20)) 
    stick.write('\n')
    #print k[1]*scale*afac/bfac**2, k[2]*scale*afac/bfac**2
stick.close()

# Set up array of zeros for spectrum
freq_min = 1.0
freq_max = 4000.0
freq_step = 1.0
spectrum = []
for i in range(0,int((freq_max - freq_min)/freq_step + 1)):
    spectrum.append([0.0,0.0,0.0])

'''
# Compute lorentzian-broadened spectrum
for v in vibs.modes:
    if v.freq > vmin and v.freq < vmax:
        for j in range(0,len(spectrum)):
            f_curr = freq_min + freq_step*j
            factor = scale * afac/( (f_curr - v.freq)**2 + bfac**2)
            spectrum[j][0] += factor * v.crs_tot[0]
            #spectrum[j][1] += factor * v.crs_real
            #spectrum[j][2] += factor * v.crs_imag

spec = open(outfilename+'.raman_lrntz_'+str(omega)+'_'+str(states),'w')
for j in range(0,len(spectrum)):
    #print freq_min + freq_step*j, spectrum[j][0]
    spec.write(string.rjust('%.1f'%(freq_min + freq_step*j),8) + string.rjust(str(spectrum[j][0]),20) + '\n')
spec.close()
'''
