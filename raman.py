#-------------------------------------------------------------------------
# Compute Raman cross-sections from changes to the INDO/CI  
# polarizabilities with geometric displacement along normal modes 
# from another level of theory
#
# By default, the Raman intensities are computed using all excited
# states in the INDO/CI output.
#
# Alternatively, the Raman intensities may be computed using a series of
# numbers of excited states, then averaging. This can be useful when the 
# system is so large that the SOS expression is not fully converged.
# 
#-------------------------------------------------------------------------

import string,sys,math
from utils.molecule import Molecule
from utils.vibrations import VibAll
import getopt

# Default values for inputs
infilename = ''
types = []
omega = 0.0
gamma = 0.1088j
disp_str = '0.01'
disp = float(disp_str)
coord = ''
nstates = 0
nstatesmin, nstatesmax, nstatesstep = None, None, None
avg = False
vmin = 350.0
vmax = 2500.0
vstep = 1.0
vwidth = 20.0
vibfile = 'nmodes.inp'

helpfile = """
raman.py -i <inputfile> -o <outputfile> -e <energy> -g <gamma> -d <displacement> -c <coordinate> -n <nrsfile> -f <minfreq> -v <maxfreq> -n <nstates> -m <min_nstates> -x <max_nstates> -s <step_nstates>

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
    -o    Base file name for output file                            Default = matches -i
    -f    File containing normal modes                              Default = nmodes.inp
    -e    Energy at which Raman intensities are computed (eV)       Default = 0.0
    -g    Lifetime (broadening) for polarizabilities (eV)           Default = 0.1088 (= 0.004 a.u.)
    -d    Displacement of geometries along vibrational modes (A)    Default = 0.01
    -c    Coordinate along which to compute Raman intensities       Default = isotropic (alternatives are x, y, z)
    -n    Number of states to include in SOS expression             Default = All states
    -v    Minimum vibrational frequency of modes to use (cm-1)      Default = 350
    -w    Maximum vibrational frequency of modes to use (cm-1)      Default = 2500
    -s    Step size for Lorentzian-broadened Raman spectrum (cm-1)  Default = 1
    -b    Line width for Lorentzian-broadened Raman spectrum (cm-1) Default = 20

To compute the Raman intensities averaged over a series of numbers of excited states
instead of one number, use the following options. If one is used, all must be used.
This option is sometimes useful for large systems where the Raman intensities are 
not fully converged with respect to the number of states in the SOS expression.
    -m    Minimum number of excited states to include in SOS expression
    -x    Maximum number of excited states to include in SOS expression
    -z    Step size in number of excited states to include in SOS expression
"""

# Parse input options
try:
    options, remainder = getopt.getopt(sys.argv[1:],"hi:o:f:e:g:d:c:v:w:s:b:m:x:z:",[
        '--help','--input=','--output=','--modefile=','--energy=','--gamma=',
        '--displacement=','--coord=','--nstates=',
        '--freqmin=', '--freqmax=', '--freqstep=','--broadening=',
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
    elif opt in ('-f','--modefile'):
        vibfile = arg
    elif opt in ('-g','--gamma'):
        gamma = complex(0.0,float(arg))
    elif opt in ('-d','--disp'):
        disp_str = arg
        disp = float(arg)
    elif opt in ('-c','--coord'):
        coord = arg
    elif opt in ('-n','--nstates'):
        nstates = int(arg)
    elif opt in ('-v','--freqmin'):
        vmin = float(arg)
    elif opt in ('-w','--freqmax'):
        vmax = float(arg)
    elif opt in ('-s','--freqstep'):
        vstep = float(arg)
    elif opt in ('-b','--broadening'):
        vwidth = float(arg)
    elif opt in ('-m','--nstatesmin'):
        nstatesmin = int(arg)
    elif opt in ('-x','--nstatesmax'):
        nstatesmax = int(arg)
    elif opt in ('-z','--nstatesstep'):
        nstatesstep = int(arg)

out = Molecule(infilename)

if nstatesmin and nstatesmax and nstatesstep:
    avg = True
    states = range(nstatesmin, nstatesmax, nstatesstep)
elif nstatesmin or nstatesmax or nstatesstep:
    if nstatesmin:
        nstates = nstatesmin
    elif nstatesmax: 
        nstates = nstatesmax
    print("Error: To average over a range of states, the min, max, and step size must all be specified")
    print("Continuing by using one value for number of states, ", nstates)
    states = [nstates]
else:
    states = [nstates]

# Constants for Lorentzian-broadened spectrum
scale = 1E+34
afac = vwidth/(2.0*math.pi)
bfac = vwidth/2.0

# Get basic info
crs_avg = []
out.read_atoms()
vibs = VibAll(vibfile=vibfile, vmin=vmin, vmax=vmax)
if avg and len(states) > 1:
    print('                 ----------Raman intensity----------')
    print('Mode  Frequency  Mean        Stdev       Stdev/Mean')
else:
    print('Mode  Frequency  Raman intensity')
for v in vibs.modes:
    if v.freq > vmin and v.freq < vmax:
        crs_list = []
        ls = ''
        for s in states:
            if s == states[0]: 
                v.norm_mode(out.atoms)
            v.alpha_slope(infilename,outfilename,disp_str,omega,gamma,s)
            v.raman_scat(coord)
            v.raman_cross(omega)
            crs_list.append(v.crs_tot[0])
            ls += str(v.crs_tot[0]*scale*afac/bfac**2).rjust(17) + ' '

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
            if mean > 0.0 and stdev/mean > 0.5:
                print(str(v.index).rjust(4) + ('%.3f' % v.freq).rjust(11) + 
                      ('%.4e' % (mean*scale*afac/bfac**2)).rjust(12) + 
                      ('%.4e' % (stdev*scale*afac/bfac**2)).rjust(12) +
                      ('%.3f' % (stdev/mean)).rjust(12))
                print('  Caution: Mode ' + str(v.index) + 
                      ' has large variation in Raman intensity for different numbers of states in SOS')
                print('  States   Raman intensity')
                for j in range(0,len(states)):
                    print(str(states[j]).rjust(8) + ('%.4e' % (crs_list[j]*scale*afac/bfac**2)).rjust(12))
            else:
                print(str(v.index).rjust(4) + ('%.3f' % v.freq).rjust(11) + 
                      ('%.4e' % (mean*scale*afac/bfac**2)).rjust(12) + 
                      ('%.4e' % (stdev*scale*afac/bfac**2)).rjust(12))
            crs_avg.append([v.index,v.freq,mean,stdev]) 
        else:
            print(str(v.index).rjust(4) + ('%.2f' % v.freq).rjust(11) 
                  + ('%.4e' % (crs_list[0]*scale*afac/bfac**2)).rjust(12))
            crs_avg.append([v.index,v.freq, crs_list[0]]) 

# Write stick spectrum
stick = open(outfilename+'.raman_stick_'+str(omega),'w')
if len(crs_avg[0]) > 3:
    stick.write('                 ----------Raman intensity----------\n')
    stick.write('Mode  Frequency  Mean        Stdev       Stdev/Mean\n')
else:
    stick.write('Mode  Frequency  Raman intensity\n')
for k in crs_avg:
    stick.write(str(k[0]).rjust(4) + ('%.3f'%k[1]).rjust(11))
    if len(k) > 3:
        stick.write(('%.4e'%(k[2]*scale*afac/bfac**2)).rjust(12) 
                    + ('%.4e'%(k[3]*scale*afac/bfac**2)).rjust(12)) 
    else: 
        stick.write(('%.4e'%(k[2]*scale*afac/bfac**2)).rjust(12))
    stick.write('\n')
    #print k[1]*scale*afac/bfac**2, k[2]*scale*afac/bfac**2
stick.close()

# Set up array of zeros for spectrum
spectrum = []
for i in range(0,int((vmax - vmin)/vstep + 1)):
    spectrum.append([0.0,0.0,0.0])

# Compute lorentzian-broadened spectrum
for v in vibs.modes:
    if v.freq > vmin and v.freq < vmax:
        for j in range(0,len(spectrum)):
            f_curr = vmin + vstep*j
            factor = scale * afac/( (f_curr - v.freq)**2 + bfac**2)
            spectrum[j][0] += factor * v.crs_tot[0]
            #spectrum[j][1] += factor * v.crs_real
            #spectrum[j][2] += factor * v.crs_imag

spec = open(outfilename+'.raman_lrntz_'+str(omega),'w')
spec.write('Frequency  Raman intensity\n')
for j in range(0,len(spectrum)):
    #print freq_min + freq_step*j, spectrum[j][0]
    spec.write(('%.3f'%(vmin + vstep*j)).rjust(9) 
               + ('%.4e'%spectrum[j][0]).rjust(12) + '\n')
spec.close()

