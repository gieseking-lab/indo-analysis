#!/usr/bin/python

#-------------------------------------------------------------------------
# Compute Raman cross-sections from changes to the INDO/CI static 
# polarizabilities with geometric displacement along DFT normal modes
#
# Replicate the math from Lasse Jensen's script ramanscat.f90
#
# Version 2.0
#
# Rebecca Gieseking, 1/6/2017
# 
#-------------------------------------------------------------------------

import string,sys,math
import constants as cnst
from mol_assign import *
from vibrations import *
import getopt

# Constants for later
types = []
omega = 0.0
usefreq = False
gamma = 0.1088j
disp_str = '0.01'
disp = float(disp_str)
coord = ''
states = 0
width = 20.0
scale = 1E+34
prog = ' '
gamma2 = 0.0
tfilename = 'types'
nrsfile = ' '

vmin = 350
vmax = 2500

helpfile = """
indo_raman.py -i <inputfile> -o <outputfile> -w <omega> -f <usefreq> -g <gamma> -b <gamma2> -t <types> -d <displacement> -c <coordinate> -s <# states> -n <nrsfile> -z <program>

    -i    Base file name for input and output (must be included)
    -o    Base file name for all output files, if different from input

    -w    Energy at which polarizabilities are computed (eV)        Default = 0.0
    -f    Include vibrational mode frequency in alpha energy        Default = False
    -g    Lifetime (broadening) for polarizabilities (eV)           Default = 0.1088j (= 0.004 a.u.)
    -d    Displacement of geometries along vibrational modes (A)    Default = 0.01
    -c    Coordinate along which to compute Raman intensities       Default = isotropic (alternatives x, y, z)
    -b    Lifetime of excited states of specific types (eV)         Default = 0.0 (= not used)
    -t    Types file, for used with -b. Only first type used.       Default = types
          Format:
             Atom    20               From atoms 1-20 to 21-max
             Attype  Ag               From Ag atoms to all other atoms
             Orbtype Ag D             From Ag D to all others
             Orb     D                From all D orbitals to all others

    -s    Number of states to include in SOS expression             Default = All states
    -n    If S factors already computed, file containing those      Default = none

    -z    Program generating output (Mopac, ADF, Qchem, etc.)       Default = Determine automatically
"""

# Parse input options
try:
    options, remainder = getopt.getopt(sys.argv[1:],"hi:o:w:f:g:b:t:d:c:s:n:z:",['--help','--input=','--output=','--omega=','--usefreq=','--gamma=','--gamma2=','--types=','--disp=','--coord=','--states=','--nrsfile=','--program='])
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
        if '.' in arg:
            arg = arg[:arg.rfind('.')]
        infilename = arg
        outfilename = arg
    elif opt in ('-o','--output'):
        if '.' in arg:
            arg = arg[:arg.rfind('.')]
        outfilename = arg
    elif opt in ('-w','--omega'):
        omega = float(arg)
    elif opt is ('-f','--usefreq='):
        if arg.lower() in ('true','t','yes','y','1'):
            usefreq = True
    elif opt in ('-g','--gamma'):
        gamma = complex(0.0,float(arg))
    elif opt in ('-b','--gamma2'):
        gamma2 = complex(0.0,float(arg))
    elif opt in ('-t','--types'):
        tfilename = arg
    elif opt in ('-d','--disp'):
        disp_str = arg
        disp = float(arg)
    elif opt in ('-c','--coord'):
        coord = arg
    elif opt in ('-s','--states'):
        states = int(arg)
    elif opt in ('-n','--nrsfile'):
        nrsfile = arg
    elif opt in ('-z','--program'):
        prog = arg.lower()

types = read_types(tfilename)

out, prog = open_mol(infilename, prog)
afac = width/(2.0*math.pi)
bfac = width/2.0

# Get basic info
out.read_atoms()
vibs = VibAll(prog=prog,nrsfile=nrsfile)
for i in range(0,len(vibs.modes)):
    v = vibs.modes[i]
    if v.freq > vmin and v.freq < vmax:
        #print i+1, v.freq
        if nrsfile == ' ':
            #print 'Mode ',v.freq
            v.norm_mode(out.atoms)
            #print v.index, v.freq, v.norm
            v.alpha_slope(outfilename,disp_str,omega,gamma,states,gamma2,types,usefreq,v.freq)
            #adiff_or = (alpha_diff[0][0] + alpha_diff[1][1] + alpha_diff[2][2]) / 3.0
            #print adiff_or.real, adiff_or.imag, abs(adiff_or)
            v.raman_scat(coord)
        v.raman_cross(omega)
        print v.crs_tot[0]*scale*afac/bfac**2

# Write stick spectrum
stick = open(outfilename+'.raman_stick_'+str(omega)+'_'+str(states),'w')
for v in vibs.modes:
    if v.freq > vmin and v.freq < vmax:
        stick.write(string.rjust('%.3f'%v.freq,10))
        for k in v.crs_tot:
            stick.write(string.rjust(str(k*scale*afac/bfac**2),20)) 
        stick.write('\n')
        #print v.crs_tot[0]*scale*afac/bfac**2, v.crs_real[0]*scale*afac/bfac**2, v.crs_imag[0]*scale*afac/bfac**2
        #print v.crs_tot[0]*scale*afac/bfac**2
stick.close()

# Set up array of zeros for spectrum
freq_min = 1.0
freq_max = 4000.0
freq_step = 1.0
spectrum = []
for i in range(0,int((freq_max - freq_min)/freq_step + 1)):
    spectrum.append([0.0,0.0,0.0])

afac = width/(2.0*math.pi)
bfac = width/2.0

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

