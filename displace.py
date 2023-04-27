#-------------------------------------------------------------------------
# Displace a geometry along a normal mode by a specified distance. 
# Requires:
#  1. A MOPAC input file for the equilibrium geometry. The keywords from
#     this file will be used for the displaced geometry.
#  2. A file containing the normal modes (by default, nmodes.inp). This
#     file must have the following format:
'''
              623.3740               623.3740              1306.9964
         ------------------    -------------------    -------------------
   O     0.00  -0.33   0.00    -0.33  -0.00  -0.00    -0.00   0.00   0.71
   C    -0.00   0.88  -0.00     0.88   0.00  -0.00     0.00  -0.00  -0.00
   O     0.00  -0.33  -0.00    -0.33  -0.00   0.00    -0.00   0.00  -0.71


             2359.8986
        -------------------
   O    -0.00  -0.00  -0.33
   C     0.00   0.00   0.88
   O    -0.00  -0.00  -0.33
'''
#    In the nmodes file, the number of lines matters, but the exact 
#    spacing on each line does not.
#
#-------------------------------------------------------------------------

import string,sys
import utils.constants as cnst
import getopt

# Default values
infilename = ''
disp = 0.01
mode = 1
modefilename = 'nmodes.inp'

helpfile = """
displace.py -i <inputfile> -o <outputfile> -m <mode> -d <displacement> -f <modesfile>

Required:
    -i    Base file name for input file

Optional:
    -o    Base file name for output file                          Default = matches -i
    -m    Integer number of vibrational mode                      Default = 1
    -d    Displacement of geometries along vibrational modes (A)  Default = 0.01
    -f    File containing normal modes                            Default = nmodes.inp

The output will be a MOPAC input file with name inputfile_mode_displacement.mop
"""

try:
    options, remainder = getopt.getopt(sys.argv[1:],"hi:o:m:d:f:",['--help','--input=','--output=','--mode=','--displacement=','--modefile='])
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
    elif opt in ('-m','--mode'):
        mode = int(arg)
    elif opt in ('-d','--displacement'):
        disp = float(arg)
    elif opt in ('-f','--modefile'):
        modefilename = arg


# Open files and prepare input
modes = open(modefilename,'r')
inp = open(infilename+'.mop','r')
new = open(outfilename+'_'+str(mode)+'_'+str(disp)+'.mop','w')

# Read the initial geometry
for i in range(0,3):
    iline = inp.readline()
    new.write(iline)

iline = inp.readline()
xyzcoord = []
while len(iline) > 2:
    line = iline.split()
    for i in [1,3,5]:
        line[i] = float(line[i])
    xyzcoord.append([line[0],line[1],line[3],line[5]])
    iline = inp.readline()

# Find and read the correct mode displacements
mline = modes.readline()
line = mline.split()
currmode = len(line)
while currmode < mode and len(mline) > 0:
    while len(mline) > 5:
        mline = modes.readline()
        #print mline, len(mline)
    mline = modes.readline()
    mline = modes.readline()
    line = mline.split()
    currmode += len(line)

currstart = currmode - len(line) + 1
mode_ind = mode - currstart

j = mode_ind*3 + 1
mline = modes.readline()
mline = modes.readline()
mode_disp = []
for i in range(0,len(xyzcoord)):
    line = mline.split()
    mode_disp.append([float(line[j]),float(line[j+1]),float(line[j+2])])
    mline = modes.readline()

# Print new input file
for i in range(0,len(xyzcoord)):
    new.write(str(xyzcoord[i][0]).ljust(4))
    for j in range(0,3):
        new.write(('%6f'%(xyzcoord[i][j+1]+disp*mode_disp[i][j])).rjust(11) + ' 0')
    new.write('\n')
    iline = inp.readline()
new.write('\n')

while len(iline) > 0:
    iline = inp.readline()
    new.write(iline)


modes.close()
inp.close()
new.close()

