# Calculate TPA cross-sections of exicted states using two-state and three-state models
# Input dipole moments in Debye, energies in eV
# Energies are converted from eV to 10^19 J (same as gamma3state.py)
# Omega output in eV, Gamma output in esu x 10^-36
# Math from Sukrit's TPA program

import string, sys, math, numpy, itertools

# Conversion factor to joules
joules = 1.6022

# INPUT SECTION
# Get the file names
if len(sys.argv)==1:
    file = raw_input("Enter the ground state file:")
else:
    file = sys.argv[1]

# Clean up the file name
if "." in file:
   file = file[:file.find(".")]

# Open input and output files
# 
# Edit - change output file to .tpa
inp = open(file+".inp",'r')
out2 = open(file+".tpa2",'w')

# Read initial data from input file
line = inp.readline()
numstates = int(line[18:21])

line = inp.readline()
omegamin = float(line[11:17])
omegamax = float(line[26:31])
omegastep = float(line[41:46])

line = inp.readline()
line = inp.readline()
onephoton = []
while "State" not in line:
 onephoton.append(int(line[:5]))
 line = inp.readline()

# Set up arrays for data to be read from input file
muge = numpy.zeros([numstates,3])
muee = numpy.zeros([numstates,numstates,3])
dmuee = numpy.zeros([numstates,3])
state = numpy.zeros([numstates])

# Read states and dipoles
# Convert state energies from eV to J x 10^19
#
# Edited - set up omega to equal state energies
omega = []
line = inp.readline()
for i in range(0,numstates):
 omega.append(float(line[:20]))
 state[i] = float(line[:20])*joules
 line = inp.readline()
line = inp.readline()

for i in range(0,numstates):
 muge[i][0] = float(line[:20])
 muge[i][1] = float(line[21:40])
 muge[i][2] = float(line[41:60])
 line = inp.readline()
line = inp.readline()
line = inp.readline()

for j in range(0,len(onephoton)):
 for i in range(0,numstates):
  muee[onephoton[j]][i][0] = float(line[:20])
  muee[onephoton[j]][i][1] = float(line[21:40])
  muee[onephoton[j]][i][2] = float(line[41:60])
  line = inp.readline()
 line = inp.readline()

for i in range(0,numstates):
 dmuee[i][0] = float(line[:20])
 dmuee[i][1] = float(line[21:40])
 dmuee[i][2] = float(line[41:60])
 line = inp.readline()


tpa = numpy.zeros([3,3,len(omega)])
loopnum = 0
tpaavg = numpy.zeros([2,len(omega)])

# END INPUT/SETUP
# START CALCULATION OF TPA CROSS-SECTION

out2.write("Energy\t\tDipolar TPA\t\tThree-state TPA\n")
#print "Energy\t\tDipolar TPA\t\tThree-state TPA"

for ws_eV in omega:
 ws = ws_eV*joules

 # Dipolar TPA
 # Set up TPA in each set of dimensions
 for i in range(0,3):
  for j in range(0,3):
   tpa[i][j][loopnum] = (muge[loopnum][i]*dmuee[loopnum][j] + muge[loopnum][j]*dmuee[loopnum][i])/(ws_eV - ws_eV/2)

 # Average over all dimensions
 for i in range(0,3):
  for j in range(0,3):
   tpaavg[0,loopnum] = tpaavg[0,loopnum] + tpa[i][i][loopnum]*tpa[j][j][loopnum] + 2*tpa[i][j][loopnum]*tpa[j][i][loopnum]

 tpaavg[0,loopnum] = tpaavg[0,loopnum]*(ws_eV/2)**2*(6.4235/600)/(1.5*joules**3)
 
 # Three-state TPA
 # Set up TPA in each set of dimensions
 for i in range(0,3):
  for j in range(0,3):
   tpa[i][j][loopnum] = (muge[onephoton[0]][i]*muee[onephoton[0]][loopnum][j] + muge[onephoton[0]][j]*muee[onephoton[0]][loopnum][i])/(omega[onephoton[0]] - ws_eV/2)

 # Average over all dimensions
 for i in range(0,3):
  for j in range(0,3):
   tpaavg[1,loopnum] = tpaavg[1,loopnum] + tpa[i][i][loopnum]*tpa[j][j][loopnum] + 2*tpa[i][j][loopnum]*tpa[j][i][loopnum]

 tpaavg[1,loopnum] = tpaavg[1,loopnum]*(ws_eV/2)**2*(6.4235/600)/(1.5*joules**3)
 if loopnum != 0:
  out2.write(str(omega[loopnum]) + "  \t" + str(tpaavg[0][loopnum]) + "    \t" + str(tpaavg[1][loopnum]) + "\n")
#  print str(omega[loopnum]) + "  \t" + str(tpaavg[0][loopnum]) + "    \t" + str(tpaavg[1][loopnum])
 loopnum = loopnum+1

inp.close()
out2.close()

