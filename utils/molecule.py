################################################################
# Molecule class                                               #
#                                                              #
# Contains functions to read from MOPAC input files and        #
# write various types of output information                    #
#                                                              #
# Information that can be read:                                #
#    method (INDO, PM7, etc.)                                  #
#    charge of system                                          #
#    atoms (elements and initial coordinates)                  #
#    norbs (# molecular orbitals)                              #
#    dipole                                                    #
#    orbs (MO coefficient matrix)                              #
#    configs (INDO only - CI configurations)                   #
#    states (INDO only - excited state energies and configs)   #
#                                                              #
# Information that can be written:                             #
#    orbs  = MO numbers and energies                           #
#    osc   = Excited state energies and oscillator strengths   #
#    exc   = Excited-state analysis for plasmon ID             #
#    sigma = Absorption spectrum                               #
#    alpha = Polarizability                                    #
#                                                              #
################################################################

import math
import utils.constants as cnst
from utils.properties import Atom, MolOrb, Config, State

# Error function
def erf(x):
    # save the sign of x
    sign = 1 if x >= 0 else -1
    x = abs(x)

    # constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    # A&S formula 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)
    return sign*y # erf(-x) = -erf(x)

# Read types file
def read_types(tfilename):
    types = []
    try:
        tfile = open(tfilename,'r')
        for line in tfile.readlines():
            print('Type   %s' % line)
            types.append(line.split())
    except IOError:
        print('%s not found; no types assigned' % tfilename)

    return types


class Molecule(object):
    def __init__(self, filename):
        try:
            file = open(filename,'r')
        except IOError:
            try:
                file = open(filename+'.out','r')
            except IOError:
                file = open(filename+'.log','r')
        self.file = file
        self.line = self.file.readline()

        # Identify what has been read
        self.method_read   = False
        self.charge_read   = False
        self.atoms_read    = False
        self.norbs_read    = False
        self.dipole_read   = False
        self.orbs_read     = False
        self.configs_read  = False
        self.states_read   = False

    ################################################
    #                                              #
    #  General read functions                      #
    #                                              #
    ################################################

    # Read to a certain message in the file
    # Use a gibberish message2 to make it useless unless specified
    def read_to(self, message, message2='qwertyuiopasdf'):
        #print('Reading to', message, message2)
        while message not in self.line and message2 not in self.line and len(self.line) > 0:
            self.line = self.file.readline()
        #print('Found message', self.line)

    # Read a specific number of lines in the file
    def read_for(self,lines):
        for i in range(0,lines):
            self.line = self.file.readline()

    ################################################
    #                                              #
    #  Read MOPAC input                            #
    #                                              #
    ################################################
    def read_method(self):
        if self.charge_read or self.atoms_read or self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)
        self.read_to('CALCULATION DONE')
        self.read_for(1)
        line = self.line.split()
        self.method = line[1]
        self.method_read = True
    
    def read_charge(self):
        if self.method_read or self.atoms_read or self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)
        
        self.read_to('CALCULATION DONE')
        self.read_to('CHARGE','**************************************')
        self.charge = int(self.line.split()[-1]) if 'CHARGE' in self.line else 0

        # Save that atoms have been read
        self.charge_read = True
        #print("Charge read")

    # Read the atomic coordinates and set up atoms
    def read_atoms(self):
        if self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

        self.read_to('CARTESIAN COORDINATES')
        self.read_for(4)
        self.atoms = []
        self.at_types = []
        
        while len(self.line) > 1:
            line = self.line.split()
            self.atoms.append(Atom([float(line[2]),float(line[3]),float(line[4])],line[1]))
            if line[1] not in self.at_types:
                self.at_types.append(line[1])
            self.line = self.file.readline()

        # Save that atoms have been read
        if len(self.atoms) > 0:
            self.atoms_read = True
            #print("Atoms read")
        else:
            print("Error: Atoms not read")

    # Read the number of molecular orbitals
    def read_norbs(self):
        if self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

        self.read_to('SHELL')
        self.read_for(1)
        line = self.line.split()
        self.homo = 0
        for i in range(1,len(line)):
            self.homo += int(line[i])


        self.read_to('TOTAL')
        line = self.line.split()
        self.numorb = 0
        for i in range(1,len(line)):
            self.numorb += int(line[i])

        # Save that norbs has been read
        if self.numorb > 0:
            self.norbs_read = True
            #print("Orbital count read")
        else:
            print("Error: Number of orbitals not read")

    # Read the ground state dipole moment into gdip
    def read_dipole(self):
        if self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

        self.read_to('Dipole moment=')
        self.gdip = float(self.line.split()[-2])
        #print(self.gdip)

        # Save that dipole has been read
        if 'Dipole moment=' in self.line:
            self.dipole_read = True
            #print("Dipole read")
        else:
            print("Error: Ground-state dipole not read")

    # Read the molecular orbitals
    def read_orbs(self, types=[], mo_min=1, mo_max=0):
        # Make sure atoms and norbs have been read
        if not self.atoms_read:
            self.read_atoms()
        if not self.norbs_read:
            self.read_norbs()
        if self.configs_read or self.states_read:
            self.file.seek(0)


        if mo_max == 0:
            mo_max = self.numorb
        self.mos = []

        #self.read_to('eigvals(Ev)')
        self.read_to('ROOT NO.')
        orb_count = 0
        aos_read = False
        orb_start = False
        orb_done = False

        # Enter main loop to get molecular orbitals
        while orb_done == False and len(self.line) > 0:
            self.read_for(2)
            line = self.line.split()
            orbline = len(line) 
            orb_count += orbline
            #print(line, orb_count, orbline)

            if mo_min <= orb_count:
                # Read the relevant orbital info

                # Set up indices for loops later on
                if orb_start == False:
                    orb_index_a = mo_min - orb_count - 1
                else:
                    orb_index_a = -orbline
                if mo_max <= orb_count:
                    orb_index_b = mo_max - orb_count
                else:
                    orb_index_b = 0
       
                #print(orb_index_a, orb_index_b)
 
                # Add each MO to the list with an energy
                for i in range(0,len(line)):
                   x = MolOrb(float(line[i]))
                   self.mos.append(x)

                self.read_for(3) 
                        
                # Loop over each AO in this set of MOs
                for i in range(0,self.numorb):
                    line = self.line.split()
                    if len(line) == 0:
                        self.line = self.file.readline()
                        line = self.line.split()
    
                    # First orbital set only - store atomic orbital types
                    if aos_read == False:
                        at = int(line[2]) - 1
                        aotype = line[0]
                        if 'S' not in aotype and 'P' not in aotype:
                            aotype = 'D' + aotype
                        #if 'xy2' in aotype: aotype = 'Dxy2'
                        self.atoms[at].add_aotype(aotype)
                        
                        # Categorize all AOs by types from type file
                        a = []
                        for j in types:
                            #for listelem in j:
                            #   print(listelem, aotype, listelem in aotype)
                            if 'Atom' in j:
                                a.append(True) if (at >= (int(j[1])-1) and at < int(j[2])) else a.append(False)
                            elif 'Attype' in j:
                                a.append(True) if self.atoms[at].elem in j else a.append(False)
                            elif 'Orb' in j:
                                a.append(True) if any(listelem in aotype for listelem in j) else a.append(False)
                            elif 'Orbtype' in j:
                                a.append(True) if (any(listelem in aotype for listelem in j) and self.atoms[at].elem in j) else a.append(False)
                            elif 'Orbat' in j: 
                                a.append(True) if (at >= (int(j[1])-1) and at < int(j[2]) and any(listelem in aotype for listelem in j)) else a.append(False)
                            else:
                                 a.append(True)

                        self.atoms[at].add_ao(a)
                        #self.atoms[at].add_aotype(aotype) 
                        #print(at, aotype, a)

                        #if 'S' or 'P' or 'D' in type:
                        #    self.atoms[at].add_ao(type)
                        #else:
                        #    self.atoms[at].add_ao('Dxy2')

                    # Store the MO coefficients
                    # WARNING - ordering of i/j is reversed from previous scripts!!! ############
                    for j in range(orb_index_a,orb_index_b):
                        self.mos[j].add_coeff(line[j])

                    self.read_for(1)

                aos_read = True
                self.read_to('ROOT NO.','Reference')
                orb_start = True
                if mo_max <= orb_count:
                    orb_done = True

            else:
                # Skip to the next section
                self.read_to('ROOT NO.','Reference')

        # Compute MO character based on AO types
        for i in self.mos:
            for m in range(0,len(types)):
                i.char.append(0.0)
                ocycle = 0
                for j in self.atoms:
                    for k in j.ao:
                        if k[m]: i.char[m] += i.coeff[ocycle]**2
                        ocycle += 1
            #print(i.char)

        # Save that orbitals have been read
        if orb_count > 0:
            self.orbs_read = True
            #print("Orbitals read", orb_count)
        else:
            print("Error: Orbitals not read")


    # Read electron configurations for CI
    def read_configs(self,types=[]):
        if self.states_read:
            self.file.seek(0)
        if not self.orbs_read:
            self.read_orbs(types)

        self.read_to('spin-adapted configurations')
        self.read_for(4)

        self.configs = []
        line = self.line.split()
        while len(self.line) > 5:
            ener = float(line[1])
            try:
                occ = int(self.line[53:57]) - 1
                vir = int(self.line[62:68]) - 1
            except ValueError:
                occ = 0
                vir = 0
            char = []
            achar = []
            bchar = []

            ochar = []
            vchar = []

            for i in range(0,len(self.mos[0].char)):
                char.append(self.mos[vir].char[i] - self.mos[occ].char[i])

                # achar = local character in orbitals included in type
                # bchar = local character in orbitals  outside of type
                achar.append(    min(self.mos[vir].char[i],self.mos[occ].char[i]))
                bchar.append(1 - max(self.mos[vir].char[i],self.mos[occ].char[i]))

                ochar.append(self.mos[occ].char[i])
                vchar.append(self.mos[vir].char[i])

                # RLG temporary edit for ligand clusters
                #char.append(self.mos[occ].char[i])
            self.configs.append(Config(occ,vir,char))
            self.configs[-1].achar = achar
            self.configs[-1].bchar = bchar
            self.configs[-1].ochar = ochar
            self.configs[-1].vchar = vchar
            self.configs[-1].change_energy(ener)

            #if len(self.configs) < 20:
            #    print(char[0], achar[0], bchar[0], abs(char[0])+achar[0]+bchar[0])

            self.line = self.file.readline()
            line = self.line.split()

        # Save that configs have been read
        if len(self.configs) > 0:
            self.configs_read = True
            #print("CI configurations read")
        else:
            print("Error: CI configurations not read")


    # Read excited-state info
    def read_states(self,types=[]):
        # Make sure configs have been read
        if not self.configs_read:
            self.read_configs(types)

        self.read_to('CI trans.')
        self.read_for(3)

        self.states = [State(0.0,0.0)]

        # Read energies and oscillator strengths
        line = self.line.split()
        while len(line) > 2:
            #num = int(line[0])
            try:
                e = float(line[1])
                o = float(line[4])
            except ValueError:
                e = float(line[2])
                o = float(line[5])
            self.states.append(State(e,o))
            if o > 1.e-5:
                pol = [float(self.line[50:60]),float(self.line[61:71]),float(self.line[72:82])]
            else:
                pol = [0.0,0.0,0.0]
            self.states[-1].calc_tdip(pol)
        
            self.line = self.file.readline()
            line = self.line.split()

        # Read CI coefficients
        if len(types) > 0 or self.configs_read == True:
            self.read_to('State')
            nstate = 0

            while 'State' in self.line:
                self.line = self.file.readline()

                ntypes = len(self.mos[0].char)
                char = []
                achar = []
                bchar = []

                # ochar & vchar for identifying character of occ/vir orbitals involved
                # Added for parametrization & auto-identification of states
                ochar = []
                vchar = []

                for i in range(0,ntypes): 
                    char.append(0.0)
                    achar.append(0.0)
                    bchar.append(0.0)

                    ochar.append(0.0)
                    vchar.append(0.0)

                while 'Config' in self.line:
                    line = self.line.split()
                    cnum = int(line[1])
                    coeff = float(line[3])
                    self.states[nstate].add_coeff(cnum, coeff)

                    if 'Config       1' in self.line:
                        while 'Config' in self.line:
                            self.line = self.file.readline()
                    else:
                        for i in range(0,ntypes):
                            char[i]  += coeff * self.configs[cnum-1].char[i]
                            achar[i] += coeff * self.configs[cnum-1].achar[i]
                            bchar[i] += coeff * self.configs[cnum-1].bchar[i]

                            ochar[i] += coeff * self.configs[cnum-1].ochar[i]
                            vchar[i] += coeff * self.configs[cnum-1].vchar[i]


                        self.line = self.file.readline()
                line = self.line.split()
                totcoeff = float(line[3])
                if totcoeff > 0:
                    for i in range(0,len(char)):
                        char[i] /= totcoeff
                        achar[i] /= totcoeff
                        bchar[i] /= totcoeff
                        ochar[i] /= totcoeff
                        vchar[i] /= totcoeff


                self.states[nstate].char = char
                self.states[nstate].achar = achar
                self.states[nstate].bchar = bchar
                self.states[nstate].ochar = ochar
                self.states[nstate].vchar = vchar

                #print(nstate, char, ochar, vchar)

                self.line = self.file.readline()
                self.line = self.file.readline()

                nstate += 1

            #for i in range(0,11):
            #    print(self.states[i].energy, self.states[i].tdip)

        # Save that states have been read
        if len(self.states) > 0:
            self.states_read = True
            #print("Excited states read")
        else:
            print("Error: Excited states not read")

    ################################################
    #                                              #
    #  Compute properties                          #
    #                                              #
    ################################################

    # Center a molecule relative to the origin
    def recenter(self):
        if not self.atoms_read:
            self.read_atoms()

        # Define a new origin
        max = [0.0,0.0,0.0]
        min = [0.0,0.0,0.0]
        origin = [0.0,0.0,0.0]
        for i in range(0,3):
            max[i] = max(self.atoms, key = lambda a: a.coord[i])
            min[i] = min(self.atoms, key = lambda a: a.coord[i])
            origin[i] = (max[i]+min[i])/2

        print(origin)
        # Displace atoms
        for a in self.atoms:
            for i in range(0,3):
                a.coord[i] -= origin[i]

    # Compute the polarizability of a system
    def compute_alpha(self, omega=0.0, gamma=0.1088j, states=0, pr=False, usefreq=False, modefreq=0.0):
        # Make sure states have been read
        #if not self.orbs_read:
        #    self.read_orbs([])
        if not self.states_read:
            self.read_states()

        alpha = [[[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
                  [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
                  [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)]]]
        #for i in range(0,len(types)*3):
        #    alpha.append([[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
        #                  [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
        #                  [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)]])

        if states == 0 or states >= len(self.states):
            states = len(self.states) - 1

        g_state = gamma
        e_diff = modefreq / cnst.ev2cm

        for k in range(1,states+1):
            s = self.states[k]
            #print(s.energy, e_diff, omega, gamma, s.energy + omega + g_state, s.energy - e_diff - omega - g_state, s.energy - omega - g_state, s.energy - e_diff + omega + g_state)
            #print(k, s.char[0], s.achar[0],s.bchar[0],s.char[0]+s.achar[0]+s.bchar[0])

            # Compute state contribution to alpha
            #a_state = [[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)]]
            
            # Following Reimers code, convert to cm-1 here
            if usefreq == False:
                a_denom_1 = (s.energy - g_state - omega) * cnst.ev2cm
                a_denom_2 = (s.energy + g_state + omega) * cnst.ev2cm
                #print('1',a_denom_1/cnst.ev2cm, a_denom_2/cnst.ev2cm)
            else:
                # Explicitly include the difference in energies between the initial and final vibronic states
                e_diff = modefreq / cnst.ev2cm

                # Version 1: Standard Kramers-Heisenberg formula
                #a_denom_1 = (s.energy          - omega - g_state) * cnst.ev2cm
                #a_denom_2 = (s.energy - e_diff + omega + g_state) * cnst.ev2cm  #Stokes line where (s.energy - e_diff) = transition energy from s to final state

                #alpha[0][i][j] += (a_num/a_denom_1 + a_num/a_denom_2) * cnst.alph

                # Version 2: Standard Kramers-Heisenberg formula with states reversed
                a_denom_1 = (s.energy          + omega + g_state) * cnst.ev2cm
                a_denom_2 = (s.energy - e_diff - omega - g_state) * cnst.ev2cm  #Stokes line where (s.energy - e_diff) = transition energy from s to final state
                #print('2',a_denom_1/cnst.ev2cm, a_denom_2/cnst.ev2cm)

            for i in range(0,3):
                for j in range(0,3):
                    a_num = s.tdip[i]*s.tdip[j]
                    alpha[0][i][j] += (a_num/a_denom_1 + a_num/a_denom_2) * cnst.alph



        if pr == True:
            for k in range(0,len(alpha)):
                atot = '%.4f'%abs((alpha[k][0][0] + alpha[k][1][1] + alpha[k][2][2]) / 3.0)
                are = '%.4f'%((alpha[k][0][0].real + alpha[k][1][1].real + alpha[k][2][2].real) / 3.0)
                aim = '%.4f'%((alpha[k][0][0].imag + alpha[k][1][1].imag + alpha[k][2][2].imag) / 3.0)
                print('Orientationally averaged polarizability')
                print('  Magnitude (au): ', atot.rjust(10))
                print('  Real (au):      ', are.rjust(10))
                print('  Imag (au):      ', aim.rjust(10))

        # Orientationally average alpha
        #a_or = complex(0.0,0.0)
        #a_or +=  (alpha[0][0] + alpha[1][1] + alpha[2][2]) / 3.0
        #print(a_or, abs(a_or))

        return alpha

    # Find the superatomic character
    def find_superatom_char(self, types, cutoff=0.4):
        # Check if type information has been calculated
        t = False
        for i in range(0,len(types)):
            type = types[i]
            if 'SUPER' in type:
                t_ind = i
                t = True
        if t == False:
            types.append(['Orbtype','Ag','S','P'])
            t_ind = len(types) - 1
            self.read_orbs(types)
            self.read_configs(types)
            self.read_states(types)
        print('t_ind',t_ind)

        # Determine whether each orbital is superatomic
        for mo in self.mos:
            if mo.char[t_ind] > cutoff:
                mo.super = True
            else:
                mo.super = False
            #print(mo.super,mo.char[t_ind])

        # Determine the superatomic character of each excited state
        # Depends only on occupied orbital characters
        for exc in self.states:
            exc.superatom_char = 0.0
            exc.interband_char = 0.0
            for conf in exc.coeff:
                # conf[0] = configuration number; conf[1] = weight
                for occ in self.configs[conf[0]-1].occ:
                    exc.interband_char += conf[1]*self.mos[occ].char[t_ind]
                    if self.mos[occ].super == True:
                        exc.superatom_char += conf[1]

    ################################################
    #                                              #
    #  Write simple output files                   #
    #                                              #
    ################################################

    def write_orbs(self, outfilename, types=[]):
        if not self.orbs_read:
            self.read_orbs(types)

        orb = open(outfilename+'.orb','w')
        orb.write("  MO#  Energy(eV)")

        for j in types:
            orb.write((j[0][:3] + ' ' + j[1]).rjust(12))

        orb.write('\n')

        for i in range(0,len(self.mos)):
           orb.write(str(i+1).rjust( 5) + ("%.4f" % self.mos[i].energy).rjust(12))
           for j in range(0,len(types)):
               orb.write(("%.4f" % self.mos[i].char[j]).rjust(12))
           orb.write('\n')
        orb.close()

    # Write oscillator strength info
    def write_osc(self, outfilename, types=[]):
        if not self.states_read:
            self.read_states(types)
        
        osc = open(outfilename+'.osc','w')
        osc.write(" Exc#  Energy(eV)      OscStr")

        for j in types:
            osc.write((j[0][:3] + ' ' + j[1]).rjust(12))

        osc.write('\n')

        nstate = len(self.states)

        for j in range(1,nstate):
            osc.write(str(j+1).rjust(5) + ("%.4f" % self.states[j].energy).rjust(12) + ("%.5f" % self.states[j].osc).rjust(12))

            for i in range(0,len(types)):
                osc.write(("%.4f" % self.states[j].char[i]).rjust(12))

            osc.write('\n')

        osc.close()

    # Write excited-state analysis
    def write_exc(self, outfilename, types=[]):
        if not self.configs_read:
            self.read_mos(types)
            self.read_configs(types)
            self.read_states(types)
        if not self.states_read:
            self.read_states(types)

        for state in self.states:
            state.find_energy_shift(self.configs, self.mos)
            state.find_collectivity()
        self.find_superatom_char(types)

        nstate = len(self.states)

        try: 
            dip = open(outfilename+'.dip','r')
            dline = dip.readline()
            for j in range(1,nstate):
                dline = dip.readline()
                line = dline.split()
                self.states[j].dipadd = float(line[3])/100.
            dipadd = True
            dip.close()
        except IOError:
            dipadd = False

        exc = open(outfilename+'.exc','w')
        ex2 = open(outfilename+'.ex2','w')
        exc.write("#Num    Energy    OscStr Superatom   Collect CoupRange Interband")
        ex2.write("#Num    Energy    OscStr Superatom   Collect CoupRange Interband")
        if dipadd == True:
            exc.write("   Dip Add")
            ex2.write("   Dip Add")
        exc.write('\n')
        ex2.write('\n')

        for j in range(1,nstate):
            exc.write(str(j+1).rjust(4) + ("%.4f" % self.states[j].energy).rjust(10) + ("%.5f" % self.states[j].osc).rjust(10))
            exc.write(("%.4f" % self.states[j].superatom_char).rjust(10) + ("%.4f" % self.states[j].collectivity).rjust(10))
            exc.write(("%.4f" % self.states[j].weighted_std).rjust(10)   + ("%.4f" % self.states[j].interband_char).rjust(10))
            if dipadd == True:
                exc.write(("%.4f" % self.states[j].dipadd).rjust(10))
            exc.write('\n')
        exc.write('\n\n\n')


#        maxnume = min(200,len(self.states))-1
#        e_max = min(self.states[maxnume].energy,7.0)

#        self.states.sort(key=lambda state: state.weighted_std, reverse=True)
        # Write only states with non-negligible transition dipoles

#        ex2.write("#Num    Energy   Osc Str Superatom   Collect CoupRange Interband")
#        if dipadd == True:
#            ex2.write("   Dip Add")
#        ex2.write('\n')
        for j in range(1,nstate):
            if self.states[j].osc > 1e-2:# and self.states[j].superatom_char > 0.1 and self.states[j].collectivity > 2:
                ex2.write(str(j+1).rjust(4) + ("%.4f" % self.states[j].energy).rjust(10) + ("%.5f" % self.states[j].osc).rjust(10))
                ex2.write(("%.4f" % self.states[j].superatom_char).rjust(10) + ("%.4f" % self.states[j].collectivity).rjust(10))
                ex2.write(("%.4f" % self.states[j].weighted_std).rjust(10)   + ("%.4f" % self.states[j].interband_char).rjust(10))
                if dipadd == True:
                    ex2.write(("%.4f" % self.states[j].dipadd).rjust(10))
                ex2.write('\n')

        exc.close()

    # Compute the absorption spectrum
    def write_sigma(self, outfilename, types=[], max_e=8.0, e_step=0.02, gamma=0.2):
        # Make sure states have been read
        if not self.states_read:
            self.read_states(types)

        step_count = int(max_e/e_step) + 1
        sigma = []
        for i in range(0,step_count):
            sigma.append([i*e_step, 0.0])
            for j in range(0,len(types)*2):
                sigma[i].append(0.0)

        for j in self.states:
            # Add the cross-section to the appropriate bins
            if j.osc > 1e-6:
                for i in range(0,step_count):
                    current_e = sigma[i][0]
        
                    # Use Lorentzian line shape
                    phi = gamma/(((current_e - j.energy)**2 + gamma**2)*math.pi)
        
                    sigma[i][1] += phi * j.osc

                    for k in range(0,len(types)):
                        if j.char[k] > 0.0: 
                            sigma[i][2 + k*2] += phi * j.osc * j.char[k]
                        else:
                            sigma[i][3 + k*2] -= phi * j.osc * j.char[k]

        sig = open(outfilename+'.sigma','w')

        sig.write('Energy(eV)         TotAbs')
        for j in types:
            sig.write(('To ' + j[0][:3] + ' ' + j[1]).rjust(15) + ('From ' + j[0][:3] + ' ' + j[1]).rjust(15))
        sig.write('\n')

        for i in range(0,step_count):
            sig.write(("%.4f"%sigma[i][0]).rjust(10)) 
            for j in sigma[i][1:]:
                sig.write(("%.6f" % j).rjust(15))
            sig.write('\n')
        sig.close()

    def write_alpha(self, outfilename, omega=0.0, gamma=0.1088j, states=0):
        alpha = self.compute_alpha(omega, gamma, states, True)

        #a_or =  (alpha[0][0][0] + alpha[0][1][1] + alpha[0][2][2]) / 3.0

        al_out = open(outfilename + '.polariz', 'w')
        al_out.write('Omega = ' + '%.4f'%omega + ' eV   Gamma = ' + '%.4f'%gamma.imag + ' eV')
        al_out.write('\n\n')

        for k in range(0,len(alpha)):
            a_or =  (alpha[k][0][0] + alpha[k][1][1] + alpha[k][2][2]) / 3.0

            al_out.write('Real polarizability directional components (au)\n')
            al_out.write('        x           y           z\n')
            print('Real polarizability directional components (au)')
            print('        x           y           z')
            coord = 'xyz'
            for i in range(0,3):
                al_out.write(coord[i] + ('%.4f'%alpha[k][i][0].real).rjust(12) 
                             + ('%.4f'%alpha[k][i][1].real).rjust(12) 
                             + ('%.4f'%alpha[k][i][2].real).rjust(12) + '\n')
                print(coord[i] + ('%.4f'%alpha[k][i][0].real).rjust(12) 
                      + ('%.4f'%alpha[k][i][1].real).rjust(12) 
                      + ('%.4f'%alpha[k][i][2].real).rjust(12))
            al_out.write('\n')
            al_out.write('Imag polarizability directional components (au)\n')
            al_out.write('        x           y           z\n')
            print('Imag polarizability directional components (au)')
            print('        x           y           z')
            for i in range(0,3):
                al_out.write(coord[i] + ('%.4f'%alpha[k][i][0].imag).rjust(12) 
                             + ('%.4f'%alpha[k][i][1].imag).rjust(12) 
                             + ('%.4f'%alpha[k][i][2].imag).rjust(12) + '\n')
                print(coord[i] + ('%.4f'%alpha[k][i][0].imag).rjust(12) 
                      + ('%.4f'%alpha[k][i][1].imag).rjust(12) 
                      + ('%.4f'%alpha[k][i][2].imag).rjust(12))
            al_out.write('\n')

        
            al_out.write('Orientationally averaged polarizability (au)\n')
            al_out.write(('%.4f'%a_or.real).rjust(12) + ('%.4f'%a_or.imag).rjust(12) + ('%.4f'%abs(a_or)).rjust(12) + '\n\n\n')

        al_out.close()
