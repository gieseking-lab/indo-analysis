################################################################
# Classes to store molecular properties                        #
# Used heavily by Molecule class                               #
#                                                              #
# Classes included here:                                       #
#   Atom                                                       #
#     coord  = xyz coordinates                                 #
#     elem   = 1-2 letter element symbol                       #
#     ao     = classification of atomic orbitals by type       #
#     mass   = atomic mass                                     #
#     aotype = list of AO types (S, PX, ...)                   #
#   MolOrb                                                     #
#     energy = MO eigenvalue in eV                             #
#     coeff  = list of AO coefficients                         #
#     occ    = occupation number (not used??)                  #
#     char   = list of MO character corresponding to each type #
#                in types file                                 #
#              Ranges from 0 to 1, one char element per type   #
#              Example: Ag SP character                        #
#   Config                                                     #
#     occ    = list of occupied MOs electrons are excited from #
#     vir    = list of virtual MOs electrons are excited to    #
#     energy = configuration energy in eV                      #
#     char   = change in character upon excitation             #
#                corresponding to each type in types file      #
#              Ranges from 0 to 1, one char element per type   #
#              Example: Change in Ag SP character upon         #
#                excitation                                    #
#   State                                                      #
#     energy = excited-state energy in eV (ground state = 0.0) #
#     osc    = oscillator strength (unitless)                  #
#     tdip   = transition dipole moment (Debye)                #
#     coeff  = list of [config_number,config_contribution]     #
#              Uses percent contributions = (CI coefficient)^2 #
#     char   = change in character upon excitation             #
#                corresponding to each type in types file      #
#              Ranges from 0 to 1, one char element per type   #
#              Example: Change in Ag SP character upon         #
#                excitation                                    #
#     superatom_char = superatomic character of the orbitals   #
#     collectivity   = transition inverse participation ratio  #
#                      Estimates the number of configurations  #
#                      contributing to an excited state        #
#                                                              #
################################################################

import math
import utils.constants as cnst

class Atom(object):
    def __init__(self, coord, elem):
        self.coord = coord
        self.elem = elem
        self.ao = []
        self.mass = cnst.atomic_mass[self.elem] 
        self.aotype = []

    # Add atomic orbitals
    def add_ao(self, type):
        self.ao.append(type)

    def add_aotype(self, aotype):
        self.aotype.append(aotype)

class MolOrb(object):
    def __init__(self, energy=0.0, occ=0.0):
        self.energy = energy
        self.coeff = []
        self.occ = occ
        self.char = []

    # Add more AO coefficients
    def add_coeff(self, new_coeff):
        self.coeff.append(float(new_coeff))

class Config(object):
    def __init__(self, occ, vir, char, energy=0.0):
        self.change_occ(occ)
        self.change_vir(vir)
        self.char = char
        self.energy = energy

    # Set energy of the configuration
    def change_energy(self,energy):
        self.energy = energy

    # Change occupied orbital(s) this configuration excites from
    def change_occ(self, occ):
        if occ is list:
            self.occ = occ
        else:
            self.occ = [occ]

    # Change virtual orbital(s) this configuration excites to
    def change_vir(self, vir):
        if vir is list:
            self.vir = vir
        else:
            self.vir = [vir]

class State(object):
    def __init__(self, energy=0.0, osc=0.0):
        self.energy = energy
        self.osc = osc
        self.coeff = []
        self.char = []
        self.tdip = [0.0,0.0,0.0]

        self.energy_shift = 0.0
        self.weighted_std = 0.0
        self.superatom_char = 0.0

        # For Param - track which excited states have been matched to reference data
        self.ref_matched = False

    def add_coeff(self, config, coeff):
        self.coeff.append([config, coeff])

    def calc_tdip(self, pol):
        for i in range(0,3):
            self.tdip[i] = pol[i]*math.sqrt(self.osc/(cnst.fosc_fact*self.energy))
        #print(self.tdip, self.energy * cnst.ev2cm)

    def calc_osc(self, tdip):
        self.tdip = tdip
        self.osc = (tdip[0]**2 + tdip[1]**2 + tdip[2]**2) * (cnst.fosc_fact*self.energy)
        #print self.energy, self.tdip

    # Find the collectivity index
    # Based on Adam's PEW script
    def find_collectivity(self):
        participation_sum = 0
        for conf in self.coeff:
            participation_sum += ((conf[1])**2)
        tipr = 1 / participation_sum
        self.collectivity = tipr


