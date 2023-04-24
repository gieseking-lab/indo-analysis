import string
import math
import utils.constants as cnst
from utils.molecule import Molecule

class VibAll(object):
    def __init__(self, vibfile='nmodes.inp', prog=' ',nrsfile=' '):
        self.modes = []
        if nrsfile == ' ':
            self.get_modes(vibfile, prog)
        else:
            self.get_nrs(nrsfile)

    # Get all vibrational information
    def get_modes(self, vibfile):
        mfile = open(vibfile,'r')
        mline = mfile.readline()
        nmodes = 0
        while len(mline) > 0:
            line = mline.split()
            curr_modes = len(line)
            for i in range(0,curr_modes):
                self.modes.append(Vib(float(line[i]),nmodes+i+1))
                # Temp edit for R6G
                #self.modes.append(Vib(float(line[i]),nmodes+i+32))
                
            mline = mfile.readline()
            mline = mfile.readline()
            while len(mline) > 10:
                line = mline.split()
                for i in range(0,curr_modes):
                    self.modes[nmodes+i].add_disp([float(line[i*3+1]),float(line[i*3+2]),float(line[i*3+3])])
                mline = mfile.readline()
            nmodes += curr_modes
            mline = mfile.readline()
            mline = mfile.readline()
        mfile.close()
        print('Modes ',len(self.modes),self.modes[0].index,self.modes[-1].index)

    # Read S factors from nrs file
    def get_nrs(self, nrsfile):
        nrs = open(nrsfile,'r')
        nline = nrs.readline()
        while 'REAL' not in nline and len(nline) > 0:
            nline = nrs.readline()
        nmodes = 0
        read = False
        while read == False:
            while '=========' not in nline and len(nline) > 0:
                line = nline.split()
                self.modes.append(Vib(float(line[0]),nmodes+1))
                self.modes[-1].s_fact_r = [float(line[2])]
                nline = nrs.readline()
                line = nline.split()
                self.modes[-1].s_fact_i = [float(line[1])]
                nline = nrs.readline()
                nline = nrs.readline()
                nmodes += 1
                #print self.modes[-1].freq,self.modes[-1].s_fact_r,self.modes[-1].s_fact_i
            while 'REAL' not in nline and 'A D F   E X I T' not in nline and len(nline) > 0:
                nline = nrs.readline()
            if len(nline) == 0 or 'A D F   E X I T' in nline:
                read = True
        nrs.close()
        print('Modes ',len(self.modes),self.modes[0].index,self.modes[-1].index)


class Vib(object):
    def __init__(self, freq, index):
        self.freq = freq
        self.index = index
        self.disp = []
        self.norm = 0.0

    def add_disp(self, atdisp):
        self.disp.append(atdisp)

    # Calculate mass-weighted normalization constant for modes
    def norm_mode(self, atoms):
        # Compute mass-weighted norm301 from Jensen code
        if len(self.disp) > 0:
            init_norm = 0.0
            for j in range(0,len(atoms)):
                for k in range(0,3):
                    init_norm += self.disp[j][k]**2 * atoms[j].mass
            init_norm = math.sqrt(init_norm)
    
            # Reset norm to distance-based??
            for j in range(0,len(atoms)):
                for k in range(0,3):
                    self.norm += self.disp[j][k]**2 / init_norm**2
            #self.norm = math.sqrt(self.norm)
            #print string.rjust(str(self.index),3), ' norm ', self.norm
        else:
            self.norm = 1.0


    # Compute the normalized derivative of alpha
    def alpha_slope(self, outfilename, disp_str, omega=0.0, gamma=0.1088j, states=0, gamma2=0.0, types=[], usefreq=False, modefreq=0.0):
        if self.norm == 0.0:
            print('Error: Must compute norms of vibrational modes before polarizability derivatives')
        
        # Matches ds in Jensen code
        ds = float(disp_str)/(math.sqrt(self.norm)*cnst.bohr2ang)
        #print 'Stepsize: ',ds,float(disp_str)/cnst.bohr2ang,self.norm

        # Compute the displaced polarizabilities
        try:
            out_minus, prog = Molecule(outfilename+'_'+str(self.index)+'_-'+disp_str+'.out')
            out_plus,  prog = Molecule(outfilename+'_'+str(self.index)+'_' +disp_str+'.out')
        except IOError:
            out_minus, prog = Molecule('mode'+str(self.index)+'-'+disp_str+'.out')
            out_plus,  prog = Molecule('mode'+str(self.index)+'+'+disp_str+'.out')

        alpha_minus = out_minus.compute_alpha(omega, gamma, states, False, gamma2, types, usefreq, modefreq)
        alpha_plus  =  out_plus.compute_alpha(omega, gamma, states, False, gamma2, types, usefreq, modefreq)

        # Print changes in excited-state properties
        out_vib = open(outfilename+'_'+str(self.index)+'_' +disp_str+'.diff','w')
        
        for i in range(1,len(out_minus.states)):
            a = out_minus.states[i]
            b = out_plus.states[i]
            out_vib.write(string.rjust('%.5f'%a.energy,10) + string.rjust('%.5f'%b.energy,10) + string.rjust('%.5f'%((a.energy-b.energy)/2/float(disp_str)),10))
            out_vib.write(string.rjust('%.5f'%math.sqrt(a.osc/cnst.fosc_fact/a.energy),10) + string.rjust('%.5f'%math.sqrt(b.osc/cnst.fosc_fact/b.energy),10)) 
            out_vib.write(string.rjust('%.3f'%((math.sqrt(a.osc/cnst.fosc_fact/a.energy) - math.sqrt(b.osc/cnst.fosc_fact/b.energy))/2/float(disp_str)),10))
            out_vib.write('\n')
        out_vib.close()

        #print alpha_minus, alpha_plus

        self.alpha_diff = []
        for k in range(0,len(alpha_minus)):
            self.alpha_diff.append([[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)]])
            for i in range(0,3):
                for j in range(0,3):
                    self.alpha_diff[k][i][j] = (alpha_plus[k][i][j] - alpha_minus[k][i][j])*cnst.bohr2ang**2/(2.0*ds)
            #print string.rjust(str(self.index),3), 'adiff', string.rjust('%.5f'%(self.alpha_diff[0][0].real),14), string.rjust('%.5f'%(self.alpha_diff[1][1].real),14), string.rjust('%.5f'%(self.alpha_diff[2][2].real),14), cnst.bohr2ang**2/(2.0*ds)
            #print self.alpha_diff

    # Compute the Raman scattering factor S
    def raman_scat(self, coord):

        asq_r = []
        asq_i = []
        gsq_r = []
        gsq_i = []
        self.s_fact_r = []
        self.s_fact_i = []

        for k in range(0,len(self.alpha_diff)):
            if 'x' in coord or 'y' in coord or 'z' in coord:
                if 'x' in coord:
                    j = 0
                elif 'y' in coord:
                    j = 1
                elif 'z' in coord:
                    j = 2
    
                asq_r.append(self.alpha_diff[k][j][j].real**2)
                asq_i.append(self.alpha_diff[k][j][j].imag**2)
        
                gsq_r.append(0.0)
                gsq_i.append(0.0)
    
            else:
                asq_r.append(((self.alpha_diff[k][0][0].real + self.alpha_diff[k][1][1].real + self.alpha_diff[k][2][2].real)/3.0)**2)
                asq_i.append(((self.alpha_diff[k][0][0].imag + self.alpha_diff[k][1][1].imag + self.alpha_diff[k][2][2].imag)/3.0)**2)
        
                gsq_r.append((6.0 * (self.alpha_diff[k][0][1].real**2 + self.alpha_diff[k][1][2].real**2 + self.alpha_diff[k][2][0].real**2) \
                      + (self.alpha_diff[k][0][0].real - self.alpha_diff[k][1][1].real)**2 \
                      + (self.alpha_diff[k][1][1].real - self.alpha_diff[k][2][2].real)**2 \
                      + (self.alpha_diff[k][2][2].real - self.alpha_diff[k][0][0].real)**2) / 2.0)
                gsq_i.append((6.0*(self.alpha_diff[k][0][1].imag**2 + self.alpha_diff[k][1][2].imag**2 + self.alpha_diff[k][2][0].imag**2) \
                      + (self.alpha_diff[k][0][0].imag - self.alpha_diff[k][1][1].imag)**2 \
                      + (self.alpha_diff[k][1][1].imag - self.alpha_diff[k][2][2].imag)**2 \
                      + (self.alpha_diff[k][2][2].imag - self.alpha_diff[k][0][0].imag)**2) / 2.0)
    
            self.s_fact_r.append(45.0*asq_r[k] + 7.0*gsq_r[k])
            self.s_fact_i.append(45.0*asq_i[k] + 7.0*gsq_i[k])
    
            #print string.rjust(str(self.index),3), string.rjust('%.3f'%self.freq,8), string.rjust('%.5f'%(math.sqrt(self.s_fact_r[k]**2 + self.s_fact_i[k]**2)),14) 

    # Compute the Raman differential cross-section
    def raman_cross(self, omega):
        tempfc = 1.0e6/(1.0 - math.exp(-cnst.exparg*self.freq/cnst.temper))
        lambda0 = omega * cnst.ev2cm
        frq4th = (lambda0 - self.freq)**4

        self.crs_real = []
        self.crs_imag = []
        self.crs_tot = []

        for k in range(0,len(self.s_fact_r)):
            self.crs_real.append(tempfc * frq4th * cnst.conver * self.s_fact_r[k] / (45.0 * self.freq))
            self.crs_imag.append(tempfc * frq4th * cnst.conver * self.s_fact_i[k] / (45.0 * self.freq))
            self.crs_tot.append(self.crs_real[k] + self.crs_imag[k])

        #print string.rjust(str(self.index),3), self.crs_real

