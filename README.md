# INDO analysis scripts for MOPAC2016
Rebecca Gieseking, 2023

In summer 2020, with the help of James Stewart, I incorporated INDO/CI into the public release of MOPAC2016. This repository contains several useful analysis scripts for INDO/CI calculations. All scripts are written in Python 3 and use only standard libraries.

## Background

The INDO/S Hamiltonian (also called ZINDO or ZINDO/S) was parametrized specifically for excited-state  properties. Although there are several variations of the INDO/S parameters available, the default parameters in MOPAC2016 are the ones from Jeffrey Reimers's CNDO/INDO code. Because the INDO/S parameters are intended for excited states, geometry optimizations were not implemented. More details on the code are available in the MOPAC2016 manual at http://openmopac.net/manual/

The INDO/S excited states are computing using configuration interaction (CI), with single excitations (CIS), single and double excitations (CISD), or using multiple reference determinants (MRCI).

## Useful links

https://github.com/openmopac/mopac - Latest open-source release of MOPAC2016

http://openmopac.net/manual/ - MOPAC2016 manual containing all keywords

http://openmopac.net/manual/INDO_Examples.html - Examples of typical INDO input and output files

https://doi.org/10.1002/jcc.26455 - Publication showing the use of this code for molecular orbitals, excited states with significant double-excitation contributions, solvatochromic shifts, and nonlinear optical properties.

## Analysis scripts

This repository contains several Python3 analysis scripts for INDO/CI. All scripts (except displace.py) require a finished MOPAC INDO/CI output file as the input for the script.

Several of these files are intended for use with a single INDO output file:

| Script | Function |
| --- | --- |
  absorption.py     |  Computes a Lorentzian-broadened absorption spectrum.
  cubegen.py        |  Computes cube files in Gaussian .cub format for visualization of molecular orbitals, transition densities, and changes in electron density.
  polarizability.py |  Computes the polarizability of the system.

Several scripts are intended to assist with vibrational spectroscopy, and so require vibrational modes from another level of theory: 

| Script | Function |
| --- | --- |
  displace.py      |   Creates a MOPAC input file for a geometry displaced along a normal mode. Requires an equilibrium geometry and normal modes from another level of theory.
  raman.py         |   Computes the frequency-specific Raman intensity of vibrational modes.  Requires outputs from the equilibrium geometry and geometries displaced along the desired normal modes.

All scripts have the following options:

| Option | Usage | Default |
| --- | --- | --- |
  -h, --help     |     Prints a help file with options specific to that script
  -i, --input    |     (Required) Base file name of the MOPAC .mop or .out file. <br /> <b>NOTE:</b> The scripts automatically attempt to remove file extensions, starting with the final "." in the file name. If your file name includes periods besides the one before the extension, use the full file name (including the extension) here. If not, either the full file name or the base without the extension works. Scripts that use MOPAC output files will append a .out extension, and scripts that use MOPAC input files will append a .mop extenstion.
  -o, --output   |     Base name of the output files (files created by the scripts).  | Uses the base name from -i

## Usage of specific scripts

### absorption.py

Computes a Lorentzian-broadened absorption spectrum. 

| Option | Usage | Default |
| --- | --- | --- |
  -e, --energy    |    Maximum energy of absorption spectrum (eV)    |   8.0
  -s, --step      |    Step size for absorption spectrum (eV)        |   0.02
  -g, --gamma     |    Gamma (Lorentzian broadening factor)          |   0.1088 eV (= 0.004 au)

Output files:

1. file.orb:   List of molecular orbital energies. 
2. file.osc:   List of excited-state energies and oscillator strengths. 
3. file.sigma: Lorentzian-broadened absorption spectrum. 


### cubegen.py

Outputs file(s) in Gaussian .cub format to enable visualization of the molecular orbitals or of 
excited-state properties. The .cub files can be read by a number of visualization packages.

| Option | Usage | Default |
| --- | --- | --- |
|  -t, --type     |    Type of cube file to calculate (required) <br /> <b>Options:</b><br /> orbital (or 0):       Wavefunction for a molecular orbital<br /> ctrans  (or 1):       Transition density for a CI configuration<br /> config  (or 2):       Change in electron density upon excitation from the SCF ground state to a CI configuration<br /> trans   (or 3):       Transition density for an excited state<br /> exc     (or 4):       Change in electron density upon excitation from the SCF ground state to an excited state
|  -m, --min      |    Minimum state/orbital to compute                 | 1
|  -x, --max      |    Maximum state/orbital to compute                 | matches -m
|  -p, --param    |    Parameter file for non-default parameters (see below) | None
|  -v, --voxelsize |   Size of voxels, in Angstroms                     | 0.25 A
|  -e, --extra     |   Size of extra space surrounding the molecule     | 3.00 A

The parameter file is only necessary if MOPAC was run using `EXTERNAL=<filename>`, and if the orbital exponents were changed to non-default values in that file (any parameter that starts with Z). The parameter file used here should be the same as the parameter file used to run MOPAC.

Output files:

file_cubetype_number.cub: 'cubetype' in the file name indicates which type of wavefunction or density is in the file (orb, ctrans, config, trans, or exc). 'number' in the file name indicates which state or orbital is in the file

### polarizability.py

Computes the real and imaginary polarizability of the system at a given frequency/photon energy, using a sum-over-states (SOS) approach. 

| Option | Usage | Default |
| --- | --- | --- |
  -e, --energy    |    Energy at which polarizabilities are computed (eV)  | 0.0
  -g, --gamma     |    Lifetime (broadening) of excited states (eV)        | 0.1088 eV (= 0.004 au)
  -s, --states    |    Number of states to include in SOS expression       | All states

Output files:

file_energy.polarizability: Real and imaginary components of the polarizability at the given energy

### displace.py

Displaces an input geometry along a normal mode by a specified distance.

Requires:

1. A MOPAC input file for the equilibrium geometry. The keywords from this file will be used for the displaced geometry. Must have a .mop file extension.
2. A file containing the normal modes (by default, nmodes.inp). This file must have the following format:
     
```
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
```

This format matches the default format ADF uses to print vibrational modes. In the nmodes file, the number of lines matters, but the exact spacing on each line does not. The frequency of the mode (in cm-1) is listed first, and the x, y, and z displacements for each atom are listed beneath the divider. The divider line needs to be present, but its content does not matter at all. The xyz displacements should be normalized such that the squares of the displacements sum to 1. If more than 3 modes are included in the file, there must be exactly 3 modes listed side-by-side. Between each block of 3 modes, there must be exactly 2 blank lines.

| Option | Usage | Default |
| --- | --- | --- |
  -m, --mode         | Integer number of vibrational mode                      | 1
  -d, --displacement | RMS displacement of geometries along vibrational modes (A). Positive values are in the positive direction of the mode coordinate, and negative values are in the negative direction of the mode coordinate.  | 0.01
  -f, --modefile     | File containing normal modes                            | nmodes.inp

Output files:

file_mode_disp.mop:    MOPAC input file displaced by disp along the indicated mode

### raman.py

Computes the Raman cross-sections based on changes in the INDO/CI sum-over-states polarizabilities upon displacement along normal modes from another level of theory. 

By default, the Raman intensities are computed using all excited states in the INDO/CI output. Alternatively, the Raman intensities may be computed using a series of numbers of excited states, then averaging. This can be useful when the system is so large that the SOS expression is not fully converged. Unless the light is resonant with a specific absorbing excited state, the SOS polarizabilities usually do not fully converge until all excited states up to ~15-20 eV are included. 

Tightening the default SCF convergence criteria by 4-6 orders of magnitude can also help give better converged Raman intensities.

Requires:

  1. A MOPAC INDO/CI output file for the equilibrium geometry.
  2. A file containing the normal modes (by default, nmodes.inp). See displace.py for file format.
  3. MOPAC INDO/CI output files for geometries displaced in the positive and negative direction for all relevant vibrational modes. 

| Option | Usage | Default |
| --- | --- | --- |
  -e, --energy       | Energy at which Raman intensities are computed (eV)     | 0.0
  -g, --gamma        | Lifetime (broadening) of excited states (eV)            | 0.1088 eV ( = 0.004 au)
  -d, --displacement | RMS displacement of geometries along vibrational modes (A). Positive values are in the positive direction of the mode coordinate, and negative values are in the negative direction of the mode coordinate.  | 0.01
  -c, --coord        | Coordinate along which to compute Raman intensities     | isotropic (alternatives are x, y, z)
  -n, --nstates      | Number of states to include in SOS expression           | All states
  -f, --freqmin      | Minimum vibrational frequency of modes to compute (cm-1). This frequency is also used as the minimum frequency of the computed Raman spectrum.    | 350 
  -v, --freqmax      | Maximum vibrational frequency of modes to compute (cm-1). This frequency is also used as the maximum frequency of the computed Raman spectrum.    | 2500 
  -s, --freqstep     | Step size for the Lorenztian-broadened Raman spectrum (cm-1) | 1
  -b, --broadening   | Line width (broadening) for Lorentzian-broadened Raman spectrum (cm-1) | 20

To compute the Raman intensities averaged over a series of numbers of excited states instead of one number, use the following options. If one is used, all must be used.

| Option | Usage | Default |
| --- | --- | --- |
  -m, --nstatemin   |  Minimum number of excited states to include in SOS expression
  -x, --nstatemax   |  Maximum number of excited states to include in SOS expression
  -z, --nstatestep  |  Step size in number of excited states to include in SOS expression

Output files:

1. file.raman_stick_energy:   Stick spectrum (Raman intensities of each mode) for a given energy of light
2. file.raman_lrntz_energy:   Lorentzian-broadened Raman spectrum for a given energy of light


