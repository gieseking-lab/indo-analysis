# INDO analysis scripts for MOPAC2016
Rebecca Gieseking, 2024

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

Note: For excited states, the ground state is defined as state 1, so the first excited state is state 2.

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

This format matches the default format ADF uses to print vibrational modes. In the nmodes file, the number of lines matters, but the exact spacing on each line does not. The frequency of the mode (in cm-1) is listed first, and the x, y, and z displacements for each atom are listed beneath the divider. The divider line needs to be present, but its content does not matter at all. If more than 3 modes are included in the file, there must be exactly 3 modes listed side-by-side. Between each block of 3 modes, there must be exactly 2 blank lines.

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
  -f, --modefile     | File containing normal modes                            | nmodes.inp
  -g, --gamma        | Lifetime (broadening) of excited states (eV)            | 0.1088 eV ( = 0.004 au)
  -d, --displacement | RMS displacement of geometries along vibrational modes (A). Positive values are in the positive direction of the mode coordinate, and negative values are in the negative direction of the mode coordinate.  | 0.01
  -c, --coord        | Coordinate along which to compute Raman intensities     | isotropic (alternatives are x, y, z)
  -n, --nstates      | Number of states to include in SOS expression           | All states
  -v, --freqmin      | Minimum vibrational frequency of modes to compute (cm-1). This frequency is also used as the minimum frequency of the computed Raman spectrum.    | 350 
  -w, --freqmax      | Maximum vibrational frequency of modes to compute (cm-1). This frequency is also used as the maximum frequency of the computed Raman spectrum.    | 2500 
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

## Examples

### Benzene: Absorption spectra, orbitals, transition densities, and polarizability

To use the scripts, we first run MOPAC2016 using the input file `benzene.mop` to get `benzene.out`.

To compute the absorption spectrum, we then use:

```
python absorption.py -i examples/benzene
```

to get the output files `benzene.orb`, `benzene.osc`, and `benzene.sigma`. 

From `benzene.osc`, the first absorption peak comes from a doubly degenerate set of states at 6.2211 eV, states 4 and 5. Looking at `benzene.out`, state 4 is primarily a linear combination of configurations 4 and 5, with a smaller contribution from configuration 33 that we will ignore. The three printed configurations sum to 96.29% of the total contributions to this state; if `WRTCONF` was set to a lower value in `benzene.mop`, more of the minor contributions would be printed.

```
State    4  6.2211     CI coeff  CI percent
    Config       4  -0.68863686  0.47422072
    Config       5  -0.68863995  0.47422498
    Config      33  -0.12015695  0.01443769
 Total coeff printed             0.96288339
```
 
Looking earlier in `benzene.out`, configuration 4 is a HOMO->LUMO transition (MO 15 to 16), and configuration 5 is a HOMO-1->LUMO+1 transition (MO 14 to 17):

```
 CI excitations=  226:        =226
 The lowest  226 spin-adapted configurations of multiplicity=  1
     sym   eV   cm**-1 -dets- dipole oscilator X FRAG ....Excitations named from first reference determinate
                       tot  #  Debye  strength
 
   1      0.000      0.  1  1  0.0000 0.000000 0  0 (    )->(       )
   2      5.802  46796.  1  1  0.0000 0.495181 1  1 (  15)->(  17   )
   3      5.802  46796.  1  1  0.0000 0.495181 1  1 (  14)->(  16   )
   4      6.483  52293.  1  1  0.0000 0.553346 1  1 (  15)->(  16   )
   5      6.483  52293.  1  1  0.0000 0.553346 1  1 (  14)->(  17   )
```

Based on this information, we can use `cubegen.py` to generate cube files of MOs 14-17:

```
python cubegen.py -i examples/benzene -t orbital -m 14 -x 17
```

This produces files `benzene_orb_14.cub`, `benzene_orb_15.cub`, `benzene_orb_16.cub`, and `benzene_orb_17.cub`. These files can be visualized using VMD, Avogadro, GaussView, and many other visualization packages. Orbitals 14 and 15 have the expected shapes for the doubly degenerate HOMO of benzene, and orbitals 16 and 17 have the expected shapes for the doubly degenerate LUMO.

We can also generate cube files for the transition densities for configuration 4 (15->16) and 5 (14->17), to produce files `benzene_ctrans_4.cub` and `benzene_ctrans_5.cub`:

```
python cubegen.py -i examples/benzene -t ctrans -m 4 -x 5
```

Based on visualizing these transition densities, configuration 4 has negative transition density on two carbons with positive y coordinates and positive transition density on two carbons with negative y coordinates. Configuration 5 has alternating positive and negative transition densities around the ring, but the largest positive and negative transition densities are on the carbon atoms with the largest +/- y coordinates.

We can also compute the transition density for excited state 4, to produce `benzene_trans_4.cub`:

```
python cubegen.py -i examples/benzene -t trans -m 4
```

Because of how the transition densities from the two configurations combine, the overall transition density for this state is negative for all carbon atoms with positive y coordinates and positive for all carbon atoms with negative y coordinates. This implies that this state should have a large transition dipole moment along the y axis, which can be seen in `benzene.out` as a large oscillator strength with a polarization along y:

```
  CI trans.  energy frequency wavelength oscillator ---------polarization---------    dipole   ------components----- 
 st.  symm.    eV      cm-1       nm      strength       x          y         z       moment     x       y       z 

   2      4.7082902    37975.    263.33  0.000000                                   0.000000   0.000  -0.000  -0.000
   3      5.3612330    43241.    231.26  0.000000                                   0.000000   0.000  -0.000   0.000
   4      6.2211294    50177.    199.29  0.872452  -0.000049   1.000000  -0.000000  0.000000   0.000  -0.000   0.000
```

We can also compute the polarizability of benzene, which depends on the frequency/energy of light at which we are computing the polarizability. MOPAC2016 outputs the static (zero-frequency) polarizability of benzene as:

```
 Polarizability (au  ) xx=   72.32 xy=    0.00 yy=   72.32 xz=    0.00 yz=   -0.00 zz=   10.64
```

We can compute the equivalent value with `polarizability.py` by setting gamma equal to zero, since MOPAC2016 does not include any broadening. At an energy of zero, the broadening makes a very small but non-zero difference to the polarizability.

```
python polarizability.py -i examples/benzene -g 0
```

The real polarizability written to `benzene_0.0000.polarizability` is nearly identical to the ones printed by MOPAC2016, and the imaginary part is zero. 

We can also compute the polarizability at an energy relatively close to the first absorption peak: 

```
python polarizability.py -i examples/benzene -e 6.0
```

At 6 eV, the output file shows that real polarizability is much larger, and the imaginary polarizability (related to the intensity of absorption) is non-negligible:

On resonance with the main absorbing state, the real polarizability is small, and the imaginary polarizability is very large:

```
python polarizability.py -i examples/benzene -e 6.2211
```

### CO2: Raman spectrum

Before starting, we need a MOPAC input with the equilibrium geometry of CO2 (`co2.mop`) and a file containing the vibrational normal modes (`co2_nmodes.inp`). These scripts compute the Raman intensities based on the numerical derivative of the polarizability with respect to displacement. That means we need to generate input files for the displaced geometries using `displace.py` for positive and negative displacements along each vibrational mode:

```
python displace.py -i examples/co2 -m 1 -d  0.01 -f examples/co2_nmodes.inp
python displace.py -i examples/co2 -m 1 -d -0.01 -f examples/co2_nmodes.inp
python displace.py -i examples/co2 -m 2 -d  0.01 -f examples/co2_nmodes.inp
python displace.py -i examples/co2 -m 2 -d -0.01 -f examples/co2_nmodes.inp
python displace.py -i examples/co2 -m 3 -d  0.01 -f examples/co2_nmodes.inp
python displace.py -i examples/co2 -m 3 -d -0.01 -f examples/co2_nmodes.inp
python displace.py -i examples/co2 -m 4 -d  0.01 -f examples/co2_nmodes.inp
python displace.py -i examples/co2 -m 4 -d -0.01 -f examples/co2_nmodes.inp
```

Each of these 8 commands generates 1 MOPAC input file. We then need to run MOPAC2016 with the equilibrium geometry input file, plus all 8 displaced geometries.

Once we have all 9 MOPAC output files, we can compute the Raman intensities, starting at zero frequency:

```
python raman.py -i examples/co2 -f examples/co2_nmodes.inp
``` 

As expected, only the asymmetric stretch at 1306.996 is Raman-active (has a Raman intensity not within rounding error of zero, as shown in `co2.raman_stick_0.0`). Because of this, there is only one peak in the Raman spectrum in `co2.raman_lrtz_0.0`, around 1307 cm-1.

We can also compute the Raman intensity at 9.0 eV, ~0.4 eV below the first absorption peak of CO2:

```
python raman.py -i examples/co2 -f examples/co2_nmodes.inp -e 9.0
```

The Raman intensity of the Raman-active mode is enhanced by > 10 orders of magnitude because of the near-resonance with the excited state. The Raman intensity very close to resonance with an absorbing state would be meaningless because light of that energy will be absorbed instead of scattered.


