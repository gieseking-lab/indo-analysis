# Examples for INDO scripts

## 1. Benzene: Absorption spectra, orbitals, transition densities, and polarizability

To use the scripts, first run MOPAC2016 using the input file `benzene.mop` to get `benzene.out`.

To compute the absorption spectrum, use:

```
python absorption.py -i examples/benzene
```

to get the output files `benzene.orb`, `benzene.osc`, and `benzene.sigma`. 

From `benzene.osc`, the first absorption peak comes from a doubly degenerate set of states at 6.2211 eV, states 4 and 5. Looking at `benzene.out`, state 4 is primarily a linear combination of configurations 4 and 5.

```
State    4  6.2211     CI coeff  CI percent
    Config       4  -0.68863686  0.47422072
    Config       5  -0.68863995  0.47422498
    Config      33  -0.12015695  0.01443769
 Total coeff printed             0.96288339
```
 
Looking earlier in `benzene.out`, configuration 4 is a HOMO->LUMO transition (MO 15 to 16), and configuration 5 is a HOMO-1->LUMO+1 transition (MO 14 to 17)

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

This produces files `benzene_orb_14.cub`, `benzene_orb_15.cub`, `benzene_orb_16.cub`, and `benzene_orb_17.cub`. These files can be visualized using Avogadro, GaussView, and many other visualization packages. We can also generate cube files for the transition densities for configuration 4 (15->16) and 5 (14->17), to produce files `benzene_ctrans_4.cub` and `benzene_ctrans_5.cub`:

```
python cubegen.py -i examples/benzene -t ctrans -m 4 -x 5
```

We can also compute the transition density for excited state 4, to produce `benzene_trans_4.cub`:

```
python cubegen.py -i examples/benzene -t trans -m 4
```

The polarizability of benzene depends on the frequency/energy of light used to compute it. MOPAC2016 outputs the static (zero-frequency) polarizability of benzene as:

```
 Polarizability (au  ) xx=   72.32 xy=    0.00 yy=   72.32 xz=    0.00 yz=   -0.00 zz=   10.64
```

We can compute the equivalent value with `polarizability.py` by setting gamma equal to zero, since MOPAC2016 does not include any broadening. At an energy of zero, the broadening makes a very small but non-zero difference to the polarizability.

```
python polarizability.py -i examples/benzene -g 0
```

The real polarizability written to `benzene_0.0000.polarizability` is nearly identical to the ones printed by MOPAC2016, and the imaginary part is zero. At an energy relatively close to the first absorption peak, the real polarizability is much larger, and the imaginary polarizability (related to the intensity of absorption) is non-negligible:

```
python polarizability.py -i examples/benzene -e 6.0
```

On resonance with the main absorbing state, the real polarizability is small, and the imaginary polarizability is very large:

```
python polarizability.py -i examples/benzene -e 6.2211
```

## 2. CO2 Raman spectrum

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

