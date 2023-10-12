################################################################
# Functions used to generate cube files                        #
#                                                              #
################################################################

import utils.constants as cnst
import math

# Get the number of voxels
def get_vnum(ats, gap, extra):
    nvox = [0,0,0]
    
    minxyz = [ats[0].coord[0]/cnst.bohr2ang, ats[0].coord[1]/cnst.bohr2ang, ats[0].coord[2]/cnst.bohr2ang]
    maxxyz = [ats[0].coord[0]/cnst.bohr2ang, ats[0].coord[1]/cnst.bohr2ang, ats[0].coord[2]/cnst.bohr2ang]
    
    for i in range(1,len(ats)):
        for j in range(0,3):
            if ats[i].coord[j]/cnst.bohr2ang > maxxyz[j]:
                maxxyz[j] = ats[i].coord[j]/cnst.bohr2ang
            elif ats[i].coord[j]/cnst.bohr2ang < minxyz[j]:
                minxyz[j] = ats[i].coord[j]/cnst.bohr2ang
    
    for j in range(0,3):
        nvox[j] = int((maxxyz[j] + 2*extra - minxyz[j])/gap)
        minxyz[j] -= extra
    return nvox, minxyz

# Build the empty voxel matrix
def build_vox(nvox):
    vox = []
    for i in range(0,nvox[0]):
        plane = []
        for j in range(0,nvox[1]):
            line = []
            for k in range(0,nvox[2]):
                line.append(0.0)
            plane.append(line)
        vox.append(plane)
    return vox

# Get range of voxels needed for a particular atom
def get_vrange(m, vdist, nvox, minxyz, gap):
    minvox = [0,0,0]
    maxvox = [nvox[0],nvox[1],nvox[2]]
    for i in range(0,3):
        minvox[i] = max(minvox[i],int((m.coord[i]/cnst.bohr2ang-vdist-minxyz[i])/gap))
        maxvox[i] = min(maxvox[i],int((m.coord[i]/cnst.bohr2ang+vdist-minxyz[i])/gap))
    return minvox, maxvox

# Calculate the displacement between the atom and voxel
def get_disp(acoord,vcoord):
    dx = acoord[0]/cnst.bohr2ang-vcoord[0]
    dy = acoord[1]/cnst.bohr2ang-vcoord[1]
    dz = acoord[2]/cnst.bohr2ang-vcoord[2]
    rad = math.sqrt(dx**2 + dy**2 + dz**2)
    # Avoid division by zero errors
    if rad < 0.0001: rad = 0.0001
    return [dx,dy,dz,rad]


# Get the r and y coefficients for an orbital at a position
def get_rycoeff(aotype, at_par,acoord,vcoord):
    disp = get_disp(acoord,vcoord)
    n, szeta, pzeta, dzeta1, dzeta2, dcoef1, dcoef2 = at_par[0], at_par[1], at_par[2], at_par[3], at_par[4], at_par[5], at_par[6]
    
    # Occupied orbital coefficients
    if 'S' in aotype:
        r = (2*szeta)**n * math.sqrt(2*szeta/math.factorial(2*n)) * disp[3]**(n-1) * math.exp(-szeta*disp[3])
        y = math.sqrt(1.0/math.pi) / 2.0
    elif 'P' in aotype:
        # s or p orbital
        r = (2*pzeta)**n * math.sqrt(2*pzeta/math.factorial(2*n)) * disp[3]**(n-1) * math.exp(-pzeta*disp[3])
        if 'X' in aotype:
            y = -math.sqrt(3.0/(4.0 * math.pi)) * disp[0]/disp[3]
        elif 'Y' in aotype:
            y = -math.sqrt(3.0/(4.0 * math.pi)) * disp[1]/disp[3]
        else:
            y = -math.sqrt(3.0/(4.0 * math.pi)) * disp[2]/disp[3]
    else:
        # Account for both d components
        r1 = (2*dzeta1)**(n-1) * math.sqrt(2*dzeta1/math.factorial(2*n-2)) * disp[3]**(n-2) * math.exp(-dzeta1*disp[3])
        r2 = (2*dzeta2)**(n-1) * math.sqrt(2*dzeta2/math.factorial(2*n-2)) * disp[3]**(n-2) * math.exp(-dzeta2*disp[3])
        r = r1*dcoef1 + r2*dcoef2
        if 'Z2' in aotype:
            y = math.sqrt(5.0/math.pi)/4.0 * (-disp[0]**2 - disp[1]**2 + 2*disp[2]**2)/disp[3]**2
        elif 'X2' in aotype:
            y = math.sqrt(15.0/math.pi)/4.0 * (disp[0]**2 - disp[1]**2)/disp[3]**2
        elif 'XY' in aotype:
            y = math.sqrt(15.0/math.pi)/2.0 * (disp[0]*disp[1])/disp[3]**2
        elif 'XY' in aotype:
            y = math.sqrt(15.0/math.pi)/2.0 * (disp[0]*disp[2])/disp[3]**2
        else:
            # YZ
            y = math.sqrt(15.0/math.pi)/2.0 * (disp[1]*disp[2])/disp[3]**2

    return r, y

# Generate the cube values for a molecular orbital
def gen_cube_orb(mo, at_types, atoms, mos, params, extra, nvox, minxyz, gap):
    print('Generating cube for molecular orbital ', mo+1)

    vox = build_vox(nvox)
    a = 0
    # Loop over atoms
    for m in atoms:
        print('  Atom ',atoms.index(m)+1)
        minvox, maxvox = get_vrange(m, extra, nvox, minxyz, gap)
        at_index = at_types.index(m.elem)
        for orb in m.aotype:
            for i in range(minvox[0],maxvox[0]):
                for j in range(minvox[1],maxvox[1]):
                    for k in range(minvox[2],maxvox[2]):
                        r, y = get_rycoeff(orb, params[at_index], m.coord, [minxyz[0]+i*gap,minxyz[1]+j*gap,minxyz[2]+k*gap])
                        vox[i][j][k] += (r*y) * mos[mo].coeff[a]
            a += 1
    return vox

# Generate the cube values for a configuration transition density
def gen_cube_ctrans(config, at_types, atoms, mos, configs, params, extra, nvox, minxyz, gap):
    print('Generating cube for configuration ', config)

    vox = build_vox(nvox)

    # Set up for single excitations only
    occ = configs[config-1].occ[0]
    vir = configs[config-1].vir[0]
    print(occ+1,vir+1)

    if occ == vir:
        print('Occupied and virtual orbital are the same; configuration does not involve a transition')
        return None

    a = 0
    # Loop over atoms
    for m in atoms:
        print('  Atom ',atoms.index(m)+1)
        minvox, maxvox = get_vrange(m, extra, nvox, minxyz, gap)
        at_index = at_types.index(m.elem)
        for i in range(minvox[0],maxvox[0]):
            for j in range(minvox[1],maxvox[1]):
                for k in range(minvox[2],maxvox[2]):
                    for orbnum1 in range(0,len(m.aotype)):
                        orba = m.aotype[orbnum1]
                        ra, ya = get_rycoeff(orba, params[at_index], m.coord, 
                                             [minxyz[0]+i*gap,minxyz[1]+j*gap,minxyz[2]+k*gap])
                        for orbnum2 in range(0,len(m.aotype)):
                            orbb = m.aotype[orbnum2]
                            rb, yb = get_rycoeff(orbb, params[at_index], m.coord, 
                                                 [minxyz[0]+i*gap,minxyz[1]+j*gap,minxyz[2]+k*gap])
                            vox[i][j][k] += (ra*ya*rb*yb) * (mos[vir].coeff[a+orbnum1] * mos[occ].coeff[a+orbnum2])
        a += len(m.aotype)
    return vox

# Generate the cube values for a configuration change in density
def gen_cube_config(config, at_types, atoms, mos, configs, params, extra, nvox, minxyz, gap):
    print('Generating cube for configuration ', config)

    vox = build_vox(nvox)

    # Set up for single excitations only
    occ = configs[config-1].occ[0]
    vir = configs[config-1].vir[0]
    print(occ+1,vir+1)

    if occ == vir:
        print('Occupied and virtual orbital are the same; configuration does not involve a transition')
        return None

    a = 0
    # Loop over atoms
    for m in atoms:
        print('  Atom ',atoms.index(m)+1)
        minvox, maxvox = get_vrange(m, extra, nvox, minxyz, gap)
        at_index = at_types.index(m.elem)
        for i in range(minvox[0],maxvox[0]):
            for j in range(minvox[1],maxvox[1]):
                for k in range(minvox[2],maxvox[2]):
                    for orbnum in range(0,len(m.aotype)):
                        orb = m.aotype[orbnum]
                        r, y = get_rycoeff(orb, params[at_index], m.coord, 
                                           [minxyz[0]+i*gap,minxyz[1]+j*gap,minxyz[2]+k*gap])
                        vox[i][j][k] += (r*y)**2 * (mos[vir].coeff[a+orbnum]**2 - mos[occ].coeff[a+orbnum]**2)
        a += len(m.aotype)
    return vox

# Generate the cube values for an excited state transition density
def gen_cube_trans(ex, at_types, atoms, mos, configs, states, params, extra, nvox, minxyz, gap):
    print('Generating cube for excited state ', ex)
    vox = build_vox(nvox)
    
    # Pre-sum over coefficients to save looping time later
    coeff = []
    a = 0
    for m in atoms:
        #print '  Atom ',out.atoms.index(m), m.aotype
        for orbnum1 in range(0,len(m.aotype)):
            coeff.append([])
            for orbnum2 in range(0,len(m.aotype)):
                coeff[-1].append(0)
                for config in states[ex-1].coeff:
                    occ = configs[config[0]-1].occ[0]
                    vir = configs[config[0]-1].vir[0]
                    if occ != vir:
                        coeff[-1][orbnum2] += (mos[vir].coeff[a+orbnum1] * mos[occ].coeff[a+orbnum2]) * config[1]
                    #print config[0], occ, vir, out.mos[vir].coeff[a+orbnum1], out.mos[occ].coeff[a+orbnum2], config[1]
            #print coeff[-1]
        a += len(m.aotype)

    a = 0
    # Loop over atoms
    for m in atoms:
        print('  Atom ',atoms.index(m)+1)
        minvox, maxvox = get_vrange(m, extra, nvox, minxyz, gap)
        at_index = at_types.index(m.elem)

        for i in range(minvox[0],maxvox[0]):
            for j in range(minvox[1],maxvox[1]):
                for k in range(minvox[2],maxvox[2]):
                    for orbnum1 in range(0,len(m.aotype)):
                        orba = m.aotype[orbnum1]
                        ra, ya = get_rycoeff(orba, params[at_index], m.coord, 
                                             [minxyz[0]+i*gap,minxyz[1]+j*gap,minxyz[2]+k*gap])
                        for orbnum2 in range(0,len(m.aotype)):
                            orbb = m.aotype[orbnum2]
                            rb, yb = get_rycoeff(orbb, params[at_index], m.coord,
                                                 [minxyz[0]+i*gap,minxyz[1]+j*gap,minxyz[2]+k*gap])
                            vox[i][j][k] += (ra*ya*rb*yb) * coeff[a+orbnum1][orbnum2]
        a += len(m.aotype)
    return vox

# Generate the cube values for an excited state change in density
def gen_cube_exc(ex, at_types, atoms, mos, configs, states, params, extra, nvox, minxyz, gap):
    print('Generating cube for excited state ', ex)
    vox = build_vox(nvox)
    
    # Pre-sum over coefficients to save looping time later
    coeff = []
    for i in range(0,len(mos)):
        coeff.append(0)
        for config in states[ex-1].coeff:
            occ = configs[config[0]-1].occ[0]
            vir = configs[config[0]-1].vir[0]
            if occ != vir:
                coeff[-1] += (mos[vir].coeff[len(coeff)-1]**2 - mos[occ].coeff[len(coeff)-1]**2) * config[1]
    #print coeff
    
    a = 0
    # Loop over atoms
    for m in atoms:
        print('  Atom ',atoms.index(m)+1)
        minvox, maxvox = get_vrange(m, extra, nvox, minxyz, gap)
        at_index = at_types.index(m.elem)

        for i in range(minvox[0],maxvox[0]):
            for j in range(minvox[1],maxvox[1]):
                for k in range(minvox[2],maxvox[2]):
                    for orbnum in range(0,len(m.aotype)):
                        orb = m.aotype[orbnum]
                        r, y = get_rycoeff(orb, params[at_index], m.coord, 
                                           [minxyz[0]+i*gap,minxyz[1]+j*gap,minxyz[2]+k*gap])
                        vox[i][j][k] += (r*y)**2 * coeff[a+orbnum]
        a += len(m.aotype)
    return vox

# Write the cube file
def write_cub(infilename, dtype, curr, vox, ats, at_types, params, nvox, minxyz, gap):
    if not vox:
        return
    if dtype == 0:
        cub = open(infilename+'_orb_'+str(curr)+'.cub','w')
        cub.write(infilename + ' orbital ' + str(curr) + '\n\n')
    elif dtype == 1:
        cub = open(infilename+'_ctrans_'+str(curr)+'.cub','w')
        cub.write(infilename + ' config_tdens ' + str(curr) + '\n\n')
    elif dtype == 2:
        cub = open(infilename+'_config_'+str(curr)+'.cub','w')
        cub.write(infilename + ' config_dens ' + str(curr) + '\n\n')
    elif dtype == 3:
        cub = open(infilename+'_trans_'+str(curr)+'.cub','w')
        cub.write(infilename + ' exc_tdens ' + str(curr) + '\n\n')
    else:
        cub = open(infilename+'_exc_'+str(curr)+'.cub','w')
        cub.write(infilename + ' exc_dens ' + str(curr) + '\n\n')

    cub.write(str(len(ats)).rjust(5) + ('%.6f'%(minxyz[0])).rjust(12) + ('%.6f'%(minxyz[1])).rjust(12) + ('%.6f'%(minxyz[2])).rjust(12) + '\n')
    cub.write(str(nvox[0]).rjust(5) + ('%.6f'%gap).rjust(12) + ('%.6f'%0.0).rjust(12) + ('%.6f'%0.0).rjust(12) + '\n')
    cub.write(str(nvox[1]).rjust(5) + ('%.6f'%0.0).rjust(12) + ('%.6f'%gap).rjust(12) + ('%.6f'%0.0).rjust(12) + '\n')
    cub.write(str(nvox[2]).rjust(5) + ('%.6f'%0.0).rjust(12) + ('%.6f'%0.0).rjust(12) + ('%.6f'%gap).rjust(12) + '\n')
    for m in ats:
        cub.write(str(cnst.atomic_symbol[m.elem]).rjust(5) + '    0.000000')
        #cub.write(str(params[at_types.index(m.elem)][0]).rjust(5) + '    0.000000')
        cub.write(('%.6f'%(m.coord[0]/cnst.bohr2ang)).rjust(12) 
                  + ('%.6f'%(m.coord[1]/cnst.bohr2ang)).rjust(12) 
                  + ('%.6f'%(m.coord[2]/cnst.bohr2ang)).rjust(12) + '\n')

    i = 0
    j = 0
    k = 0
    for n in range(0,nvox[0]*nvox[1]*nvox[2]):
        cub.write(('%.6f'%vox[i][j][k]).rjust(13))
        if n%6 == 5:
            cub.write('\n')

        # Cycle indices
        k += 1
        if k >= nvox[2]:
            k = 0
            j += 1
        if j >= nvox[1]:
            j = 0
            i += 1
    cub.close()
