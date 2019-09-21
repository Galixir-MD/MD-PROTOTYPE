import numpy as np
import math
import sys


def distance_pbc(box, x1, x2):
    dx = []
    for m in range(3):
        dx.append(x1[m] - x2[m])
        while (dx[m] > box[m]/2):
            dx[m] -= box[m]
        while (dx[m] <= -box[m]/2):
            dx[m] += box[m]
    return dx

def calculate_bonded_forces(box, coords, elems, bonds, force, bond_length, force_constant):
    """
    A function to calculate bonded forces base on the bond length and force constant
    """
    energy = 0
    for bond in bonds:
        dx = distance_pbc(box, coords[bond[0]], coords[bond[1]])
        dx2 = np.dot(dx, dx)
        dx1 = math.sqrt(dx2)
        bond_name = elems[bond[0]] + '-' + elems[bond[1]]
        if(not bond_name in bond_length or not bond_name in force_constant):
            sys.exit("quitting program because of unknown %s bond"%bond_name)
        ddx = dx1 - bond_length[bond_name]
        energy += 0.5 * force_constant[bond_name]*ddx**2
        # force on one vector
        force_a = -force_constant[bond_name] * ddx / dx1
        # on xyz vectors
        for m in range(3):
            force[bond[0]][m] += force_a * dx[m]
            force[bond[1]][m] -= force_a * dx[m]

    return energy, force

def calculate_nonbonded_forces(box, coords, elems, exclude, force, sigma, epsilon, charge):
    """
    A function to calculate nonbonded interaction:
    coulomb's law: F = k * (q1 * q2) / r^2, where k = 1/(4pi * epslion_0), where epslion_0 is the electrical
                   permittivity of free space (vacuum), 8.8541487817 * 10 ^ -12
    van der waals: LJ potential, 4 * epslion_ij * [(sigma/r)^12-(sigma/r)^6]
    """
    energy = 0
    k_coulomb = 138.1
    # half of the box size should be larger than the cutoff
    cutoff = 0.9 * min(0.5*box[0], 0.5*box[1], 0.5*box[2])
    for i in range(len(coords)):
        for j in range(i+1, len(coords)):
            if (not j in exclude[i]):
                dx = distance_pbc(box, coords[i], coords[j])
                dx2 = np.dot(dx, dx)
                dx1 = math.sqrt(dx2)
                if (dx2 < cutoff**2):
                    # calculate coulomb energy and force
                    coulomb_energy = charge[elems[i]] * charge[elems[j]] * k_coulomb / dx1    
                    coulomb_force = -coulomb_energy / dx1
                    # calculate vdw energy and force
                    epsilon_ij = math.sqrt(epsilon[elems[i]]*epsilon[elems[j]])
                    if (epsilon_ij > 0):
                        sigma_ij = 0.5*(sigma[elems[i]] + sigma[elems[j]])
                        sr = sigma_ij / dx1
                        sr6 = sr**6
                        sr12 = sr**12
                        vdw_energy = 4 * epsilon_ij * (sr12 - sr6)
                        vdw_force = 4 * epsilon_ij * ((sr6*sr)**2 - sr6*sr**2)
                    else:
                        vdw_energy, vdw_force = 0, 0
                    energy += (coulomb_energy + vdw_energy) 
                    for m in range(3):
                        dfm = (coulomb_force + vdw_force) * dx[m]
                        force[i][m] += dfm
                        force[j][m] -= dfm
    
    return energy, force

def calculate_forces(box, coords, elems, bonds, exclude, ff):

    # initialize forces
    force = [[0.0, 0.0, 0.0] for i in range(len(coords))]

    bonded_energy, force = calculate_bonded_forces(box, coords, elems, bonds, force, ff["bond_length"], ff["bond_force_const"])
    nonbonded_energy, force = calculate_nonbonded_forces(box, coords, elems, exclude, force, ff["sigma"], ff["epsilon"], ff["charge"])

    return bonded_energy + nonbonded_energy, force