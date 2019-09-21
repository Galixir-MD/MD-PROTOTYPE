import math

def sorted_pair(a, b):
    """
    A function to return a list of a and b sorted
    """

    if(a > b):
        return [b, a]
    else:
        return [a, b]

def make_angles(bonds):
    """
    A function to find the new_bonds base on shared atoms
    
    :type  bonds: list
    :param bonds: A list of bonds from pdb file

    return new_bonds: a list of bonds that have shared atoms. 
    """

    new_bonds = []
    N = len(bonds)
    for i in range(N):
        bonds[i] = sorted_pair(bonds[i][0], bonds[i][1])

    for i in range(N):
        bond = bonds[i]
        for j in range(i+1, N):
            new_bond = bonds[j]
            if(bond[0] == new_bond[0] and bond[1] != new_bond[1]):
                new_bonds.append(sorted_pair(bond[1], new_bond[1]))
            elif(bond[1] == new_bond[0] and bond[0] != new_bond[1]):
                new_bonds.append(sorted_pair(bond[0], new_bond[1]))
            elif (bond[0] == new_bond[1] and bond[1] != new_bond[0]):
                new_bonds.append(sorted_pair(bond[1], new_bond[0]))
            elif (bond[1] == new_bond[1] and bond[0] != new_bond[0]):
                new_bonds.append(sorted_pair(bond[0], new_bond[0]))
    
    return sorted(new_bonds)

def make_exclusion(N, bonds):
    """
    A function to generate exclusion list for MD

    :type  N: int
    :param N: number of atoms

    :type  bonds: list
    :param bonds: a list of bonds

    return exclusion_list: a list of atoms, base on the connection. 
    """

    exclusion_list = [[] for i in range(N)]
    for bond in bonds:
        a1, a2 = bond[0], bond[1]
        exclusion_list[a1].append(a2)
        exclusion_list[a2].append(a1)

    return exclusion_list

def get_masses(elems, mass):
    """
    A function to return mass for all elements in the system

    :type  elems: list
    :param elems: a list of 
    """
    masses = []
    for ele in elems:
        if ele in mass:
            masses.append(mass[ele])
        else:
            print("No mass for element '%s'"%ele)

    return masses

def get_temperature(natoms, ekinetic):
    """
    A function to compute T using Ekin = 3/2 * natoms * kB * T
    """
    kB = 0.00831415
    return (2 * ekinetic) / (3 * natoms * kB)

def compute_lambda_T(T, T_reference, time_step, tau_T):
    """
    A function to compute Berendsen temperature coupling, GROMACS 5.1 manual, Eqn. 3.45
    """
    if (T == 0 or tau_T == 0):
        return 1
    
    return math.sqrt(1 + (time_step/tau_T)*(T_reference/T -1))

def put_in_box(box, resnr, coords):
    """
    A function to put coordinates back into back
    """
    N      = len(coords)
    cgcm   = []
    old    = -1
    invres = []
    for i in range(N):
        if (resnr[i] != old):
            cgcm.append([ 0.0, 0.0, 0.0 ])
            invres.append([])
            old = resnr[i]
        for m in range(3):
            cgcm[len(cgcm)-1][m] += coords[i][m]
        invres[len(invres)-1].append(i)
    N = len(cgcm)
    for i in range(N):
        for m in range(3):
            cgcm[i][m] /= len(invres[i])
    for i in range(N):
        for m in range(3):
            if (cgcm[i][m] > box[m]):
                for k in invres[i]:
                    coords[k][m] -= box[m]
            if (cgcm[i][m] <= 0):
                for k in invres[i]:
                    coords[k][m] += box[m]