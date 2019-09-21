#!/usr/bin/env python3

import sys, argparse
from parameter_reader import read_md_params, validate_md_params, read_pdb, read_force_field 
from utilities import make_angles, make_exclusion, get_masses, get_temperature, compute_lambda_T, put_in_box
from md_forces import calculate_forces
from md_integrator import integrate

def parse_args():
    """
    parse commandline arguments
    """

    parser = argparse.ArgumentParser(
        description="running md with user defined molecule and force field",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-c',
        '--coordinates',
        dest='coordinates',
        type=str,
        default=None,
        help='name of the pdb file (must end in .pdb)'
    )

    parser.add_argument(
        '-traj',
        '--trajectory',
        dest='trajectory',
        type=str,
        default='traj.pdb',
        help='name of the output files that contains trajectory of md run, by default it saves in traj.pdb'
    )

    parser.add_argument(
        '-w',
        '--outcoords',
        dest='outcoords',
        type=str,
        default=None,
        help='name of the pdb file for writing and restarting'
    )

    parser.add_argument(
        '-p',
        '--parameters',
        dest='parameters',
        type=str,
        default=None,
        help='name of the parameter file.txt that contains #steps, time_step, T, tau, output frequency'
    )

    parser.add_argument(
        '-ff',
        '--forcefield',
        dest='forcefield',
        type=str,
        default=None,
        help='force field parameters for MD simulation'
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    print("args: ", args)
    # check commandline arguments
    assert args.coordinates.endswith('.pdb')
    if not args.parameters:
        sys.exit('Please provide MD parameters for the simulation!')
    if not args.forcefield:
        sys.exit('Please provide ff parameters for the simulation!')

    # get md parameters
    md_params = read_md_params(args.parameters)
    validate_md_params(md_params)
    print("md_params: ", md_params)

    # get structure related informations from pdb
    [box, coords, atom_names, res_names, res_numbers, elems, bonds] = read_pdb(args.coordinates)
    print("bonds: ", bonds)
    print("box: ", box)
    print("res_name: ", res_names)
    print("res_numbers: ", res_numbers)

    # add angles
    bond_original = [bond for bond in bonds]
    new_bonds = make_angles(bonds)
    print("new_bonds: ", new_bonds)

    # generate intramolecular exclusion list
    exclusion_list = make_exclusion(len(coords), new_bonds)
    print("exclusion_list", exclusion_list)

    # generate random velocity, currently set all to zero
    velocities = [[0.0, 0.0, 0.0] for i in range(len(coords))]

    # get the force field parameters
    force_field_parameters = read_force_field(args.forcefield)
    print("force_field_parameters: ", force_field_parameters)

    masses = get_masses(elems, force_field_parameters["mass"])

    print("masses: ", masses)

    # Initial Temperature coupling factor
    lambda_T = 1.0

    # loop over MD steps
    for step in range(int(md_params['number-of-steps'])):
        # compute force
        epotential, forces = calculate_forces(box, coords, elems, bonds, exclusion_list, force_field_parameters)

        # integrate
        time_step = float(md_params["time-step"])
        [ekinetic, coords, velocities] = integrate(box, coords, velocities, forces, masses, time_step, lambda_T)

        # compute temperature
        T = get_temperature(len(coords), ekinetic)

        # update coupling factor
        temperature, tau_T = float(md_params["temperature"]), float(md_params["tau-T"])
        lambda_T = compute_lambda_T(T, temperature, time_step, tau_T)

        # put the coordinates back in box
        put_in_box(box, res_names, coords)

        print("Step: %5d Epotential %10.3f Ekinetic %10.3f Etotal %10.3f T %7.2f lambda %.2f"%(
            step, epotential, ekinetic, ekinetic+epotential, T, lambda_T
        ))