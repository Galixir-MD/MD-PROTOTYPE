#!/usr/bin/env python3

import sys, argparse
from parameter_reader import read_md_params, validate_md_params, read_pdb, read_force_field, write_pdb_frame 
from utilities import make_angles, make_exclusion, get_masses, get_temperature, compute_lambda_T, put_in_box
from md_forces import calculate_forces
from md_integrator import integrate
import time

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

    # get structure related informations from pdb
    [box, coords, atom_names, res_names, res_numbers, elems, bonds] = read_pdb(args.coordinates)

    # add angles
    bond_original = [bond for bond in bonds]
    new_bonds = make_angles(bonds)

    # generate intramolecular exclusion list
    exclusion_list = make_exclusion(len(coords), new_bonds)

    # generate random velocity, currently set all to zero
    velocities = [[0.0, 0.0, 0.0] for i in range(len(coords))]
    
    # get the force field parameters
    force_field_parameters = read_force_field(args.forcefield)

    masses = get_masses(elems, force_field_parameters["mass"])

    # Initial Temperature coupling factor
    lambda_T = 1.0

    # Open the trajectory file
    # outputfile = open(args.trajectory, "w", encoding='utf-8')

    # loop over MD steps
    with open(args.trajectory, 'w') as fh:
        for step in range(int(md_params['number-of-steps'])):
            # compute force
            time_force_start = time.time()
            epotential, forces = calculate_forces(box, coords, elems, new_bonds, exclusion_list, force_field_parameters)
            print('time to calculate force: ', time.time() - time_force_start)

            # integrate
            time_step = float(md_params["time-step"])
            time_integrate_start = time.time()
            [ekinetic, coords, velocities] = integrate(box, coords, velocities, forces, masses, time_step, lambda_T)
            print('time to integrate md: ', time.time() - time_integrate_start)

            # compute temperature
            time_compute_temp = time.time()
            T = get_temperature(len(coords), ekinetic)
            print('time to compute temperature: ', time.time() - time_compute_temp)

            # update coupling factor
            temperature, tau_T = float(md_params["temperature"]), float(md_params["tau-T"])
            lambda_T = compute_lambda_T(T, temperature, time_step, tau_T)
            # put the coordinates back in box
            put_in_box(box, res_names, coords)

            print("Step: %5d Epotential %10.3f Ekinetic %10.3f Etotal %10.3f T %7.2f lambda %.2f"%(
                step, epotential, ekinetic, ekinetic+epotential, T, lambda_T
            ))

            if (step % int(md_params['output-frequency']) == 0):
                write_pdb_frame(fh, step, box, coords, atom_names, res_names, res_numbers, elems, None)

    if (args.outcoords):
        with open(args.outcoords, 'w') as fh:
            write_pdb_frame(fh, step, box, coords, atom_names, res_names, res_numbers, elems, bond_original)