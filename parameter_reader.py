import sys

def read_md_params(filename):
    """
    A function to read MD parameters
    :type  filename: str
    :param filename: the name of the file contains md parameters

    return md_params which contains md parameters in a python dict
    """

    md_params = {}
    with open(filename, 'r') as fh:
        # use stack generator expressions to get rid of blank lines
        lines = (line.rstrip() for line in fh)
        lines = list(line for line in lines if line)
        for line in lines:
            comment_position = line.find('#')
            try:
                name, params = line[:comment_position].split()[:]
                md_params[name] = params
            except:
                sys.exit('Errors in reading md parameters from file %s'%filename)

    return md_params

def validate_md_params(md_params):
    """
    A function to validate MD paramters
    :type  md_params: dict
    :param md_params: a dictionary contains "number-of-steps", "time-step", "temperature",                      "output-frequency", "tau-T"
    """
    all_keys = ["number-of-steps", "time-step", "temperature", "output-frequency", "tau-T"]
    for key in md_params.keys():
        if key in all_keys:
            continue
        else:
            sys.exit('parameter %s for md is need to execute the simulation'%key)

def read_pdb(filename):
    """
    A function to read pdb data

    :type  filename: str
    :param filename: name of the file that contains structure information of input molecule

    return [box, coords, atom_name, res_name, res_number, elem, bond]
    """
    box = []
    coords = []
    atom_name = []
    res_name = []
    res_number = []
    elem = []
    bond = []
    with open(filename, 'r') as fh:
        try:
            for line in fh.readlines():
                if (line.find("ATOM") == 0 or line.find("HETATM") == 0):
                    x = float(line[31:37])
                    y = float(line[38:45])
                    z = float(line[46:53])
                    coords.append([x, y, z])
                    atom_name.append(line[12:16])
                    res_name.append(line[17:20])
                    res_number.append(int(line[22:27])-1)
                    if (len(line) >= 77):
                        elem.append(line[76:78].strip())
                    else:
                        elem.append("")
                elif (line.find("CRYST1") == 0):
                    box.append(float(line[7:15]))
                    box.append(float(line[16:24]))
                    box.append(float(line[25:33]))
                elif (line.find("CONECT") == 0):
                    bond.append([int(line[7:12])-1, int(line[13:18])-1])
        except:
            print('line: %s is not able to be processed!'%line)    

    return [box, coords, atom_name, res_name, res_number, elem, bond]

def read_force_field(filename):

    """
    A function to return force field params

    :type  filename: str
    :param filename: name of the file that contains force field information

    return ff: a dictionary that contains all force field parameters for MD simulations
    """

    bond_length = {}
    bond_force_constants = {}
    sigma = {}
    epsilon = {}
    mass = {}
    charge = {}
    ff = {}

    with open(filename, 'r') as fh:
        for line in fh.readlines():
            comment_position = line.find('#')
            params = line[:comment_position].split()
            key = params[0].strip()
            Nparam = len(params)

            if (key.find("sigma") == 0 and Nparam == 3):
                sigma[params[1].strip()] = float(params[2].strip())
            elif (key.find("epsilon") == 0 and Nparam == 3):
                epsilon[params[1].strip()] = float(params[2].strip())
            elif (key.find("bond") == 0 and Nparam == 5):
                bond  = params[1] + "-" + params[2]
                bond2 = params[2] + "-" + params[1]
                bond_length[bond] = float(params[3])
                bond_length[bond2] = float(params[3])
                bond_force_constants[bond] = float(params[4])
                bond_force_constants[bond2] = float(params[4])
            elif (key.find("mass") == 0 and Nparam == 3):
                mass[params[1].strip()] = float(params[2].strip())
            elif (key.find("charge") == 0 and Nparam == 3):
                charge[params[1].strip()] = float(params[2].strip())
            else:
                print("Unknown keyword '%s' in %s" % ( key, filename ) )
        
        ff["bond_length"]       = bond_length
        ff["bond_force_const"]  = bond_force_constants
        ff["sigma"]             = sigma
        ff["epsilon"]           = epsilon
        ff["charge"]            = charge
        ff["mass"]              = mass

    return ff