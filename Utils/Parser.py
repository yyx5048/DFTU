import re
import os
import numpy as np

def pw_parser(fname = "dftu.out"):
    """
    A simple parser to grep the results from the converged PW.x
    Args: (Str) output file name of the PW calculations
    Return: (Dict) A dictionary of possible calculations
    """
    # Section identifiers
    if not os.path.isfile(fname):
        print(os.getcwd())
        raise FileNotFoundError("The calculation didn't start!!")

    parse_dict = {}

    _PW_NBND = '     number of Kohn-Sham states='

    #--critical warnings
    _PW_IONIC_MAX = 'The maximum number of steps has been reached.'
    _PW_ELECTRON_MAX = 'convergence NOT achieved after'
    _PW_DAMP_MAX = 'iterations completed, stopping'
    _PW_CRASH = '%%%%%%%%%%%%%%'
    _PW_WALLTIME_MAX = 'Maximum CPU time exceeded'

    Status = "DONE"

    indexes = {_PW_NBND: [],
               _PW_IONIC_MAX : [],
               _PW_ELECTRON_MAX : [],
               _PW_DAMP_MAX : [],
               _PW_CRASH : [],
               _PW_WALLTIME_MAX : [],
               }

    with open(fname) as f:
        pwo_lines = f.readlines()
        for idx, line in enumerate(pwo_lines):
            for identifier in indexes:
                if identifier in line:
                    indexes[identifier].append(idx)

    nbnd = re.sub(r"\s+", "", pwo_lines[indexes[_PW_NBND][0]].split('=')[1])

    parse_dict['nbnd'] = int(nbnd)

    indexes.pop(_PW_NBND) #remain indexes only contain the failed message

    for er in indexes:
        if len(indexes[er]) != 0:
            Status = "FAILED"
            print(os.getcwd())
            raise ValueError("\n!!!\nCalculation failed!!!\n!!!")

    parse_dict['status'] = Status
    return parse_dict

def hp_parser(fname = "pwscf.Hubbard_parameters.dat"):
    """
    A simple parser to grep the results from the converged PW.x
    Args: (Str) output file name of the PW calculations
    Return: (Dict) A dictionary of possible calculations
    """
    # Section identifiers
    if not os.path.isfile(fname):
        print(os.getcwd())
        raise FileNotFoundError("The calculation failed!!!")

    res_dict = {}
    dftu_parse_dict = {}

    _HEADER = '=-------------------------------------------------------------------='

    #--critical warnings
    Status = "DONE"

    indexes = {_HEADER: [],
               }

    with open(fname) as f:
        hp_lines = f.readlines()
        for idx, line in enumerate(hp_lines):
            for identifier in indexes:
                if identifier in line:
                    indexes[identifier].append(idx)

    headers = indexes[_HEADER]
    if len(headers) != 2:
        Status = "FAILED"
        raise ValueError("The header number is incorrect!!")

    for line in hp_lines[headers[0]+5: headers[1]-1]:
        l = line.split()
        if l[2] not in dftu_parse_dict.keys():
            dftu_parse_dict[l[2]] = [float(l[-1])]

        else:
            dftu_parse_dict[l[2]].append(float(l[-1]))

    for k in dftu_parse_dict:#average the Hubbard U
        dftu_parse_dict[k] = np.mean(np.array(parse_dict[k]))

    res_dict['dftu'] = dftu_parse_dict
    res_dict['status'] = Status

    return parse_dict
