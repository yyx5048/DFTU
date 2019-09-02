import re
import os

def pw_parser(fname = "dftu.out"):
    """
    A simple parser to grep the results from the converged PW.x
    Args: (Str) output file name of the PW calculations
    Return: (Dict) A dictionary of possible calculations
    """
    # Section identifiers
    if not os.path.isfile("dftu.out"):
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
