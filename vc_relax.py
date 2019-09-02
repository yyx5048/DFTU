import os
import sys

import logging
from pymatgen.core import structure
from pymatgen.ext.matproj import MPRester
from pymatgen.io.pwscf import PWInput
from pymatgen.io.ase import AseAtomsAdaptor

import numpy as np
import subprocess
import pandas as pd

from Base.Pseudos import Pseudos
from Base.PBS_submit import pbs_submit
import Base.QE_input_settings as qe_input
from ase.io import read, write

#==============================================================================#
#  VC-relax from materials project relaxed structures                          #
#==============================================================================#

def record_mp_data(chemform, entry_id):
    """
    Record the mp raw data in the calculation root mk_folder

    Args: (Str) chemical formula
          (Structure) MP_ID associated with the MP structure

    Return:
    """
    if not os.path.isfile('mp_raw_data.dat'):

        with open('mp_raw_data.dat' , 'w') as f:
            f.write('Formula    MP_ID\n')

    with open('mp_raw_data.dat','a') as f:
                f.write(chemform + "    " + entry_id + '\n')
    return

def read_calculation_list(fname="./Data/Final_DMREF_Materials_List.csv",start = 200, end = 205):#200
    """
    Obtaion the list of structure to be calculated

    Args: (Str) Excel file name
          (Int) Start index
          (Int) End index
    Return: (numpy array) A list of chemical names (Str) read from DMREF list
    """
    mat_list = pd.read_csv(fname)
    return mat_list['Formula'].values[start:end]


def get_structure_from_mp(formula) -> str:
    """
    Obtain a crystal from the MP database via the API.

    Args: formula (str): A formula

    Returns: (Structure) the lowest energy structure on the Convex hull form MP
    database and its material_id
    """
    m = MPRester("N7AIm1s2v43BQ6FT")
    entries = m.get_entries(formula, inc_structure="final")
    #-- returned the computed final structure (with vasp relaxation)
    if len(entries) == 0:
        raise ValueError("This crystal structure with formula {} has been "
              "removed from the current version of MP DB".format(formula))
    elif len(entries) > 1:
        logging.warning("{} structures with formula {} found in MP, however only "
                     "the lowest energy one is returned".format(len(entries), formula))

    min_e_mp = min(entries, key = lambda e: e.energy_per_atom)

    return min_e_mp.structure, min_e_mp.entry_id

def mk_folder(chem_form) -> str:
    """
    Args: Chemical formula for creating folder.
    Returns:
    """
    if not os.path.isdir(chem_form):
        os.mkdir(chem_form)
        os.chdir(chem_form)
        os.mkdir('vc_relax')
        os.chdir('vc_relax')
        os.mkdir('tmp')

    else:
        raise FileExistsError("Duplicate folder created")

    return

def input_generator(S):
    """
    Using ASE write functions for generate raw structures
    Args: (Pymatgen Structure)
    Returns: Input files for Quantum_espress vc_relax
    """
    pseudos = Pseudos(S)
    VC_input = qe_input.VC_settings()
    ase_S = AseAtomsAdaptor.get_atoms(S)
    write("dftu.in", ase_S, format = "espresso-in", \
          pseudopotentials=pseudos, input_data = VC_input, kspacing=0.04)

def main():

    #Read DMREF materials list CSV
    DMREF_mat_list = read_calculation_list()

    for idx, mat in enumerate(DMREF_mat_list):
        print("\n###")
        print("nubmer {} materials in the DMREF.csv data list is {}".format(idx,mat))

        try:#sequence of a vc-relax operation

            mp_structure, mp_id = get_structure_from_mp(mat)

            record_mp_data(mat, mp_id)

            mk_folder(mat)#--create folder

            input_generator(mp_structure)

            pbs_submit(1)

            os.chdir('./../../')

        except (ValueError, FileExistsError) as e:
            print(str(e))
            continue

if __name__=="__main__":

    main()
