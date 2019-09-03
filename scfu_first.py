import os
import sys

import re
from pymatgen.core import structure
from pymatgen.ext.matproj import MPRester
from pymatgen.io.ase import AseAtomsAdaptor

import numpy as np
import subprocess
import pandas as pd

from Base.Pseudos import Pseudos
from Base.PBS_submit import pbs_submit
import Base.QE_input_settings as qe_input
from ase.io import read, write
from Utils.whether_DFTU import *
from Utils.Parser import pw_parser
from Utils.Hubbard import initialize_hubbard, insert_hubbard_block, reorder

#==============================================================================#
#  Scf calculations with DFT+U initialized (U = 1e-8)                          #
#==============================================================================#
calc_root = os.getcwd()

def read_calculation_list(fname="./Data/Final_DMREF_Materials_List.csv",start = 0, end = 250):
    """
    Obtaion the list of structure to be calculated

    Args: (Str) Excel file name
          (Int) Start index
          (Int) End index
    Return: (numpy array) A list of chemical names (Str) read from DMREF list
    """
    mat_list = pd.read_csv(fname)
    return mat_list['Formula'].values[start:end], start

def mk_scfu_folder(chem_form) -> str:
    """
    Making folder fo SCFU calculations, the validation of vc-relax is also
    performed.

    Args: Chemical formula for creating folder.
    Returns:
    """
    if os.path.isdir(chem_form):

        whether_DFTU(chem_form)#--Check do we need DFTU

        os.chdir(chem_form)
        os.chdir('vc_relax')

        res_dict = pw_parser()

        if res_dict['status'] == "DONE":
            os.chdir('../')
            os.mkdir('first_scfu')
            os.chdir('first_scfu')
            os.mkdir('tmp')
    else:
        raise FileNotFoundError("Folder is not created.")

    return res_dict

def get_relaxed_structures():
    """
    Obtain relaxed structure from vc-relax calculations
    (Validation of vc-relax is performed in mk_folder func.)

    Args:
    Returns: (ASE_Structure)
    """
    ase_S = read("../vc_relax/dftu.out", format='espresso-out', results_required = False)

    return ase_S


def ase_input_generator(ase_S, nbnd):
    """
    Using ASE write functions for generate input
    Args: (ASE Structure)
    Returns: Input files for Quantum Espresso first SCFU calculations
    """

    SCF_input = qe_input.scfu()

    pymat_S = AseAtomsAdaptor.get_structure(ase_S)

    hubbard_u_list = initialize_hubbard(pymat_S)

    #--update the SYSTEM card
    SCF_input['SYSTEM']['nbnd'] = nbnd
    SCF_input['SYSTEM'].update(insert_hubbard_block(hubbard_u_list))

    converted_pymat_S = reorder(pymat_S, hubbard_u_list)

    pseudos = Pseudos(converted_pymat_S)

    converted_ase_S = AseAtomsAdaptor.get_atoms(converted_pymat_S)

    write("dftu.in", converted_ase_S, format = "espresso-in", \
          pseudopotentials=pseudos, input_data = SCF_input, kspacing=0.04)

def main():

    #Read DMREF materials list CSV
    DMREF_mat_list, ini_idx = read_calculation_list(start = 200, end = 250)

    for idx, mat in enumearte(DMREF_mat_list):

        print("\n###")
        print("nubmer {} materials in the DMREF.csv data list is {}".format(ini_idx+idx,mat))

        try:

            results = mk_scfu_folder(mat)#--grep vc_realx results and create folder

            ase_structure = get_relaxed_structures()

            ase_input_generator(ase_structure, results['nbnd'])

            pbs_submit("pw",1)

            os.chdir(calc_root)

        except (ValueError, FileNotFoundError) as e:
            print(str(e))
            os.chdir(calc_root)
            continue

if __name__=="__main__":

    main()
