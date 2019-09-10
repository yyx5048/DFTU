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
from Base.Elements import transition_metal_elements

from ase.io import read, write
from Utils.whether_DFTU import *
from Utils.Parser import hp_parser, pw_parser
from Utils.Hubbard import initialize_hubbard, insert_hubbard_block, reorder

#==============================================================================#
#  Final Scf calculations with DFT+U                                           #
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

def mk_final_scfu_folder(chem_form,dftu_type = "All") -> str:
    """
    Making folder fo SCFU calculations, the validation of vc-relax is also
    performed.

    Args: Chemical formula for creating folder.
    Returns:
    """
    if dftu_type == "All":
        hp_folder_name = "final_scf_all"
    elif dftu_type == "TM_only":
        hp_folder_name = "final_scf_TM"
    else:
        raise ValueError("Invalid DFTU type")

    if os.path.isdir(chem_form):

        whether_DFTU(chem_form, dftu_type = dftu_type)#--Check do we need DFTU

        os.chdir(chem_form)
        os.chdir('hp')

        res_dict = hp_parser()

        if res_dict['status'] == "DONE":
            os.chdir('../')

            if not os.path.isdir(hp_folder_name):
                os.mkdir(hp_folder_name)
                os.chdir(hp_folder_name)
                os.mkdir('tmp')
            else:
                raise FileExistsError("Duplicate folder created")
    else:
        raise FileNotFoundError("Folder did not create.")

    return res_dict

def get_scfu_structures():
    """
    Obtain relaxed structure from vc-relax calculations
    (Validation of vc-relax is performed in mk_folder func.)

    Args:
    Returns: (ASE_Structure)
    """
    ase_S = read("../first_scfu/dftu.in", format='espresso-in')
    nbnd = pw_parser(fname = "../first_scfu/dftu.out")['nbnd']
    return ase_S, nbnd


def ase_input_generator(ase_S, hubbard_list, nbnd, dftu_type = "All"):
    """
    Using ASE write functions for generate input
    Args: (reordered ASE Structure) ase_S
          (Int) number of bands in the calculation
          (Str) dftu_type ("All" or "TM_only")
          (Dict) Dictionary contains the hubbard DFTU list

    Returns: Input files for Quantum Espresso final SCFU calculations
    """

    SCF_input = qe_input.scfu()

    pymat_S = AseAtomsAdaptor.get_structure(ase_S)
    #-- still needs to be reordered depends on the dftu_type

    if dftu_type == "all":
        pass

    elif dftu_type == "TM_only":
        hubbard_list_iterator = hubbard_list.copy() #--avoid size change during iteration
        for kind in hubbard_list_iterator:
            if kind not in transition_metal_elements():
                hubbard_list.pop(kind)

    #--update the SYSTEM card with smearing settings
    #SCF_input['SYSTEM']['occupations'] = 'smearing'
    #SCF_input['SYSTEM']['smearing'] = 'mv'
    #SCF_input['SYSTEM']['degauss'] = 0.005

    SCF_input['SYSTEM']['nbnd'] = nbnd

    SCF_input['SYSTEM'].update(insert_hubbard_block(hubbard_list))

    converted_pymat_S = reorder(pymat_S, hubbard_list)

    pseudos = Pseudos(converted_pymat_S)

    converted_ase_S = AseAtomsAdaptor.get_atoms(converted_pymat_S)

    write("dftu.in", converted_ase_S, format = "espresso-in", \
          pseudopotentials=pseudos, input_data = SCF_input, kspacing=0.04)

def main():

    #Read DMREF materials list CSV
    dftu_type = "All"

    DMREF_mat_list, ini_idx = read_calculation_list()

    for idx, mat in enumerate(DMREF_mat_list):

        print("\n###")
        print("nubmer {} materials in the DMREF.csv data list is {}".format(idx+ini_idx,mat))

        try:

            results = mk_final_scfu_folder(mat, dftu_type = dftu_type)

            ase_structure, nbnd = get_scfu_structures()

            ase_input_generator(ase_structure, results['dftu'], nbnd, dftu_type)

            pbs_submit("pw",1)

            os.chdir(calc_root)

        except (ValueError, FileExistsError, FileNotFoundError) as e:
            print(str(e))
            os.chdir(calc_root)
            continue

if __name__=="__main__":

    main()
