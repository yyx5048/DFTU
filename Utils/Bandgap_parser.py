import os
import pandas as pd
import xml.etree.ElementTree as ET
import numpy as np

from Utils.whether_DFTU import *
from Base.Elements import transition_metal_elements
from Base.Elements import dftu_elements

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

def xml_parser():
    """
    Parse the band gap from the pwscf.xml

    Args:
    Return: (Float) Band gap in unit of ev.
    """

    xml_file = "tmp/pwscf.xml"

    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Initialize large values for CBM and VBM
    CBM = 100000
    VBM = -100000

    # Iterate through each ks_energies for the eigenvalues and occupation state
    for ks in root.iter('ks_energies'):
        for eig in ks.iter('eigenvalues'):#-for spin-polarized calculation
            eigen = np.array([s for s in eig.text.split()]).astype('float')
        for occ in ks.iter('occupations'):
            occup = np.array([s for s in occ.text.split()]).astype('float')

        # Parse the eigenvalues and occupations
        # When the occupation is less than 10^-4 consider that state to be empty.
        for i in range(len(occup)):
            if occup[i+1] <= 1e-4:
                CBMt = eigen[i+1]
                VBMt = eigen[i]
                break

        # Update the CBM and VBM values
        if CBMt < CBM:
            CBM = CBMt
        if VBMt > VBM:
            VBM = VBMt

    # Convert values to eV
    bandgap = (CBM - VBM) * 27.211396132

    return bandgap


def main():

    calc_root = os.getcwd()

    dftu_type = "All"

    DMREF_mat_list, ini_idx = read_calculation_list()

    #--Record DFTU Band gap-----------------------------------------------------
    # - If no DFTU elements in the compound, None wiil be used for DFTU_all_gap
    #   and DFTU_all_gap, depends on which atoms the Hubbard_U is applied on.
    #
    #--Error_types:
    #              Error_1: vc-relax didn't finish
    #              Error_2: hp.x didn't finish (either scfu_first or hp)
    #              Error_3: final_scf_all didn't finish
    #              Error_4: final_scf_TM didn't finish
    #
    #===========================================================================

    DFT_gap = []
    DFTU_all_gap = []#--Apply Hubbard_U on all possible elements
    DFTU_TM_gap = []#--Apply Hubbard_U on elements with d/f electrons.

    for idx, mat in enumerate(DMREF_mat_list):

        print("\n###")
        print("nubmer {} materials in the DMREF.csv data list is {}".format(idx+ini_idx,mat))

        try:
            whether_DFTU(mat, dftu_type = dftu_type)#--Check do we need DFTU

        except ValueError as e:
            print(str(e))
            if dftu_type = "All":
                DFTU_all_gap.append(None), DFTU_TM_gap.append(None)
            elif dftu_type = "TM_only":
            os.chdir(calc_root)
            continue





if __name__=="__main__":

    main()
