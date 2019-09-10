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

def xml_parser(PATH):
    """
    Parse the band gap from the pwscf.xml

    Args: (Str) PATH to the calcualtion contained dictionary
    Return: (Float) Band gap in unit of ev.
    """

    xml_file = PATH + "/tmp/pwscf.xml"

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

    DMREF_mat_list, ini_idx = read_calculation_list(start = 200, end = 220)

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
    DFTU_all_elements = []#--Does it contain any DFTU_all elements
    DFTU_TM_elements = []#--Does it contain any DFTU_TM elements
    material_list = []

    for idx, mat in enumerate(DMREF_mat_list):

        #--initialize
        vc_finished = True
        dftu_all = True
        dftu_TM = True

        print("\n###")
        print("nubmer {} materials in the DMREF.csv data list is {}".format(idx+ini_idx,mat))
        material_list.append(mat)

        try:#Do we need DFTU for all calculations?
            whether_DFTU(mat, dftu_type = "All")
            DFTU_all_elements.append(True)
        except ValueError:
            print("No DFTU elements for DFTU type All for "+ mat + ".")
            dftu_all = False
            DFTU_all_elements.append(False)
            pass

        try:#Do we need DFTU for all calculations?
            whether_DFTU(mat, dftu_type = "TM_only")
            DFTU_TM_elements.append(True)
        except ValueError:
            print("No DFTU elements for DFTU type TM_only for "+ mat + ".")
            dftu_TM = False
            DFTU_TM_elements.append(False)
            pass

        try: #Parse vc-relax (DFT-gap), which should be available for all materials
            dft_bg = xml_parser(calc_root + "/" + mat + "/vc_relax")
            print("DFT band gap Parsed...")
            DFT_gap.append(dft_bg)
        except FileNotFoundError:
            print("VC-relax calculation failed!!")
            print(calc_root + "/" + mat + "/vc_relax")
            vc_finished = False
            DFT_gap.append(None), DFTU_all_gap.append(None), DFTU_TM_gap.append(None)
            continue#BE CAREFUL OVER HERE!!!

        if dftu_all: #vc_relax is finished using CONTINUE

            try:
                dftu_all_bg = xml_parser(calc_root + "/" + mat + "/final_scf_all")
                print("DFTU_all band gap Parsed...")
                DFTU_all_gap.append(dftu_all_bg)
            except FileNotFoundError:#--Failed
                print("DFTU_all calculation failed!!")
                DFTU_all_gap.append(None)
                pass
        else:
            DFTU_all_gap.append(dft_bg)#--this is dft band gap due to no U elements

        if dftu_TM:

            try:
                dftu_tm_bg = xml_parser(calc_root + "/" + mat + "/final_scf_TM")
                print("DFTU_TM band gap Parsed...")
                DFTU_TM_gap.append(dftu_tm_bg)
            except FileNotFoundError:
                print("DFTU_TM calculation failed!!")
                DFTU_TM_gap.append(None)
                pass
        else:
            DFTU_TM_gap.append(dft_bg)#--this is dft band gap due to no U elements

    #ZIP the results here
    print(len(material_list))
    print(len(DFT_gap))
    print(len(DFTU_all_elements))
    print(len(DFTU_all_gap))
    print(len(DFTU_TM_elements))
    print(len(DFTU_TM_gap))

    results = list(zip(material_list, DFT_gap, DFTU_all_elements, DFTU_all_gap,
                       DFTU_TM_elements, DFTU_TM_gap))

    cols = ["Materials", "DFT-gap", "DFTU_All_elements", "DFTU_All_bandgap",
            "DFTU_TM_elements", "DFTU_TM_bandgap"]

    df = pd.DataFrame(results, columns = cols)
    df.to_csv("Bandgap.csv")

if __name__=="__main__":

    main()
