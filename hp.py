import os
import sys

from pymatgen.core import structure
from pymatgen.ext.matproj import MPRester

import subprocess
import pandas as pd

from Base.PBS_submit import pbs_submit
import Base.QE_input_settings as qe_input
from Utils.whether_DFTU import *
from Utils.Parser import pw_parser

#==============================================================================#
#  HP.x calculations                                                           #
#==============================================================================#
calc_root = os.getcwd()

def read_calculation_list(fname="./Data/Final_DMREF_Materials_List.csv",start = 201, end = 202):#200
    """
    Obtaion the list of structure to be calculated

    Args: (Str) Excel file name
          (Int) Start index
          (Int) End index
    Return: (numpy array) A list of chemical names (Str) read from DMREF list
    """
    mat_list = pd.read_csv(fname)
    return mat_list['Formula'].values[start:end], start

def mk_hp_folder(chem_form) -> str:
    """
    Making folder fo HP.X calculations, the validation of SCFU_first is
    performed.

    Args: Chemical formula for creating folder.
    Returns:
    """
    if os.path.isdir(chem_form):

        whether_DFTU(chem_form)#--Check do we need DFTU

        os.chdir(chem_form)
        os.chdir('first_scfu')

        res_dict = pw_parser()

        if res_dict['status'] == "DONE":
            os.chdir('../')
            os.mkdir('hp')
            os.chdir('hp')
            os.mkdir('tmp')
    else:
        raise FileNotFoundError("Folder did not created.")

    return res_dict

def copy_scfu_folder(src = "../first_scfu/*", dest = "./"):
    """
    Copy over the first_scfu folder for HP.X calculation. This enforce the
    first_scfu folder not touched, slower but guarantee the data integrity.
    """
    f = open("rysnc.sh",'w')
    f.write("""
        #!/bin/bash

        # SETUP OPTIONS
        export SRCDIR="%s"
        export DESTDIR="%s"
        export THREADS="4"

        # RSYNC DIRECTORY STRUCTURE
        rsync -zr -f"+ */" -f"- *" $SRCDIR/ $DESTDIR/ \
        # FIND ALL FILES AND PASS THEM TO MULTIPLE RSYNC PROCESSES
        cd $SRCDIR  &&  find . ! -type d -print0 | xargs -0 -n1 -P$THREADS -I%% rsync -az %% $DESTDIR/%%"""%(src,dest))
    f.close()
    subprocess.call(shlex.split("chmod +x rysnc.sh"))
    subprocess.call(shlex.split("sh rysnc.sh"))
    return

def main():

    #Read DMREF materials list CSV
    DMREF_mat_list = read_calculation_list()

    for idx, mat in enumerate(DMREF_mat_list):
        print("\n###")
        print("nubmer {} materials in the DMREF.csv data list is {}".format(idx,mat))

        try:#sequence of a hp.x operation

            mp_structure, mp_id = get_structure_from_mp(mat)

            mk_hp_folder(mat)#--create folder

            copy_scfu_folder()

            qe_input.hp()

            pbs_submit("hp",1)

            os.chdir(calc_root)

        except (ValueError, FileExistsError) as e:
            print(str(e))
            os.chdir(calc_root)
            continue
