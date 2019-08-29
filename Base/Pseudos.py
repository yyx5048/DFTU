import os
import sys

import warnings
from pymatgen.core import structure
from pymatgen.ext.matproj import MPRester

def Pseudos(structure) -> structure:
    """
    GBRV Pseudopotentials

    Args: (pymatgen.core.structure.Structure) pymatgen structure

    Return: (dictionary) pseudopotentials
    """
    pseudo = {"Ba": "ba_pbe_v1.uspp.F.UPF", "Pb":"pb_pbe_v1.uspp.F.UPF",
              "O": "o_pbe_v1.2.uspp.F.UPF", "Ag": "ag_pbe_v1.uspp.F.UPF",
              "Fe": "fe_pbe_v1.2.uspp.F.UPF", "Na": "na_pbe_v1.uspp.F.UPF",
              "Si": "si_pbe_v1.uspp.F.UPF", "Al": "al_pbe_v1.uspp.F.UPF",
              "F": "f_pbe_v1.2.uspp.F.UPF", "Nb":" nb_pbe_v1.uspp.F.UPF",
              "Sn": "sn_pbe_v1.uspp.F.UPF", "As" : "as_pbe_v1.uspp.F.UPF",
              "Ga": "ga_pbe_v1.uspp.F.UPF", "Ni":"ni_pbe_v1.2.uspp.F.UPF",
              "S": "s_pbe_v1.2.uspp.F.UPF", "Au" : "au_pbe_v1.uspp.F.UPF",
              "Ge": "ge_pbe_v1.uspp.F.UPF", "N": "n_pbe_v1.2.uspp.F.UPF",
              "Sr": "sr_pbe_v1.uspp.F.UPF", "Ba" : "ba_pbe_v1.uspp.F.UPF",
              "O": "o_pbe_v1.2.uspp.F.UPF",
              "Ta" : "ta_pbe_v1.uspp.F.UPF", "Be" : "be_pbe_v1.2.uspp.F.UPF",
              "Hf" : "hf_pbe_v1.uspp.F.UPF", "Os": "os_pbe_v1.2.uspp.F.UPF",
              "Tc": "tc_pbe_v1.uspp.F.UPF" , "Bi": "bi_pbe_v1.uspp.F.UPF",
              "Hg" : "hg_pbe_v1.uspp.F.UPF", "Pb" : "pb_pbe_v1.uspp.F.UPF",
              "Te" : "te_pbe_v1.uspp.F.UPF", "B": "b_pbe_v1.01.uspp.F.UPF",
              "H" : "h_pbe_v1.uspp.F.UPF", "Pd": "pd_pbe_v1.2.uspp.F.UPF",
              "Ti" : "ti_pbe_v1.uspp.F.UPF", "Br" : "br_pbe_v1.uspp.F.UPF",
              "In" : "in_pbe_v1.uspp.F.UPF", "P" : "p_pbe_v1.uspp.F.UPF",
              "Tl" : "tl_pbe_v1.2.uspp.F.UPF", "Ca" : "ca_pbe_v1.uspp.F.UPF",
              "I" : "i_pbe_v1.uspp.F.UPF", "Pt" : "pt_pbe_v1.uspp.F.UPF",
              "V" : "v_pbe_v1.uspp.F.UPF", "Cd" : "cd_pbe_v1.uspp.F.UPF",
              "Ir" : "ir_pbe_v1.2.uspp.F.UPF", "Rb": "rb_pbe_v1.uspp.F.UPF",
              "W": "w_pbe_v1.2.uspp.F.UPF", "Cl" : "cl_pbe_v1.2.uspp.F.UPF",
              "K" : "k_pbe_v1.uspp.F.UPF", "Re" : "re_pbe_v1.2.uspp.F.UPF",
              "Y" : "y_pbe_v1.uspp.F.UPF", "Co" : "co_pbe_v1.2.uspp.F.UPF",
              "La" : "la_pbe_v1.uspp.F.UPF", "Rh" : "rh_pbe_v1.2.uspp.F.UPF",
              "Zn" : "zn_pbe_v1.uspp.F.UPF", "C" : "c_pbe_v1.2.uspp.F.UPF",
              "Li" : "li_pbe_v1.uspp.F.UPF", "Ru" : "ru_pbe_v1.2.uspp.F.UPF",
              "Zr" : "zr_pbe_v1.uspp.F.UPF", "Cr" : "cr_pbe_v1.2.uspp.F.UPF",
              "Mg" : "mg_pbe_v1.uspp.F.UPF", "Sb" : "sb_pbe_v1.uspp.F.UPF",
              "Cs" : "cs_pbe_v1.uspp.F.UPF", "Mn" : "mn_pbe_v1.2.uspp.F.UPF",
              "Sc" : "sc_pbe_v1.uspp.F.UPF", "Cu" : "cu_pbe_v1.2.uspp.F.UPF",
              "Mo" : "mo_pbe_v1.uspp.F.UPF", "Se" : "se_pbe_v1.uspp.F.UPF"}

    pseudo_atm = list(structure.symbol_set)

    pseudo_name = [pseudo[k] for k in pseudo_atm]

    pseudo_dict = dict(zip(pseudo_atm,pseudo_name))

    return pseudo_dict
