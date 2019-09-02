import re
from pymatgen.core import structure
from Base.Elements import dftu_elements

def initialize_hubbard(structure):
    """
    Initialize the hubbard dictionary based on the existing DFTU elements in
    Quantum Espresso.
    Args: (Pymatgen structure)
    Return: hubbard_u dictionary with 1e-8 as initial value
    """
    u_list = []
    val = 1e-08
    dftu_list = dftu_elements()#--load DFTU list

    atom_list = list(structure.symbol_set)
    #--More robust, consider remove all the space, then digits, then capital letters
    #atom_list = re.findall('[A-Z][^A-Z]*', re.sub("\d+", "",re.sub(" ", "" ,structure.formula)))

    for site in atom_list:
        if site in dftu_list:
            u_list.append(site)

    val_list = [val]*len(u_list)

    return dict(zip(u_list,val_list))

def insert_hubbard_block(hubbard_u):
    """
    Insert the hubbard block (Hubbard_U(1) = 1.0e-8)

    Args: (Pymatgen structure)
    Return: hubbard_u dictionary with 1e-8 as initial value
    """
    hubbard_block = {}
    for idx, (k, v) in enumerate(hubbard_u.items()):
        _key = "Hubbard_U({})".format(idx+1)
        hubbard_block[_key] = v

    return hubbard_block

def reorder(structure, hubbard_u):
    """
    Create a copy of the structure but with elements in the order of hp.x, where
    all Hubbard atoms appear first in the atomic positions card.

    Args: structure: (Pymatgen structure) StructureData node
          hubbard_u: (Dict) a dictionary with the Hubbard U kinds and values

    Return: Reordered copy of structure base on the kind
    """
    reordered = structure.copy()
    reordered.clear()

    sites = structure.sites
    hubbard_kinds = list(hubbard_u.keys())
    #hubbard_kinds.sort(reverse=True)
    #--list of reversed hubbard atoms,, due to reverse pop in the while statement

    ordered_sites = []

    while hubbard_kinds:

        hubbard_kind = hubbard_kinds.pop()

        hubbard_sites = []
        remaining_sites = []

        condition = lambda s: s.species_string == hubbard_kind
        hubbard_sites = [s for s in sites if condition(s)]
        remaining_sites = [s for s in sites if not condition(s)]

        ordered_sites.extend(hubbard_sites)
        sites = remaining_sites

    # Extend the current site list with the remaining non-hubbard sites
    ordered_sites.extend(sites)

    for site in ordered_sites:
        reordered.append(site.species_string, site.frac_coords,
                         coords_are_cartesian = False)

    return reordered
