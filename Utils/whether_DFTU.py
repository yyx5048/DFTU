import re
from Base.Elements import dftu_elements
from Base.Elements import transition_metal_elements

def whether_DFTU(chem_form, dftu_type = "All") -> str:
    """
    Verify if DFTU is needed based on the dftu_elements

    Args: (Str) chemical formula
    """
    if dftu_type == "All":
        dftu_atm = dftu_elements()

    elif dftu_type == "TM_only":
        dftu_atm = transition_metal_elements()

    else:
        raise ValueError("Invalid DFTU type")

    atm = re.findall('[A-Z][^A-Z]*', re.sub("\d+", "",re.sub(" ", "" ,chem_form)))

    if not any(x in dftu_atm for x in atm):
        raise ValueError("No DFTU elements in {}".format(chem_form))

    return
