from pymatgen.io.ase import AseAtomsAdaptor

def perturb_select_atoms():
    """
    Modify the hp.x input file by adding perturb_only_atom(i).

    Args:
    Returns: (dftu.in) updated with perturb_only_atom(i).
    """
    perturb_idx = []
    perturb_atm = []
    ase_s = read("../first_scfu/dftu.in", format = "espresso-in")

    pymat_s = AseAtomsAdaptor.get_structure(ase_s)
    for idx, site in enumerate(pymat_s):
        if site.specie.symbol in dftu_elements() and not in perturb_atm:
            pertub_atm.append(site.specie.symbol)
            perturb_idx.append(idx)
    #--remove the last line and adding back the
    with open("dftu.in","r+") as foo:
        foo.seek(0,os.SEEK_END)
        pos = fil.tell() - 1
        while pos > 0 and foo.read(1) != "\n":
            pos -= 1
            foo.seek(pos, os.SEEK_SET)
        if pos > 0:
            foo.seek(pos, os.SEEK_SET)
            foo.truncate()

        #TODO: check the writing format.
        for k in perturb_idx:
            f.write("perturb_only_atom({}) = .true.,\n".format(k))
        f.write("/\n")
    return
