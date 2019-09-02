def VC_settings():
    """
    Quantum Espresso basic settings for vc-relax calculations
    Args:
    Return: default input dictionary
    """
    pseudo_dir =('/gpfs/group/ixd4/default/Pseudos/GBRV_USPP_PBE_UPF_format')

    control = {"calculation":"vc-relax",
               "pseudo_dir" : pseudo_dir,
                 "verbosity" : 'high',
                 "restart_mode" : "from_scratch",
                 "wf_collect" : True,
                 "nstep" : 200,
                 "outdir" : "./tmp",
                 "max_seconds" : 172800}

    system = {"ecutwfc" : 90,
              "ecutrho" : 720,
              "occupations" : "smearing",
              "smearing" : "mv",
              "degauss" : 0.005}

    electrons = {"diagonalization" : "david",
                 "conv_thr" :1.0e-8,
                 "mixing_beta" : 0.50,
                 "electron_maxstep" : 250,
                 "mixing_mode" : "plain"}

    ions = {'ion_dynamics' : 'bfgs'}

    cell = {'cell_dynamics' : 'bfgs'}

    input_dict = {"CONTROL" : control,
                  "SYSTEM" : system,
                  "ELECTRONS" : electrons,
                  "IONS" : ions,
                  "CELL" : cell,
                  }

    return input_dict

def scfu():
    """
    Quantum Espresso basic settings for scf DFT+U calculations
    Args: (Int) number of bands (For fixed occupuations)
    Return: (Dict) default input dictionary
    """
    pseudo_dir =('/gpfs/group/ixd4/default/Pseudos/GBRV_USPP_PBE_UPF_format')

    control = {"calculation":"scf",
               "pseudo_dir" : pseudo_dir,
                 "verbosity" : 'high',
                 "restart_mode" : "from_scratch",
                 "wf_collect" : True,
                 "nstep" : 200,
                 "outdir" : "./tmp",
                 "max_seconds" : 172800}

    system = {"ecutwfc" : 90,
              "ecutrho" : 720,
              "occupations" : "fixed",
              "lda_plus_u": True,
              "lda_plus_u_kind":0,
              "U_projection_type" : 'ortho-atomic',
              "degauss" : 0.00}

    electrons = {"diagonalization" : "david",
                 "conv_thr" :1.0e-8,
                 "mixing_beta" : 0.50,
                 "electron_maxstep" : 250,
                 "mixing_mode" : "plain"}

    ions = {'!ion_dynamics' : 'bfgs'}

    cell = {'!cell_dynamics' : 'bfgs'}

    input_dict = {"CONTROL" : control,
                  "SYSTEM" : system,
                  "ELECTRONS" : electrons,
                  "IONS" : ions,
                  "CELL" : cell,
                  }

    return input_dict

def hp():
    """
    Quantum Espresso basic settings for HP calculations
    Args:
    Return:
    """

    with open('uscf.in','w') as f:

        f.write("""&inputUscf
prefix = 'pwscf',
outdir = '/tmp',
nq1 = 1 , nq2 = 1 , nq3 = 1,)
iverbosity = 2,
niter_ph = 150,
alpha_mix(1) = 0.1
find_atpert = 2,\n
/\n""")

    return
