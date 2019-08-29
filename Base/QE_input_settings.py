def VC_settings():
    """
    Quantum Espresso basic settings for DFT+U calculations
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
