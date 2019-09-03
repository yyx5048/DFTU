import os
import sys
import glob

import warnings
from pymatgen.core import structure
from pymatgen.ext.matproj import MPRester

from pymatgen.io.pwscf import PWInput

import numpy as np
import os
import shutil
import subprocess

def get_structure_from_mp(formula) -> str:
    """
    Obtain a crystal from the MP database via the API.

    Args: formula (str): A formula

    Returns: (Structure) the lowest energy structure on the Convex hull form MP
    database
    """
    m = MPRester("N7AIm1s2v43BQ6FT")
    entries = m.get_entries(formula, inc_structure="final")
    #-- returned the computed final structure (with vasp relaxation)

    if len(entries) == 0:
        raise ValueError("This crystal structure with formula {} has been\
                          removed from the current version of MP DB".format(formula))
    elif len(entries) > 1:
        warnings.warn("{} structures with formula {} found in MP, however only\
                       the lowest one is returned in this case".format(len(entries), formula))

    return min(entries, key = lambda e: e.energy_per_atom).structure


k_pts_grid = (2,2,2)
pseudo_dir ='/storage/home/qjc5019/Pseudopotentials/Pseudos/GBRV_USPP_Library/GBRV_USPP_PBE_UPF_format'

control = {"calculation":"scf",
	   "pseudo_dir" : pseudo_dir,
	   "verbosity" : 'high',
	   "restart_mode" : "from_scratch",
	   "wf_collect" : True,
	   "nstep" : 200,
	   "outdir" : "./tmp",
	   "!max_seconds" : 172800}



electrons = {"diagonalization" : "david",
	     "conv_thr" :1.0e-8,
	     "mixing_beta" : 0.50,
	     "electron_maxstep" : 250,
             "!mixing_mode" : "local-TF"}

ions = {'ion_dynamics' : 'bfgs'}

cell = {'cell_dynamics' : 'bfgs'}

elements_w_dftu = ['Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                   'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                   'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
                   'La', 'Sr','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',
                   'Ho','Er','Tm','Yb','Lu','Th','Pa','U', 'Np','Pu','Am',
                   'Cm','Bk','Cf','Es','Fm','Md','No','Lr','C', 'N', 'O', 'Al',
                   'Si', 'As' , 'Ga', 'In']


def main():

    #TODO pandas read in the calculations


    for structure_chem_formte_complete_Uscf_file(structure,dftu_elements_in_struc):
        qe_file = PWInput(structure,pseudo = pseudo_dict, control = control, system = system,
                          electrons=electrons,ions=ions,cell=cell,kpoints_grid = k_pts_grid,
                          kpoints_shift = (0,0,0))

        qe_file.write_file("dftu.in")

        shutil.copy('/gpfs/group/ixd4/default/tutorials/pymatgen_tutorials/restart_dftu.py', './')

        create_pbs_file()

        subprocess.call(['qsub','%s.dftu.pbs'%control['prefix'].replace('(','').replace(')','_')])
        os.chdir('./../') in chem_forms:

        system = {"ecutwfc": 90,
	   "ecutrho": 720,
	   "occupations": "fixed",
	   "!smearing" : "mv",
	   "!degauss" : 0.002,
	   "lda_plus_u": True,
           "lda_plus_u_kind":0,
           "U_projection_type" : 'ortho-atomic',
           "nosym" : True}

        if not os.path.isdir(structure_chem_form):
            os.mkdir(structure_chem_form)
        os.chdir(structure_chem_form)

        structure, formation_energy, band_gap = find_structure_bandgap(structure_chem_form)

        k = open('matproj_bandgap.txt','w')
        k.write("The materials project band gap for this material is: %s"%(band_gap))
        k.close()

        #oxidation_states = find_oxidations(structure)

        #oxidation_state = {"Zn":2,"O":-2}
        #print( oxidation_states)

        #for ox in oxidation_states:
        #structure.add_oxidation_state_by_element(ox)

        system['nbnd'] = len(structure.sites)*8
        control["prefix"] = structure_chem_form

        dftu_elements_in_struc = []
        #print ( structure_chem_form)
        at_no = 1
        for el in sorted(structure.symbol_set):
            if el in elements_w_dftu:
                #print (at_no)
                system['Hubbard_U(%s)'%(at_no)] = 1e-8
                if el not in dftu_elements_in_struc:
                    dftu_elements_in_struc.append(el)
            at_no +=1



        create_complete_Uscf_file(structure,dftu_elements_in_struc)
        qe_file = PWInput(structure,pseudo_lib = 'GBRV_US',control = control,system = system,electrons=electrons,ions=ions,cell=cell,kpoints_grid = k_pts_grid, kpoints_shift = (0,0,0))
        qe_file.write_file("scf.in")
        #create_environ_file(slab)
        shutil.copy('/gpfs/group/ixd4/default/tutorials/pymatgen_tutorials/restart_dftu.py', './')
        create_pbs_file()
        subprocess.call(['qsub','%s.dftu.pbs'%control['prefix'].replace('(','').replace(')','_')])
        os.chdir('./../')

    return

def create_complete_Uscf_file(structure,dftu_elements_in_struc):
    f = open('uscf.in','w')

    ### I found that skipping atoms and doing an entire calculation at once works better
    ### than splitting up the atoms. You should check out the HUBBARD docs tho

    f.write("""&inputUscf
   prefix = '%s',
   outdir = '%s',
   nq1 = 1 , nq2 = 1 , nq3 = 1,\n"""%(control['prefix'],control['outdir']))
    f.write("""   iverbosity = 2,
   niter_ph = 150,
   search_method = 2,\n""")

    f.write('/\n')

    f.close()

    return

def create_pbs_file():
    f = open('%s.dftu.pbs'%control['prefix'].replace('(','').replace(')','_'),'w')

    f.write("""#PBS -l nodes=1:ppn=20
               #PBS -l walltime=15:00:00:00
               #PBS -l pmem=10gb
               #PBS -j oe
               #PBS -A ixd4_d_g_lc_default

               cd $PBS_O_WORKDIR


               echo " "
               echo "Job started on `hostname` at `date`"
               echo " "

               mpirun -np 20 /gpfs/group/ixd4/default/software/espresso-5.0.2_Uscf_26.07.2017/bin/pw.x -in scf.in > scf.out

               mpirun -np 20 /gpfs/group/ixd4/default/software/espresso-5.0.2_Uscf_26.07.2017/bin/Uscf.x -in uscf.in > uscf.out

               /gpfs/group/ixd4/default/software/anaconda3/bin/python restart_dftu.py

               echo " "
               echo "Job Ended at `date`"
            """)
    f.close()

    return

def create_Uscf_file(atom_no):
    f = open('uscf.at{}.in'.format(atom_no),'w')

    f.write("""&inputUscf
               prefix = {},
               outdir = {}',
               nq1 = 1, nq2 = 1, nq3 = 1,
               do_one_only({}) = .true.,
               iverbosity = 2,
               niter_ph = 150,
               /""".format(control['prefix'],control['outdir'],atom_no))
    return



def find_structure_bandgap(chemical_formula):
    API_ID = "N7AIm1s2v43BQ6FT"
    from pymatgen.ext.matproj import MPRester
    # Gets the structure from the Materials Project Interface
    # I am using my api id for MPRester. Others should get their own API ID from
    # the materials project dashboard
    with MPRester(API_ID) as m:
        #Generates a list of all structures with the chemical formula given noting formation energy
        data = m.query(chemical_formula,['material_id','formation_energy_per_atom','band_gap'])
        #Selecting the structure with the lowest formation energy (i.e. most stable)
        material_id = ''
        formation_energy = 0.0
        band_gap = 0.0
        for entry in data:
            if entry['formation_energy_per_atom'] < formation_energy:
                formation_energy = entry['formation_energy_per_atom']
                material_id = entry['material_id']
                band_gap = entry['band_gap']
        if material_id == '':
            print("There is no structure with a formation energy below 0.0 eV")
            return 0.0, 0.0, 0.0
        #Downloads the structure with lowest formation energy from materials project website
        #print(material_id)
        structure = m.get_structure_by_material_id(material_id)#,conventional_unit_cell = True)
        return structure, formation_energy, band_gap

if __name__=='__main__':
    main()
