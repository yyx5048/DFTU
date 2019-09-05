'''
The purpose of this script is to:
-Retreive chemical info from Materials Project Page
-Convert Info into Quantum Espresso input file & submit to server
-Check final output to ensure that it's converged
-Calculate Hubbard U param.
-Iterate until Hubbard U converges w/in .01 eV
-Extract material properties
'''


import sys, os, subprocess, glob
import numpy as np
import xml.etree.ElementTree as ET
from ase.io import read, write


# Global variables
elements_w_dftu = ['Sc', 'Ti',  'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                   'La', 'Hf', 'Ta',  'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
                   'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
                   'Th', 'Pa',  'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']

def main():

    # Grab all the files in the current folder
    files = glob.glob("*")

    # Check the status of the calculation
    if "Status.txt" in files:
        # Check which file is currently going
        file_status = Check_Status()
        try:
            iter_num = int(file_status.split("_")[-1])
        except:
            iter_num = file_status.split("_")[-1]

        # Parse input file for pertinent information
        if "VCR" in file_status:
            input_dict = Parse_Input('VCRelax',iter_num)
        elif "HP" in file_status:
            input_dict = Parse_Input('SCF',iter_num)
        elif "SCF" in file_status:
            input_dict = Parse_Input("SCF","Final")

        # Check job status
        Job_Done, input_dict = Error_Check(file_status,input_dict)

        # Chemical formula name
        chem_formula = glob.glob("*_VCRelax_1.in")[0].split("_")[0]

        if Job_Done == True:
            # Check to see if the previous calculation was a VCRelax or HP
            if "VCRelax" in file_status:
                # Parse the XML file for pertinent information
                Parse_XML(file_status.split("_")[0],input_dict)
                # Write the HP file feeding in the correct information
                input_dict['CONTROL']['outdir'] = './HP_tmp'
                Write_HP(chem_formula,iter_num,input_dict)
                # Extract data from the VCRelax
                data = read("%s_%s.out"%(chem_formula,file_status),format='espresso-out')
                # Modify the atoms that need a Hubbard U
                # Check if iter_num equals 1 otherwise the Hubbard U should be set from the read in values.
                if iter_num == 1:
                    input_dict = Set_Hub_U(data.get_chemical_symbols(),input_dict)
                # Grab the list of pseudopotentials
                pseudos = Pseudos(data.get_chemical_symbols())
                # Write the SCF file for the Hubbard U
                write("%s_SCF_%s.in"%(chem_formula,iter_num),data,format='espresso-in',pseudopotentials=pseudos,input_data=input_dict,kspacing=0.04)
                # Write the pbs file
                Submit_pbs_in(chem_formula,'U',iter_num,input_dict)
                # Update the status file
                lines = open('Status.txt','r').readlines()
                for i in range(len(lines)):
                    if file_status in lines[i]:
                        lines[i] = "%s Done\n"%(file_status)
                lines.append("HP_%s Running"%(file_status.split("_")[-1]))
                temp = open('Status.txt','w')
                for item in lines:
                    temp.write(item)
                temp.close()

            elif "HP" in file_status:
                # Extract data from the VCRelax
                data = read("%s_%s_%s.out"%(chem_formula,"VCRelax",iter_num),format='espresso-out')
                # Increase the iteration number
                iter_num += 1
                # Read in the Hubbard U parameters from data file and update input dictionary.
                input_dict = Read_Hubbard_U(input_dict)
                input_dict['CONTROL']['outdir'] = './HP_tmp'
                input_dict['CONTROL']['calculation'] = 'scf'
                # Grab the list of pseudopotentials
                pseudos = Pseudos(data.get_chemical_symbols())
                # Write the SCF file for the Hubbard U
                write("%s_SCF_Final.in"%(chem_formula),data,format='espresso-in',pseudopotentials=pseudos,input_data=input_dict,kspacing=0.04)
                # Write the pbs file
                Submit_pbs_in(chem_formula,'SCF',iter_num,input_dict)
                # Update the status file
                lines = open('Status.txt','r').readlines()
                for i in range(len(lines)):
                    if file_status in lines[i]:
                        lines[i] = "%s Done\n"%(file_status)
                lines.append("SCF_Final Running")
                temp = open('Status.txt','w')
                for item in lines:
                    temp.write(item)
                temp.close()

            elif "SCF" in file_status:
                # Parse the XML file for pertinent information
                Parse_XML(file_status.split("_")[0],input_dict)
                # Update the status file
                lines = open('Status.txt','r').readlines()
                for i in range(len(lines)):
                    if file_status in lines[i]:
                        lines[i] = "%s Done\n"%(file_status)
                temp = open('Status.txt','w')
                for item in lines:
                    temp.write(item)
                temp.close()
                print("Final SCF calculation finished correctly.")

        elif Job_Done == 'Restart' and "VCR" in file_status:
            # Extract data from the VCRelax
            data = read("%s_%s.out"%(chem_formula,file_status),format='espresso-out')
            # Grab the list of pseudopotentials
            pseudos = Pseudos(data.get_chemical_symbols())
            # Write the input file
            write("%s_%s.in"%(chem_formula,file_status),data,format='espresso-in',pseudopotentials=pseudos,input_data=input_dict,kspacing=0.04)
            # Write the pbs file
            Submit_pbs_in(chem_formula,'VCR',int(file_status.split("_")[-1]),input_dict)

        elif Job_Done == 'Restart' and "HP" in file_status:
            # Extract data from the SCF
            data = read("%s_SCF_1.in"%(chem_formula),format='espresso-in')
            input_dict['CONTROL']['calculation'] = 'scf'
            input_dict['CONTROL']['temp_dir'] = './HP_tmp'
            # Grab the list of pseudopotentials
            pseudos = Pseudos(data.get_chemical_symbols())
            # Write the input file
            write("%s_SCF_1.in"%(chem_formula),data,format='espresso-in',pseudopotentials=pseudos,input_data=input_dict,kspacing=0.04)
            # Write the pbs file
            Submit_pbs_in(chem_formula,'U',int(file_status.split("_")[-1]),input_dict)

        else:
            print("Unknown error has occurred. You will need to check the files.")
            sys.exit()

    else:
        Initialize_Structure()

def Initialize_Structure():

    from pymatgen.io.ase import AseAtomsAdaptor

    #pass in structure as first variable in command line
    try:
        chem_formula = sys.argv[1]
    except:
        print("Please type the chemical formula after the python script.")
        chem_formula = input()

    if not os.path.isdir(chem_formula):
        os.mkdir(chem_formula)
    os.chdir(chem_formula)

    #get chem_formula from materials project page and returns structure info
    structure, e_form, bandgap = Mat_Proj_Struct(chem_formula)

    #convert to from pymatgen to ase
    ase_data = AseAtomsAdaptor()
    ase_data = ase_data.get_atoms(structure)

    # Order the atoms with the transition metals first
    U_atoms = []
    No_U_atoms = []
    temp = []
    # Creat a list of chemical symbols and indices
    for i,j in enumerate(ase_data.get_chemical_symbols()):
        temp.append([j,i])
    # Sort the atoms alphabetically by the atom type
    temp.sort()
    for atom in temp:
        if atom[0] in elements_w_dftu:
            U_atoms.append(atom[1])
        else:
            No_U_atoms.append(atom[1])
    # Reorder the ase Atom type by the atoms with U and without U.
    ase_data = ase_data[U_atoms+No_U_atoms]

    pseudos = Pseudos(ase_data.get_chemical_symbols())
    QE_input = QE_Settings()
    QE_input['SYSTEM']['nat'] = len(ase_data.get_chemical_symbols())

    write("%s_VCRelax_1.in" %(chem_formula), ase_data, format = "espresso-in", pseudopotentials=pseudos, input_data = QE_input, kspacing=0.04)
    Submit_pbs_in(chem_formula, "VCR", 1,QE_input)

    f = open("Status.txt",'w')
    f.write("VCRelax_1 Running")
    f.close()

    os.chdir('../')

def Check_Status():

    # Open the status file and read which calculation is currently going
    with open('Status.txt','r') as infile:
        for line in infile:
            temp = line.split()
            if temp[-1] == "Done":
                continue
            elif temp[-1] == "Running":
                return temp[0]


def Parse_Input(filetype='VCRelax',iternum=1):

    input_file = glob.glob("*_%s_%s.in"%(filetype,iternum))[0]

    input_dict = {}

    parse = False

    with open(input_file,'r') as infile:
        for line in infile:
            if "&" in line and not parse:
                key = line.strip().lstrip("&")
                input_dict[key] = {}
                parse = True

            elif "/\n" in line:
                parse = False

            elif parse:
                temp = line.strip().rstrip(",").split("=")
                if 'd-' in temp[1]:
                    temp[1] = temp[1].replace('d','e')
                elif '.true.' in temp[1]:
                    temp[1] = True

                if temp[1] == True:
                    input_dict[key][str(temp[0].strip())] = temp[1]
                else:
                    try:
                        input_dict[key][str(temp[0].strip())] = int(temp[1])
                    except:
                        try:
                            input_dict[key][str(temp[0].strip())] = float(temp[1])
                        except:
                            input_dict[key][str(temp[0].strip())] = str(temp[1].strip().strip("'"))

    return input_dict

def Pseudos(atype):
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
              "Hf" : "hf_pbe_plus4_v1.uspp.F.UPF", "O": "o_pbe_v1.2.uspp.F.UPF",
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

    pseudo_out = {}
    #narrow down pseudo list
    for atom in set(atype):
        pseudo_out[atom] = pseudo[atom]

    return pseudo_out


def QE_Settings():

    QE_Input = {"CONTROL" : {"calculation":"vc-relax",
                             "pseudo_dir": "/gpfs/group/ixd4/default/Pseudos/GBRV_USPP_PBE_UPF_format",
                             "verbosity" : "high",
                             "wf_collect" : True,
                             "nstep" : 200,
                             "outdir" : "./VCR_tmp",
                             "prefix" : "USCF",
                             "etot_conv_thr" : 1e-5,
                             "forc_conv_thr" : 1e-4 },
                "SYSTEM" : {"ecutwfc" : 90,
                            "ecutrho" : 720,
                            "occupations" : "smearing",
                            "smearing" : "mv",
                            "degauss" : 0.005 },
                "ELECTRONS" : {"conv_thr" : 1.0e-6,
                               "mixing_beta" : 0.70,
                               "electron_maxstep" : 100,
                               "mixing_mode" : "local-TF"},
                "IONS" : {"!ion_dynamics" : "damp",
                          "upscale"       : 10000},
                "CELL" : {"!cell_dofree" : "xy"} }
    return QE_Input

def Mat_Proj_Struct(chemical_formula):

    from pymatgen.ext.matproj import MPRester

    #ID num for accessing Materials Project (Julian's)
    API_ID = "7ThaNCc9ib1oRHz9DlI"

    # Gets the structure from the Materials Project Interface
    with MPRester(API_ID) as m:
        #Generates a list of all structures with the chemical formula given noting formation energy
        data = m.query(chemical_formula,['material_id','formation_energy_per_atom','band_gap'])
        #Selecting the structure with the lowest formation energy (i.e. most stable)
        material_id = ''
        formation_energy = 0.0
        band_gap = 0.0
        for entry in data:
            if entry['formation_energy_per_atom'] <= formation_energy:
                formation_energy = entry['formation_energy_per_atom']
                material_id = entry['material_id']
                band_gap = entry['band_gap']
        if material_id == '':
            print("There is no structure with a formation energy below 0.0 eV")
            return 0.0, 0.0, 0.0
        #Downloads the structure with lowest formation energy from materials project website
        #print(material_id)
        structure = m.get_structure_by_material_id(material_id)#,conventional_unit_cell = True)

    # For Nicole
    #ASCII_Bird()

    f = open("Bandgap.txt",'a')
    f.write("""Materials Project
Material ID: %s
Bandgap: %s
"""%(material_id,band_gap))
    f.close()
    return structure, formation_energy, band_gap

#write and submit a PBS file to the query
def Submit_pbs_in(chem_formula, calc_type, iternum, input_dict, n = 1, ppn = 20, walltime = 48, queue = 'open'):

    # Specify filename based on calculation type
    filename = ''
    if calc_type == "VCR":
        filename = chem_formula + "_VCRelax_" + str(iternum)
    elif calc_type == "U":
        SCF = chem_formula + "_SCF_" + str(iternum)
        filename = chem_formula + "_HP_" + str(iternum)
    elif calc_type == "SCF":
        filename = chem_formula + "_SCF_Final"

    # Change requested resources based on size of system
    if input_dict['SYSTEM']['nat'] > 25:
        walltime = 96
        queue = 'ixd4_e_g_bc_default'

    f=open("%s.pbs" %(filename), 'w')
    f.write("""#!/bin/bash
#PBS -l nodes={0}:ppn={1}
#PBS -l pmem=5gb
#PBS -l walltime={2}:00:00
#PBS -j oe
#PBS -A {3}

module purge
module use /gpfs/group/dml129/default/sw/modules
module load intel/2018
module load python

cd $PBS_O_WORKDIR

""".format(n,ppn,walltime,queue))

    # Line to append for VCRelax calculations
    if calc_type == "VCR" or calc_type == "SCF":
        f.write("""mpirun -np 20 /gpfs/group/ixd4/default/software/qe_6.4.1-environ_1.1/qe_6.4.1/bin/pw.x -ndiag 1 -in {0}.in > {0}.out

""".format(filename))

    # Lines to append for a Hubbard U calculation
    elif calc_type == 'U':
        f.write("""mpirun -np 20 /gpfs/group/ixd4/default/software/qe_6.4.1-environ_1.1/qe_6.4.1/bin/pw.x -ndiag 1 -in {0}.in > {0}.out

mpirun -np 20 /gpfs/group/ixd4/default/software/qe_6.4.1-environ_1.1/qe_6.4.1/bin/hp.x -ndiag 1 -in {1}.in > {1}.out

""".format(SCF,filename))

    f.write("python ~/work/Codes/bandgap_U.py >> python_output.txt")

    f.close()

    subprocess.call("qsub %s.pbs"%(filename), shell = True)

#check to see if file finished correctly, return false if changes need to be made
def Error_Check(file_status,input_dict):

    filename = glob.glob("*_%s.out"%(file_status))[0]

    Job_Done = False

    if "VCRelax" in filename:
        with open(filename) as infile:
            for line in infile:

                # Electron density did not converge with fixed occupations
                # Restart simulation with smearing

                # Check by hand if maximum iterations were reached with smearing turned on
                if "Iteration #%3d"%(input_dict["ELECTRONS"]["electron_maxstep"]) in line or "convergence NOT acheived after" in line and input_dict['SYSTEM']['occupations'] == 'smearing':
                    print("Maximum iterations reached with smearing turned on. Check by hand.")
                    sys.exit()

                elif "Iteration #%3d"%(input_dict["ELECTRONS"]["electron_maxstep"]) in line or "convergence NOT acheived after" in line:
                    print("Reached Maximum electron iterations")
                    input_dict['CONTROL']['calculation'] = 'scf'
                    input_dict['SYSTEM']['occupations'] = 'smearing'
                    input_dict['SYSTEM']['smearing'] = 'mv'
                    input_dict['SYSTEM']['degauss'] = 0.005
                    Job_Done = 'Restart'

                elif 'iterations completed, stopping' in line:
                    print("Maximum number of iterations reached in Wentzcovitch Damped Dynamics.")
                    sys.exit()

                elif "The maximum number of steps has been reached." in line:
                    print("Maximum number of ionic/electronic relaxation has been reached.")
                    sys.exit()

                # Take the final coordinates and restart the calculation
                elif "Maximum CPU time exceeded" in line:
                    print("Maximum CPU time exceeded")
                    Job_Done = "Restart"

                # Will come back and fill in with possible errors as we get them.
                elif "%%%%%%%%%%" in line:
                    print("Fatal error %%%%% \n Manual Edit \n")
                    sys.exit()

                elif "JOB DONE." in line:
                    print("NO ERRORS")
                    Job_Done = True

    elif 'HP' in filename:
        files = glob.glob("*")
        with open(filename) as infile:
            for line in infile:

                if "%%%%%" in line:
                    print("Simulation found an error. Check by hand.")
                    sys.exit()

                if "WARNING: The Fermi energy shift is too big!" in line:
                    input_dict['SYSTEM']['occupations'] = 'fixed'
                    del input_dict['SYSTEM']['degauss']
                    del input_dict['SYSTEM']['smearing']
                    Job_Done = 'Restart'

                elif "JOB DONE." in line and '%s.Hubbard_parameters.dat'%(input_dict['CONTROL']['prefix']) in files:
                    print("NO ERRORS")
                    Job_Done = True

    elif 'SCF' in filename:
        with open(filename) as infile:
            for line in infile:
                if "JOB DONE." in line:
                    Job_Done = True

    return Job_Done, input_dict

def Write_HP(chem_formula,iter_num,QE_Input):

    f = open("%s_HP_%s.in"%(chem_formula,iter_num), 'w')
    #what is &inputhp
    f.write("""&inputhp
prefix = '%s',
outdir = '%s',
nq1 = 2, nq2 = 2, nq3 = 2,
iverbosity = 2
alpha_mix(1) = 0.05
/"""%(QE_Input["CONTROL"]["prefix"], QE_Input["CONTROL"]['outdir']))
    f.close()

def Set_Hub_U(chem_sym,QE_Inputs):

    #allow hub_U calculations
    QE_Inputs["SYSTEM"]["lda_plus_U"] = True

    unique_atoms = []
    for atom in chem_sym:
        if atom not in unique_atoms:
            unique_atoms.append(atom)

    # Check to see if any atoms are in the Hubbard U list
    check = False
    for atom in unique_atoms:
        if atom in elements_w_dftu:
            check = True
    if not check:
        print("No atoms in this system need a Hubbard U.")
        sys.exit()

    for i, atom in enumerate(unique_atoms):
        if atom in elements_w_dftu: #NAME OF DICTIONARY W TRANSITION METALS
            QE_Inputs["SYSTEM"]["Hubbard_U(%s)"%(i+1)]=1e-8 #?

    QE_Inputs['CONTROL']['calculation'] = 'scf'

    return QE_Inputs

def Read_Hubbard_U(input_dict):

    Hubbard_dict = {}
    parse = False

    with open('%s.Hubbard_parameters.dat'%(input_dict['CONTROL']['prefix']),'r') as infile:
        for line in infile:
            if "Hubbard U (eV)" in line:
                parse = True
            elif line in ['\n','\r\n'] or '--------' in line:
                parse = False
            elif parse:
                Hub_data = line.split()
                if Hub_data[2] not in Hubbard_dict.keys():
                    Hubbard_dict[Hub_data[2]] = [float(Hub_data[-1])]
                else:
                    Hubbard_dict[Hub_data[2]].append(float(Hub_data[-1]))
            else:
                continue

    keys = Hubbard_dict.keys()

    input_dict['SYSTEM']['lda_plus_u'] = True

    for i,key in enumerate(keys):
        print(Hubbard_dict[key])
        input_dict['SYSTEM']['Hubbard_U(%s)'%(i+1)] = np.round(np.average(Hubbard_dict[key]),decimals=2)

    return input_dict

def Parse_XML(calc_type,input_dict):

    # Get the file path to the xml file
    xml_file = "%s/%s.xml"%(input_dict['CONTROL']['outdir'],input_dict['CONTROL']['prefix'])

    # Get the xml information
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Initialize large values for CBM and VBM
    CBM = 100000
    VBM = -100000

    # Iterate through each ks_energies for the eigenvalues and occupation state
    for ks in root.iter('ks_energies'):
        for eig in ks.iter('eigenvalues'):
            eigen = np.array([s for s in eig.text.split()]).astype('float')
        for occ in ks.iter('occupations'):
            occup = np.array([s for s in occ.text.split()]).astype('float')

        # Parse the eigenvalues and occupations to compare the CBM and VBM values from previous kpoints
        # When the occupation reaches less than 10^-6 consider that state to be empty.
        for j in range(len(occup)):
            if occup[j+1] <= 1e-4:
                CBMt = eigen[j+1]
                VBMt = eigen[j]
                break

        # Update the CBM and VBM values
        if CBMt < CBM:
            CBM = CBMt
        if VBMt > VBM:
            VBM = VBMt

    # Convert values to eV
    CBM = CBM*27.211396132
    VBM = VBM*27.211396132

    # Print the information to a text file
    f = open("Bandgap.txt",'a')
    f.write("""%s
CBM: %s eV
VBM: %s eV
Bandgap: %s eV
"""%(calc_type,CBM,VBM,(CBM-VBM)))
    f.close()

def ASCII_Bird():

    f = open("Bandgap.txt",'a')
    f.write(r"""                             _..._
                            \_.._ `'-.,--,
                             '-._'-.  `\a\\
                                 '. `_.' (|
                                   `7    ||
                                   /   .' |
                                  /_.-'  ,J
                                 /         \
                                ||   /      ;
                     _..        ||  |       |  /`\.-.
                   .' _ `\      `\  \       |  \_/__/
                  /  /e)-,\       '. \      /.-` .'\
                 /  |  ,_ |        /\ `;_.-'_.-'`\_/
                /   '-(-.)/        \_;(((_.-;
              .'--.   \  `       .(((_,;`'.  \
             /    `\   |   _.--'`__.'  `\  '-;\
           /`       |  /.-'  .--'        '._.'\\
         .'        ;  /__.-'`             |  \ |
       .'`-'_     /_.')))                  \_\,_/
      / -'_.'---;`'-)))
     (__.'/   /` .'`
      (_.'/ /` /`
        _|.' /`
       ` __.'|
      .-'||  |
         \_`/

""")
    f.close()

if __name__ == "__main__":
    main()
