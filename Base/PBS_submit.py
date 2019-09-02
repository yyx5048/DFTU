import subprocess

def create_vc_pbs_file(nodes):
    """
    Generate PBS file for submission

    Args: (Int) nodes: number of nodes, using default 20 ppn.
    Return:
    """
    with open('calculation.pbs','w') as f:

        f.write("""#PBS -l nodes={}:ppn=20
#PBS -l walltime=48:00:00
#PBS -l pmem=10gb
#PBS -j oe
#PBS -A ixd4_c_g_sc_default

cd $PBS_O_WORKDIR

module purge
module use /gpfs/group/dml129/default/sw/modules
module load quantum_espresso/6.4.1

echo " "
echo "Job started on `hostname` at `date`"
echo " "

mpirun -np {} pw.x -ndiag 1 -in dftu.in > dftu.out

echo " "
echo "Job Ended at `date`"
            """.format(nodes,nodes*20))

    return

def create_hp_pbs_file(nodes):
    """
    Generate PBS file for submission

    Args: (Int) nodes: number of nodes, using default 20 ppn.
    Return:
    """
    with open('calculation.pbs','w') as f:

        f.write("""#PBS -l nodes={}:ppn=20
#PBS -l walltime=48:00:00
#PBS -l pmem=10gb
#PBS -j oe
#PBS -A ixd4_c_g_sc_default

cd $PBS_O_WORKDIR

module purge
module use /gpfs/group/dml129/default/sw/modules
module load quantum_espresso/6.4.1

echo " "
echo "Job started on `hostname` at `date`"
echo " "

mpirun -np {} hp.x -in dftu.in > dftu.out

echo " "
echo "Job Ended at `date`"
            """.format(nodes,nodes*20))

def pbs_submit(typ, nodes):
    """
    Submit PBS file

    Args: (Str) typ: calculation types, "pw" or "hp"
          (Int) nodes: number of nodes, using default 20 ppn.
    Return:
    """
    if typ == "pw":
        create_vc_pbs_file(nodes)
        subprocess.call(['qsub','calculation.pbs'])
    elif typ == "hp":
        create_hp_pbs_file(nodes)
        subprocess.call(['qsub','calculation.pbs'])

    else:
        raise ValueError("Invalid calculation type passed")

    return
