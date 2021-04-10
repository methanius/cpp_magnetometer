#!/bin/bash
#SBATCH --job-name=SPARSE_PQS_OU
#SBATCH --partition=q40,q36,q28,q24
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=30:00:00


echo "========= Job started  at `date` =========="
VARIABLES=('10000' '50' '0.04' '1000')

# Go to the directory where this job was submitted
cd $SLURM_SUBMIT_DIR

# copy inputdata and the executable to the scratch-directory
cp autocorrelation_OU.py /scratch/$SLURM_JOB_ID
cp -r ../sparse_OU_library/ /scratch/$SLURM_JOB_ID
cp ../ehrenfest_chain.py /scratch/$SLURM_JOB_ID

# change directory to the local scratch-directory:
cd /scratch/$SLURM_JOB_ID
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
python3 autocorrelation_OU.py ${VARIABLES[@]} 

# copy home the outputdata:
cp *pkl /home/normann/programs/autocorrelation_OU/data/

echo "========= Job finished at `date` =========="
#
