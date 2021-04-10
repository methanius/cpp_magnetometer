#!/bin/bash
#SBATCH --job-name=SPARSE_PQS_OU
#SBATCH --partition=q40,q36,q28,q24
#SBATCH --mem=310G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=55:00:00


echo "========= Job started  at `date` =========="

K = 10
name = printf("K_%d", K)
mkdir /home/normann/programs/sparse_OU/data/$name


# Go to the directory where this job was submitted
cd $SLURM_SUBMIT_DIR

# copy inputdata and the executable to the scratch-directory
cp strength_OU.py /scratch/$SLURM_JOB_ID
cp -r ../sparse_OU_library/ /scratch/$SLURM_JOB_ID
cp ../ehrenfest_chain.py /scratch/$SLURM_JOB_ID


# change directory to the local scratch-directory:
cd /scratch/$SLURM_JOB_ID
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
python3 strength_OU.py $K 

# copy home the outputdata:
cp *.pkl /home/normann/programs/sparse_OU/data/$name

echo "========= Job finished at `date` =========="
#
