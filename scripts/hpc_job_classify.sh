#!/bin/bash -l

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=60G
#SBATCH --job-name="de_novo_gtdb"
#SBATCH --mail-user=s.paragkamian@hcmr.gr
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --requeue

start=`date +%s`

# set the environment
module purge # unloads all previous loads

module load anaconda3/default
conda activate /home1/s.paragkamian/software/conda_envs/gtdb-2-3-2

# run the script
/home1/s.paragkamian/bacillus/scripts/de_novo_gtdb.sh

module purge

# summary of job
end=`date +%s`
runtime=$((end-start))
echo "Job ID: " $SLURM_JOB_ID
echo "Job name: " $SLURM_JOB_NAME
echo $runtime "in seconds" 
echo $((runtime/60)) "in minutes" 
