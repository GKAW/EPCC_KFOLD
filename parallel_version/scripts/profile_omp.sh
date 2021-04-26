#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=OMP_Profile
#SBATCH --time=0:20:0
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1

# Replace [budget code] below with your budget code (e.g. t01)
#SBATCH --account=tc027
# We use the "standard" partition as we are running on CPU nodes
#SBATCH --partition=standard
# We use the "standard" QoS as our runtime is less than 4 days
#SBATCH --qos=standard

# Load any required modules
module load mpt
module load nvidia/cuda-11.2

# Change to the submission directory
cd $SLURM_SUBMIT_DIR

# Set the number of threads to the CPUs per task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Launch the parallel job
#   Using 36 threads per node
#   srun picks up the distribution from the sbatch options
srun --cpu-bind=cores nsys profile --stats=true --force-overwrite=true -o ../../profiles/omp-report ./../src_omp/Kfold.x -i ../../examples/rbs/rbs_0000.start -o ../../outputs/rbs_0000.traj -l ../../outputs/rbs_0000.log -s 1234 -n 1 -t 1000
