#!/bin/bash

#SBATCH --job-name=alignparams
#SBATCH --time=7-00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --output=sm.slurm.out
#SBATCH --error=sm.slurm.err
#SBATCH --partition=io

# Dependencies
module load apptainer/1.1.5-suid

# clobber broken apptainer config from module
export APPTAINER_BINDPATH=""

# Application specific commands:
export TMPDIR="${JOBDIR}"

printf "JOBDIR: %s\n" "${JOBDIR}"
printf "LOCALDIR: %s\n" "${LOCALDIR}"
printf "MEMDIR: %s\n" "${MEMDIR}"
printf "TMPDIR: %s\n" "${TMPDIR}"

snakemake \
    --profile profiles/petrichor_memdir \
    --retries 0 \
    --keep-going \
    --cores 64 \
    --local-cores 2 \
    output/010_captus-align/NUC.NT.mafft_auto.gap0.0.cov0.0
