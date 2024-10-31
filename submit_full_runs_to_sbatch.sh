#!/bin/bash

#SBATCH --job-name=tryparams
#SBATCH --time=7-00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g
#SBATCH --output=full_run.slurm.out
#SBATCH --error=full_run.slurm.err
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
    --retries 1 \
    --keep-going \
    --cores 64 \
    --local-cores 2 \
    -s rerun_full_dataset.smk \
    output/090_rerun-full-dataset/025_trimal-processed/NUC.NT.mafft_auto.gap0.9.cov0.4/kept.tar
