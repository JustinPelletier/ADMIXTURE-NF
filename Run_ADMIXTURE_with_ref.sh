#!/bin/bash

#SBATCH --job-name=IBD_segment
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=24:00:00
#SBATCH --output=Run_ADMIXTURE_with_ref.out

set -euo pipefail

module load nextflow

export NXF_OFFLINE=true


nextflow run ADMIXTURE_with_ref.nf \
    -c ADMIXTURE_with_ref.config \
    -work-dir $(pwd)/log #\
    #-resume
