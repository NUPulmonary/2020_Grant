#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --mem=48G
#SBATCH --ntasks-per-node=12
#SBATCH --mail-user=rogangrant2022@u.northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name="200604_SCRIPT_batch9_bcl2fastq"

module purge
cd /projects/b1038/Pulmonary/rgrant/script_bulk/200604_NB501488_0362_AHJ3T3BGXF
module load bcl2fastq/2.19.1 
bcl2fastq --loading-threads 4 --writing-threads 4 --processing-threads 4
