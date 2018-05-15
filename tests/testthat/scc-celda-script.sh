#!/bin/bash -l

# Specify the version of R to be used $module avail
module load R/3.4.3
module load texlive
module load pandoc/2.2.1

#$ -l h_rt=24:00:00                   # Specify the hard time limit for the job
#$ -N celda_c_dev_profile             # Give job a name
#$ -j y                               # Merge the error and output streams into a single file

Rscript -e "options(keep.source=TRUE); 'run_profvis_DEV.R'"
