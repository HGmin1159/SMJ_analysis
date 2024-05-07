#!/bin/bash
#SBATCH -t 2-00:00:00
#SBATCH --mem=16g 
#SBATCH --mail-user=mhg@ad.unc.edu
#SBATCH --mail-type=end,fail,begin
#SBATCH --job-name=Setup_Environment_Sangmi
#SBATCH --export=ALL

module add anaconda
# conda create -n R_envs r-essentials r-base
conda activate R_envs
conda install -c conda-forge r-seurat

R -e "if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')"
R -e "BiocManager::install(version = '3.16')"
R -e "BiocManager::install('SingleR')"
R -e "BiocManager::install('scater')"
