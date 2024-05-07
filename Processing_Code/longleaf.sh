#!/bin/bash
#SBATCH -t 2-00:00:00
#SBATCH -n 16
#SBATCH --mem=300g 
#SBATCH --mail-user=mhg@ad.unc.edu
#SBATCH --mail-type=end,fail,begin
#SBATCH --job-name=D1
#SBATCH --export=ALL

cd /proj/hyejunglab/singlecell/Sangmi/Processing_Code

module add anaconda
conda activate R_envs

# Rscript B1.Annot_SC001.R
# Rscript B2.Annot_SC002.R
# Rscript C1.Integrating_SC001.R
# Rscript C2.Integrating_SC002.R
# Rscript C3.Integrating_All.R
Rscript D1.celltype_analysis.R
