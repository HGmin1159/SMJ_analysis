#!/bin/bash
#SBATCH -t 2-00:00:00
#SBATCH --mem=500g 
#SBATCH --mail-user=mhg@ad.unc.edu
#SBATCH --mail-type=end,fail,begin
#SBATCH --job-name=Download_file_Sangmi
#SBATCH --export=ALL

data_url=http://ngs.c2b2.columbia.edu/ngs/release/singleCell/240118_SANGMI_CAIJING_1_HUMAN_10X/
data_dir=/proj/hyejunglab/singlecell/Sangmi/Data1

mkdir -p "$data_dir"
cd "$data_dir"

wget -r -np -nH --cut-dirs=4 -R "index.html*" "$data_url"

tar -xvf /proj/hyejunglab/singlecell/Sangmi/Data1/SC001/analysis/240118_SANGMI_CAIJING_1_HUMAN_10X-SC001-cellranger-count-default/SC001_cellranger_count_outs.tar
tar -xzvf /proj/hyejunglab/singlecell/Sangmi/Data1/SC001/analysis/240118_SANGMI_CAIJING_1_HUMAN_10X-SC001-cellranger-count-default/SC001_cellranger_count_info.tar.gz
tar -xvf /proj/hyejunglab/singlecell/Sangmi/Data1/SC002/analysis/240118_SANGMI_CAIJING_1_HUMAN_10X-SC002-cellranger-count-default/SC002_cellranger_count_outs.tar
tar -xzvf /proj/hyejunglab/singlecell/Sangmi/Data1/SC002/analysis/240118_SANGMI_CAIJING_1_HUMAN_10X-SC002-cellranger-count-default/SC002_cellranger_count_info.tar.gz

gunzip /proj/hyejunglab/singlecell/Sangmi/Data/GSE208672_Seurat_allsamples.rds.gz