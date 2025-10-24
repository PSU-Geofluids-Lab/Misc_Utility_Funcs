#!/bin/bash
#SBATCH -p urseismo -A tolugboj_lab
#SBATCH -n 1
#SBATCH --mem=40G
#SBATCH -t 30:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=user@ur.rochester.edu
#SBATCH -o /dev/null

module unload python/2.7.12
source /scratch/tolugboj_lab/softwares/anaconda/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate instaseis

srun --output=job_status/job_${3}_${4}_%j.out python Single_Station_Data_Download_Workflow.py "$1" "$2" "$3" "$4"


