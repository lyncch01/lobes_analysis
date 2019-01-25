#!/bin/bash -l
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --partition=workq
#SBATCH --account=mwaeor
#SBATCH -J procImg_093
#SBATCH --out=prImg.o%A
#SBATCH --error=prImg.e%A

module load aegean
module load MWA_Tools/mwa-sci


python process_images.py
