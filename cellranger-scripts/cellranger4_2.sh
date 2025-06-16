#!/bin/bash
#SBATCH --job-name=CR-pool1_2
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=4:00:00
#SBATCH -p bigmem

/home/users/chensj16/oak/software/cellranger-8.0.1/cellranger multi \
  --id=pool4_2 --csv=/home/users/chensj16/s/tmp/yimiao/cfg-pool4_2.csv \
  --localcores=$SLURM_CPUS_PER_TASK \
  --localmem=250
