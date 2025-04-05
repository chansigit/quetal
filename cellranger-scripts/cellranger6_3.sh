#!/bin/bash
#SBATCH --job-name=CR-pool6_3
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=7:40:00
#SBATCH -p bigmem

/home/users/chensj16/oak/software/cellranger-8.0.1/cellranger multi \
  --id=pool6_3 --csv=/home/users/chensj16/s/tmp/yimiao/cfg-pool6_3.csv \
  --localcores=$SLURM_CPUS_PER_TASK \
  --localmem=250