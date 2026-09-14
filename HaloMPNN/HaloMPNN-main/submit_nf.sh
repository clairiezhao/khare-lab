#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --requeue
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --output=out/%x.slurm_out.%j-%2t.%a
#SBATCH --mem=11G
#SBATCH --time=1:00:00
#SBATCH --job-name=mpnnretrain

mkdir -p out
export NXF_LOG_FILE=out/.nextflow.log

# Initialize shared conda environment
. /projects/community/anaconda/2022.10/bd387/etc/profile.d/conda.sh
conda activate /projects/f_sdk94_1/conda/envs/nextflow

nextflow run main.nf -profile slurm --input_fasta test.fasta