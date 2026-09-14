#!/bin/bash
#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=out/%x.slurm_out.%j-%2t.%a
#SBATCH --mem=11G
#SBATCH --time=1:00:00
#SBATCH --job-name=mpnnretrain

mkdir -p out
export NXF_LOG_FILE=out/.nextflow.log

# Initialize shared conda environment
. /projects/community/anaconda/2022.10/bd387/etc/profile.d/conda.sh
conda activate /projects/f_sdk94_1/conda/envs/nextflow

nextflow run main.nf -profile local --input_fasta test.fasta -resume --num_epochs 5