#!/bin/bash
#SBATCH --job-name=lasso
#SBATCH --output=logs/log.out
#SBATCH --error=logs/log.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

conda init --all
conda activate /cluster/home/julian/.conda/envs/pmlb

INPUT_NOTEBOOK="notebooks/6-common_genes_models.ipynb"
OUTPUT_NOTEBOOK="notebooks/6-common_genes_models_exc.ipynb"

jupyter nbconvert --to notebook \
                  --execute "$INPUT_NOTEBOOK" \
                  --output "$OUTPUT_NOTEBOOK" \
                  --ExecutePreprocessor.timeout=-1

echo "Done, saved to $OUTPUT_NOTEBOOK"
