#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm-nf-meta-%j.out
set -euo pipefail

# Load the Nextflow Conda environment
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/nextflow

# Point to the main workflow file/dir
WORKFLOW=main.nf    # Eventually, this should point to a GitHub repo

# Report
date
echo

# Run the workflow
nextflow run $WORKFLOW -ansi-log false -resume "$@"

# Report
echo
date
