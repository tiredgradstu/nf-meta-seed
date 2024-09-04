#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm-nf-meta-%j.out
set -euo pipefail

# Load the Nextflow Conda environment
module load miniconda3/24.1.2-py310
conda activate /users/PAS2693/zenmckenzie14/miniconda3/envs/nextflow

# Point to the main workflow file/dir
WORKFLOW=git@github.com:zenmck/nf-meta-seed.git   
# Constants
#[[ -n $SLURM_JOB_ACCOUNT ]] && osc_account=$(echo "$SLURM_JOB_ACCOUNT" | tr "[:lower:]" "[:upper:]") //ask about this
WORKFLOW=/fs/ess/PAS2693/jelmer/meta_pipeline
WORKDIR=/fs/scratch/PAS2693/nf-meta
OUTDIR=results/nf-meta

# Report
echo
date
echo -e "\n# Starting Nextflow run with Nextflow base call:"
echo "nextflow run $WORKFLOW -ansi-log false -resume -work-dir $WORKDIR" 
echo -e "\n# ... and with pipeline parameters:"
echo "--outdir $OUTDIR $*"
echo -e "\n==========================================\n"

# Run the workflow
nextflow run $WORKFLOW \
    -ansi-log false \
    -resume \
    --outdir "$OUTDIR" \
    -work-dir "$WORKDIR" \
    "$@"

# Report
echo
date

# Run the workflow
nextflow run $WORKFLOW -ansi-log false -resume "$@"
