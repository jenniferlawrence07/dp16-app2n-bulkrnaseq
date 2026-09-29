#!/bin/bash
#SBATCH --job-name=FeatureCountsAPPHyp
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jela6408@colorado.edu                        # Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=8                                      # Number of CPU (tasks)
#SBATCH --time=02:00:00                                 # Time limit hrs:min:sec
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --output=/scratch/alpine/jela6408/APPDp16_Hyp/counts/%x_%j.out
#SBATCH --error=/scratch/alpine/jela6408/APPDp16_Hyp/counts/%x_%j.err

echo Job started at `date +"%T %a %d %b %Y"`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

module load R

Rscript /scratch/alpine/jela6408/APPDp16_Hyp/scripts/featureCounts_Hyp.R
