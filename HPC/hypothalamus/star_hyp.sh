#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jela6408@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=01:00:00
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --output=/scratch/alpine/jela6408/APPDp16_Hyp/star_%A_%a.out
#SBATCH --error=/scratch/alpine/jela6408/APPDp16_Hyp/star_%A_%a.err
#SBATCH --array=0-23

echo Job started at `date +"%T %a %d %b %Y"`
module load star/2.7.10b

# define directories based on your current path
FASTQ_DIR="/scratch/alpine/jela6408/APPDp16_Hyp/fastp"
OUT_DIR="/scratch/alpine/jela6408/APPDp16_Hyp/star"

# move to the fastq directory to grab the filenames
cd ${FASTQ_DIR}

# create an array of all unique sample prefixes
SAMPLES=($(ls *_R1_trimmed.fastq.gz | sed 's/_R1_trimmed.fastq.gz//'))

# select the specific sample assigned to this specific array task
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "Currently processing sample: ${SAMPLE}"

# run STAR alignment

STAR --genomeDir /pl/archive/Link_Lab_archive/griffioen/grc39/ \
--readFilesIn ${FASTQ_DIR}/${SAMPLE}_R1_trimmed.fastq.gz ${FASTQ_DIR}/${SAMPLE}_R2_trimmed.fastq.gz \
--readFilesCommand zcat \
--runThreadN 12 ${SLURM_CPUS_PER_TASK} \
--runMode alignReads \
--outFileNamePrefix ${OUT_DIR}/${SAMPLE}_ \
--outSAMtype BAM Unsorted

echo "Job finished at $(date +"%T %a %d %b %Y")"