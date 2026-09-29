#!/bin/bash
#SBATCH --job-name=AppDp16_hyp_fastp_array
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jela6408@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --array=0-23
#SBATCH --output=/scratch/alpine/jela6408/APPDp16_Hyp/fastp/%x_%A_%a.out
#SBATCH --error=/scratch/alpine/jela6408/APPDp16_Hyp/fastp/%x_%A_%a.err

echo "Job started at $(date)"
echo "This is array task ID: $SLURM_ARRAY_TASK_ID"

# load Modules
module load gcc

export PATH=$PATH:/projects/gonzalpk/
module --ignore-cache load intel/16.0.3 openmpi/1.10.2 boost/1.61.0 python/2.7.11
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/curc/sw/intel/17.0.0/tbb/lib/intel64/gcc4.7

# setup Paths
raw_dir="/pl/archive/Link_Lab_archive/JLawrence/APPDp16_Hyp/"
out_dir="/scratch/alpine/jela6408/APPDp16_Hyp/fastp"

# gather all Forward read (R1) files matching the Hoeffer prefix into a bash array
r1_files=($(ls ${raw_dir}/*_R1_001.fastq.gz))

# Select the specific R1 file for this job core
r1_file=${r1_files[$SLURM_ARRAY_TASK_ID]}

# automatically deduce the corresponding Reverse read (R2) file name
r2_file="${r1_file/_R1_001.fastq.gz/_R2_001.fastq.gz}"

# extract just the clean prefix (e.g., Hoeffer_31_S19_L005) for naming outputs
r1_basename=$(basename "$r1_file")
sample_name="${r1_basename/_R1_001.fastq.gz/}"

echo "Processing Sample Name: $sample_name"
echo "Forward Read: $r1_file"
echo "Reverse Read: $r2_file"

# RUN FASTP
fastp \
    -i "$r1_file" \
    -I "$r2_file" \
    -o "${out_dir}/${sample_name}_R1_trimmed.fastq.gz" \
    -O "${out_dir}/${sample_name}_R2_trimmed.fastq.gz" \
--detect_adapter_for_pe \
    --trim_poly_g \
    --cut_mean_quality 20 \
    --cut_front \
    --cut_tail \
    --length_required 50 \
    --thread 4 \
    --html "${out_dir}/${sample_name}.fastp.html"

echo "Job finished at $(date)"