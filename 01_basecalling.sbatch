#!/bin/bash
#SBATCH --job-name=rumen_18s_basecalling
#SBATCH --partition=aoraki_gpu
#SBATCH --gpus-per-node=1                         # Request 1 GPU
#SBATCH --cpus-per-task=2                         # At least 2 CPUs per GPU
#SBATCH --mem=64GB                                # Adjust based on memory requirements
#SBATCH --time=128:00:00 
#SBATCH --output=output/logs/rumen_18s_basecalling.out
#SBATCH --error=output/logs/rumen_18s_basecalling.log

# Load the required module for GPU support (if needed)
module load cuda

# Change to the scratch directory for work
cd /weka/health_sciences/bms/biochemistry/dearden_lab/sarahinwood

mkdir -p rumen_18s_basecalling

singularity exec --nv containers/dorado_v0.9.5.sif dorado \
    basecaller \
    sup \
    --recursive /projects/health_sciences/bms/biochemistry/dearden_lab/sarahinwood/emma-rumen-18s/data/emma_18s_barcoding \
    --kit-name SQK-NBD114-24 \
    > rumen_18s_basecalling/all.bam

mkdir -p /projects/health_sciences/bms/biochemistry/dearden_lab/sarahinwood/emma-rumen-18s/output/01_basecalled

# Once the job is done, move results to the final destination
mv rumen_18s_basecalling/all.bam /projects/health_sciences/bms/biochemistry/dearden_lab/sarahinwood/emma-rumen-18s/output/01_basecalled/