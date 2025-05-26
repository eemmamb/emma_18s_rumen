#!/bin/bash
#SBATCH --job-name=rumen_18s_demux
#SBATCH --cpus-per-task=50
#SBATCH --mem=64GB
#SBATCH --time=128:00:00 
#SBATCH --output=output/logs/rumen_18s_demux.out
#SBATCH --error=output/logs/rumen_18s_demux.log

singularity exec containers/dorado_v0.9.5.sif dorado \
    demux \
    --output-dir output/02_demuxed \
    --no-classify \
    --threads 50 \
    --emit-fastq \
    output/01_basecalled/all.bam

for file in output/02_demuxed/550f7d57-7e3e-4ebe-a1e4-37d432de2edf_SQK-NBD114-24_barcode*.fastq; do
    filename=$(basename "$file")  # get just the filename
    newname="${filename#550f7d57-7e3e-4ebe-a1e4-37d432de2edf_SQK-NBD114-24_}"  # strip the prefix
    mv "$file" "output/02_demuxed/$newname"  # move to new name in same folder
done