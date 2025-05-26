#!/bin/bash
#SBATCH --job-name=vsearch_clustering
#SBATCH --cpus-per-task=50
#SBATCH --mem=64GB
#SBATCH --time=128:00:00 
#SBATCH --output=output/logs/vsearch_clustering.log
#SBATCH --error=output/logs/vsearch_clustering.log

# Define input/output
INPUT_DIR="output/04_chopper"
OUTPUT_DIR="output/05_vsearch_clustering"
mkdir -p "$OUTPUT_DIR"

# List of files (adjust filenames as needed)
FILES=("barcode01.fastq" "barcode02.fastq" "barcode03.fastq" "barcode04.fastq")

# an illumina workflow would typically do a dereplication step - collapsing identical reads before any clustering is done but with higher error rates
# with nanopore error rates and longer reads (i.e. noisier reads) we don't know that what should be identical sequence will actually match 100%

## here we cluster noisy reads to make rough/draft OTUs that we can then polish to reduce sequencing errors

for FILE in "${FILES[@]}"; do
    BASENAME=$(basename "$FILE" .fastq)
    echo "VSEARCH clustering $FILE"
    singularity exec containers/vsearch_2.30.0--hd6d6fdc_0.sif \
        vsearch \
        --cluster_unoise "$INPUT_DIR/$FILE" \
        --minsize 1 \
        --threads 50 \
        --relabel "${BASENAME}_" \
        --centroids "$OUTPUT_DIR/${BASENAME}_otus.fa"
done
# have to rename output with relabel so that draft OTU names aren't the same as read names for running racon

echo "Concatenating clustered OTU's" 

cat output/05_vsearch_clustering/barcode*_otus.fa > output/05_vsearch_clustering/all_otus.fa
