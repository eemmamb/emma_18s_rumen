#!/bin/bash
#SBATCH --job-name=chopper
#SBATCH --cpus-per-task=50
#SBATCH --mem=64GB
#SBATCH --time=128:00:00 
#SBATCH --output=output/logs/chopper.log
#SBATCH --error=output/logs/chopper.log

# Define input/output
INPUT_DIR="output/02_demuxed"
OUTPUT_DIR="output/04_chopper"
mkdir -p "$OUTPUT_DIR"

# List of files (adjust filenames as needed)
DEMUXED_FILES=("barcode01.fastq" "barcode02.fastq" "barcode03.fastq" "barcode04.fastq")

## chopper - filter reads based on length and quality
for FILE in "${DEMUXED_FILES[@]}"; do
    BASENAME=$(basename "$FILE" .fastq)
    echo "Chopper filtering $FILE"
    singularity exec containers/chopper_0.1.0.sif \
        /usr/local/bin/chopper \
        --input "$INPUT_DIR/$FILE" \
        --threads 50 \
        --quality 15 \
        --minlength 430 \
        --maxlength 530 \
        > "$OUTPUT_DIR/${BASENAME}_filtered.fastq"
done

FILTERED_FILES=("barcode01_filtered.fastq" "barcode02_filtered.fastq" "barcode03_filtered.fastq" "barcode04_filtered.fastq")

## now trim primer sequence off reads - couldn't tell whether this happened before or after above filtering
for FILE in "${FILTERED_FILES[@]}"; do
    BASENAME=$(basename "$FILE" .fastq)
    echo "Chopper cropping primer seq $FILE"
    singularity exec containers/chopper_0.1.0.sif \
        /usr/local/bin/chopper \
        --input "$OUTPUT_DIR/$FILE" \
        --threads 50 \
        --headcrop 20 \
        --tailcrop 19 \
        > "$OUTPUT_DIR/${BASENAME}_cropped.fastq"
done

## min and max length in this script reflect expected amplicon size +/- 50 bp BEFORE head and tail crop to remove primers
## porechop should filter reads based on length/qual and then crop the reads that passed filters

## fwd primer - P-SSU-316F: GCTTTCGWTGGTAGTGTATT
## rev primer - GIC758R:    CAACTGTCTCTATKAAYCG