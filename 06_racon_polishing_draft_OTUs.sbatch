#!/bin/bash
#SBATCH --job-name=06_racon_polishing
#SBATCH --cpus-per-task=80
#SBATCH --mem=64GB
#SBATCH --time=128:00:00 
#SBATCH --output=output/logs/06_racon_polishing.log
#SBATCH --error=output/logs/06_racon_polishing.log

## here we polish the noisy draft OTUs for each sample separately to reduce error rates

SAMPLES=("barcode01" "barcode02" "barcode03" "barcode04")

READ_INDIR="output/04_chopped"
INDIR="output/05_vsearch_clustering"
OUTDIR="output/06_racon_polishing"

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"



### Alignment of each samples filtered reads to the same sample's OTUs ####
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Running minimap on $SAMPLE"
    singularity exec containers/minimap2_2.29.sif \
        minimap2 \
        -ax map-ont \
        -t 80 \
        --secondary=no \
        "$INDIR/${SAMPLE}_otus.fa" "$READ_INDIR/${SAMPLE}.fastq" \
        > "$OUTDIR/${SAMPLE}.sam"
done
echo "mapping finished"



#### racon to use the mapped reads to polish OTU sequences (making them more accurate) ####
echo "Starting Racon polishing"
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Running racon on $SAMPLE"
    singularity exec containers/racon_1.5.0.sif \
        racon \
        "$READ_INDIR/${SAMPLE}.fastq" \
        "$OUTDIR/${SAMPLE}.sam" \
        "$INDIR/${SAMPLE}_otus.fa" \
        -t 80 \
        > "$OUTDIR/${SAMPLE}_polished_otus.fa"
done
echo "polishing finished"
# first input file contains sequences used for correction
# second input file contains overlaps between sequences and target sequences
# third input file contains sequences which will be corrected



#### vsearch - sorting sample OTUs by abundance(size), and renaming them to include the barcode0X in their name ####
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Running VSEARCH on $SAMPLE"
    singularity exec containers/vsearch_2.30.0--hd6d6fdc_0.sif \
        vsearch \
        --sortbysize "$OUTDIR/${SAMPLE}_polished_otus.fa" \
        --sample "$SAMPLE" \
        --output "$OUTDIR/${SAMPLE}_relabelled.fasta"
done
echo "relabelling finished"


#### cat polished and renamed OTUs into a single file ####
cat "$OUTDIR"/*_relabelled.fasta > "$OUTDIR/polished_relabelled_combined_OTUs.fa"
