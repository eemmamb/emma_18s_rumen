#!/bin/bash
#SBATCH --job-name=nanoplot_postchopper
#SBATCH --output=output/logs/nanoplot_postchopper.log
#SBATCH --error=output/logs/nanoplot_postchopper.log
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

# Activate virtual environment
source venv/bin/activate

# Input files
INPUT1="output/04_chopper/barcode01.fastq"
INPUT2="output/04_chopper/barcode02.fastq"
INPUT3="output/04_chopper/barcode03.fastq"
INPUT4="output/04_chopper/barcode04.fastq"

# Output base
OUTBASE="output/04_chopper/nanoplot_post_chopper"
mkdir -p "$OUTBASE"

# Run NanoPlot for each sample with a different line/point color
echo "Running NanoPlot on barcode01"
NanoPlot --fastq "$INPUT1" --outdir "$OUTBASE/sample1" --threads 2 \
 --color "cornflowerblue" --plots dot kde --N50 --title "Sample 1"

echo "Running NanoPlot on barcode02"
NanoPlot --fastq "$INPUT2" --outdir "$OUTBASE/sample2" --threads 2 \
  --color "lightgreen" --plots dot kde --N50 --title "Sample 2"

echo "Running NanoPlot on barcode03"
NanoPlot --fastq "$INPUT3" --outdir "$OUTBASE/sample3" --threads 2 \
  --color "orange" --plots dot kde --N50 --title "Sample 3"

echo "Running NanoPlot on barcode04"
NanoPlot --fastq "$INPUT4" --outdir "$OUTBASE/sample4" --threads 2 \
  --color "hotpink" --plots dot kde --N50 --title "Sample 4"

# ## all samples together
 echo "Running NanoPlot on all samples together"
 cat "$INPUT1" "$INPUT2" "$INPUT3" "$INPUT4" > "$OUTBASE/all_samples.fastq"
NanoPlot --fastq "$OUTBASE/all_samples.fastq" --outdir "$OUTBASE/all_samples" --threads 2 \
  --color "violet" --plots dot kde --N50 --title "All Samples" 

echo "NanoPlot runs completed."
