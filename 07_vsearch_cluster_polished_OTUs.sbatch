#!/bin/bash
#SBATCH --job-name=07_cluster_polished_OTUs
#SBATCH --cpus-per-task=80
#SBATCH --mem=64GB
#SBATCH --time=128:00:00 
#SBATCH --output=output/logs/07_cluster_polished_OTUs.log
#SBATCH --error=output/logs/07_cluster_polished_OTUs.log

## here we perform the re-cluster, this time using polished draft OTUs

mkdir -p output/07_cluster_polished_OTUs/cluster_size

#### cluster polished OTUS from all samples that have been relabelled and combined into a single file - cluster based on seq-id ####
## doing this using centroid-based clustering ##


## cluster with percent ID = 95%
echo "Running VSEARCH cluster_size with ID = 95%"
    singularity exec containers/vsearch_2.30.0--hd6d6fdc_0.sif \
    vsearch \
    --cluster_size output/06_racon_polishing/polished_relabelled_combined_OTUs.fa \
    --id 0.95 \
    --threads 80 \
    --relabel OTU_ \
    --sizeout \
    --otutabout output/07_cluster_polished_OTUs/cluster_otus_ID_95.tsv \
    --centroids output/07_cluster_polished_OTUs/cluster_otus_ID_95.fa
echo "Finished VSEARCH cluster_size with %ID = 95"


## cluster with percent ID = 97%
echo "Running VSEARCH cluster_size with ID = 97%"
    singularity exec containers/vsearch_2.30.0--hd6d6fdc_0.sif \
    vsearch \
    --cluster_size output/06_racon_polishing/polished_relabelled_combined_OTUs.fa \
    --id 0.97 \
    --threads 80 \
    --relabel OTU_ \
    --sizeout \
    --otutabout output/07_cluster_polished_OTUs/cluster_otus_ID_97.tsv \
    --centroids output/07_cluster_polished_OTUs/cluster_otus_ID_97.fa
echo "Finished VSEARCH cluster_size with %ID = 97"

## cluster with percent ID = 99%
echo "Running VSEARCH cluster_size with ID = 99%"
    singularity exec containers/vsearch_2.30.0--hd6d6fdc_0.sif \
    vsearch \
    --cluster_size output/06_racon_polishing/polished_relabelled_combined_OTUs.fa \
    --id 0.99 \
    --threads 80 \
    --relabel OTU_ \
    --sizeout \
    --otutabout output/07_cluster_polished_OTUs/cluster_otus_ID_99.tsv \
    --centroids output/07_cluster_polished_OTUs/cluster_otus_ID_99.fa
echo "Finished VSEARCH cluster_size with %ID = 99"

## could then do chimera detection and then re-cluster if we are worried about chimeras