#!/bin/bash
#SBATCH --job-name=08_OTU_taxonomic_annotation_blastn
#SBATCH --account=inwsa89p
#SBATCH --cpus-per-task=80
#SBATCH --mem=64GB
#SBATCH --output=output/logs/08_OTU_taxonomic_annotation_blastn.log
#SBATCH --error=output/logs/08_OTU_taxonomic_annotation_blastn.log

mkdir -p output/08_OTU_taxonomic_annotation
mkdir -p output/08_OTU_taxonomic_annotation/blastn

echo "blastn search for final otus"
singularity exec containers/blast_2.16.0.sif \
    blastn \
    -outfmt "6 std staxids salltitles" \
    -db data/nt_2023_03_20/nt \
    -num_threads 80 \
    -query output/07_cluster_polished_OTUs/cluster_otus_ID_97.fa \
    -out output/08_OTU_taxonomic_annotation/blastn/blastn_cluster_otus_ID_97.outfmt6