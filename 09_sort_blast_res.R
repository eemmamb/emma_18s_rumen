library(data.table)
library(tidyr)
library(dplyr)

blast_res <- fread("output/08_OTU_taxonomic_annotation/blastn/blastn_cluster_otus_ID_97.outfmt6")

setnames(blast_res, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("query", "subject", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))

blast_res$otu_no <- tstrsplit(blast_res$query, ";", keep=1)
blast_res$otu_no <- tstrsplit(blast_res$otu_no, "_", keep=2)
blast_res$otu_no <- as.numeric(blast_res$otu_no)
blast_res$otu_size <- tstrsplit(blast_res$query, ";", keep=2)
blast_res$otu_size <- tstrsplit(blast_res$otu_size, "=", keep=2)
blast_res$otu_size <- as.numeric(blast_res$otu_size)

##################################
## NOT filtering out uncultured ##
##################################

##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(blast_res, otu_no, evalue, -bit_score, -`%_identical_matches`)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues_no_filtering <- blast_res[,.SD[which.min(evalue)], by=otu_no]
fwrite(min_evalues_no_filtering, "output/08_OTU_taxonomic_annotation/blastn/best_res_per_otu_no_filtering.csv")

##############################
## filtering out uncultured ##
##############################

blast_res_filtered <- subset(blast_res, !grepl("Uncultured", blast_res$annotation))

##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(blast_res_filtered, otu_no, evalue, -bit_score, -`%_identical_matches`)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues_filtered <- blast_res_filtered[,.SD[which.min(evalue)], by=otu_no]
fwrite(min_evalues_filtered, "output/08_OTU_taxonomic_annotation/blastn/best_res_per_otu_uncultured_removed.csv")

fwrite(list(min_evalues_filtered$taxid), "output/08_OTU_taxonomic_annotation/blastn/filtered_hits_taxids.txt")

############# pause here and run scripts/09_taxonomy_for_blast_res.sbatch #############

taxonomy <- fread("output/08_OTU_taxonomic_annotation/blastn/formatted_taxonomy.tsv")
taxonomy <- taxonomy[,c(1,3)]
unique_taxonomy <- unique(taxonomy)

min_evalues_filtered_w_taxonomy <- merge(min_evalues_filtered, unique_taxonomy,
                                         by.x="taxid", by.y="V1")

phyloseq_taxonomy_table <- min_evalues_filtered_w_taxonomy[,c(2,17)]
phyloseq_taxonomy_table$otu <- paste("OTU_", phyloseq_taxonomy_table$otu_no, sep="")
phyloseq_taxonomy_table <- phyloseq_taxonomy_table[,c(3,2)]

separated_taxonomy_phyloseq <- phyloseq_taxonomy_table %>% 
  separate(V3, into = paste0("Rank", 1:13),
            sep = ";",
            fill = "right") %>%
  mutate(across(everything(), ~ trimws(.)))

phylo_ranks <- separated_taxonomy_phyloseq %>%
  transmute(
    OTU = otu,
    Domain = Rank2,
    Supergroup = Rank3,
    Superphylum = Rank4,
    Phylum = Rank5,
    Subphylum = Rank6,
    Class = Rank7,
    Subclass = Rank8,
    Order = Rank9,
    Family = Rank10,
    Genus = Rank11,
    Species = Rank12,
    Subspecies = Rank13)

phylo_ranks[is.na(phylo_ranks)] <- "unclassified"

fwrite(phylo_ranks, "output/09_for_phyloseq/tax_table_for_phyloseq.csv")
