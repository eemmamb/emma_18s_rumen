library(phyloseq)
library(data.table)
library(microbiome)
library(viridis)

# Step 1: Convert count table
phyloseq_count_table <- fread("output/07_cluster_polished_OTUs/cluster_otus_ID_97.tsv")
otu_counts <- as.matrix(phyloseq_count_table[, -1, with = FALSE])  # exclude first column (OTU IDs)
rownames(otu_counts) <- phyloseq_count_table[[1]]  # set rownames to OTU IDs

otu_table_ps <- otu_table(otu_counts, taxa_are_rows = TRUE)

# Step 2: Convert taxonomy table
phyloseq_tax_table <- fread("output/09_for_phyloseq/tax_table_for_phyloseq.csv")
tax_mat <- as.matrix(phyloseq_tax_table[, -1, with = FALSE])  # exclude OTU column
rownames(tax_mat) <- phyloseq_tax_table[[1]]  # set rownames to OTU IDs

tax_table_ps <- tax_table(tax_mat)

# Step 3: sample data
# Extract sample names from your count table
sample_names <- colnames(phyloseq_count_table)[-1]  # skip first column with OTU IDs
# Create a simple metadata data.frame
sample_metadata <- data.frame(SampleID = sample_names, row.names = sample_names)
sample_metadata$sample_name <- paste0("Sample ", seq_len(nrow(sample_metadata)))
# Convert to phyloseq sample_data object
sample_data_ps <- sample_data(sample_metadata)


# Step 3: Combine into phyloseq object
phyloseq_object <- phyloseq(otu_table_ps, tax_table_ps, sample_data_ps)

###############
# composition # https://microbiome.github.io/tutorials/Composition.html
###############

# transform to compositional
ps_comp <- microbiome::transform(phyloseq_object, transform="compositional")
## aggregate to class (87.5% ASVs annotated) - below this only 66% annotated
## here you can pick different levels of taxonomy
ps_plot <- ps_comp %>%
  microbiome::aggregate_rare(level="Genus", detection=0.01, prevalence=0.1) %>%
  microbiome::aggregate_taxa(level="Genus")

microbiome::plot_composition(ps_plot,
                             otu.sort = "abundance",
                             x.label="sample_name") +
  geom_col(width=1)+ # removes white gap between bars
  scale_fill_viridis(discrete=T) +
  scale_y_continuous(label = scales::percent)+
  labs(fill = "Genus")+
  theme_bw()
