# Import libraries
library(data.table)
library(dplyr)
library(readxl)

# 1. Read SNPs 
snps <- read_excel("SNPs.xlsx")

# 2. Read Pharos + OMIM filtered gene table 
pharos_omim_genes <- read_excel("pharos.xlsx")$Symbol

# 3. Read the pathway enrichment file and extract key pathway gene set members
asthma_sets <- read_excel("GWAS_Enriched_Pathways.xlsx")
genes_in_pathways <- unique(unlist(strsplit(paste(asthma_sets$genes, collapse = ":"), split = "[:,\\s]+")))

# 4. Read eQTL prioritization results 
top_eqtl <- read_excel("Top_Prioritized_eQTLs.xlsx")$gene 

# 5. SNP nearest gene mapping 
multi_map_genes <- snps$nearestGene %>% unique()

# ---- Combine & Prioritize ----

# Make a list of all sources
gene_sources <- list(
  SNPs = multi_map_genes,
  Pharos_OMIM = pharos_omim_genes,
  Pathways = genes_in_pathways,
  eQTL = top_eqtl
)

# Convert to long format (gene, source)
all_genes <- stack(gene_sources) %>%
  rename(gene = values, source = ind) %>%
  distinct()

# Count how many sources support each gene
gene_counts <- all_genes %>%
  group_by(gene) %>%
  summarise(supporting_sources = n(), sources = paste(source, collapse = ";"), .groups = "drop")

# Choose final candidates = genes supported by ≥2 sources 
final_candidates <- gene_counts %>% filter(supporting_sources >= 2)

# Save results
fwrite(final_candidates, "Top_Candidate_Genes.xlsx")

