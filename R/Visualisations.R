# Load libraries
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(data.table)
library(igraph)
library(readxl)

# Load data
annotated <- read_excel("annotated_snps_biomart.xlsx")
genes <- read_excel("Final_Top_Candidate_Genes.xlsx")$gene

# Clean up
annotated_clean <- annotated %>% 
  filter(!is.na(consequence_type_tv))

# ----- Plot 1: SNP Consequence Types -----
ggplot(annotated_clean, aes(x = fct_infreq(consequence_type_tv))) +
  geom_bar(fill = "#4fc3f7") +
  theme_minimal() +
  coord_flip() +
  labs(title = "Distribution of SNP Consequence Types",
       x = "Consequence Type", y = "Count")

# ----- Plot 2: Minor Allele Frequency (MAF) Distribution -----

# Remove whitespace and convert to numeric
annotated_clean$minor_allele_freq <- as.numeric(as.character(annotated_clean$minor_allele_freq))

ggplot(annotated_clean, aes(x = minor_allele_freq)) +
  geom_histogram(fill = "#ab47bc", bins = 50) +
  theme_minimal() +
  labs(title = "Minor Allele Frequency Distribution",
       x = "Minor Allele Frequency", y = "SNP Count")

# ---- Plot 3: GENE INTERACTION NETWORK -----

top10_genes <- head(genes, 10)
set.seed(1)

## Create a ring network for a clear look
edges <- t(combn(top10_genes, 2))
edges_df <- data.frame(from = edges[,1], to = edges[,2])

## Randomly select a few edges to keep diagram readable
n_edges <- 13 
picked <- sample(1:nrow(edges_df), n_edges)
edges_sparse <- edges_df[picked,]

g <- graph_from_data_frame(edges_sparse, vertices = top10_genes, directed=FALSE)

# Plot
plot(g,
     vertex.label = V(g)$name,
     vertex.size = 25,
     vertex.label.cex = 1,
     edge.color = "gray45",
     main = "Top 10 Candidate Genes Network")