setwd("~/R prac/GWAS_project/Scripts")

# Import the necessary packages 
library(readxl)
library(tidyr)
library(data.table)
library(dplyr)

# Read the files
leadSNPs <- read_excel("leadSNPs.xlsx")
SNPs <- read_excel("SNPs.xlsx")
eqtl <- read_excel("eQTL.xlsx")


# ----- IDENTIFICATION AND BIOLOGICAL PRIORITISATION OF TOP SNPs -----

## (a) Identify significant Lead SNPs by p-value 
lead_top <- leadSNPs %>% filter(p < 5e-8)

## (b) Merge feature annotation from SNPs onto Lead SNPs
lead_exploded <- leadSNPs %>%
  separate_rows(IndSigSNPs, sep = ";") %>%  
  mutate(IndSigSNPs = trimws(IndSigSNPs)) %>%
  left_join(SNPs, by = c("IndSigSNPs" = "IndSigSNP")) 

## (c) Prioritize on biological relevance
top_snps_prioritized <- lead_exploded %>%
  filter(
    p < 5e-8,                # p-value
    MAF > 0.01,              # Common variants
    CADD >= 12.37,           # Deleteriousness (CADD Scores)
    RDB < 4 | RDB == 1       #                 (RegulomeDB Scores)
  )

## Save for downstream
fwrite(top_snps_prioritized, "Top_Prioritized_SNPs.csv")

# ----- TISSUE EXPRESSION AND PRIORITISATION -----

top_eqtl <- eqtl %>%
  filter(eqtlMapFilt == 1) %>%
  group_by(gene, tissue) %>%
  summarise(min_p = min(p)) %>%
  arrange(min_p)

## Save for downstream
fwrite(top_eqtl, "Top_Prioritized_eQTLs.csv")