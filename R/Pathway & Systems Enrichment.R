setwd("~/R prac/GWAS_project/Scripts")

# Import the necessary libraries
library(readxl)
library(data.table)
library(tidyr)
library(dplyr)
library(stringr)

# ----- PATHWAYS AND SYSTEMS ENRICHMMENT -----

gs <- read_excel("GS.xlsx")

# Make the normalized description column
desc_norm <- tolower(gsub("_", " ", paste(gs$Category, gs$GeneSet, sep = ": ")))

# 3. Terms to search for (asthma/airway-related)
priority_terms <- c("asthma", "airway", "lung", "eosinophil", "th2", "allergy", "il13", "il4", "il5")
pattern <- paste(priority_terms, collapse = "|")

# 4. Find matches
matches <- grepl(pattern, desc_norm, ignore.case = TRUE) & gs$adjP < 0.05

# 5. Subset your table
asthma_sets <- gs[matches, ]

# 6. Print and export
print(asthma_sets)
fwrite(asthma_sets, "GWAS_Enriched_Pathways.csv")