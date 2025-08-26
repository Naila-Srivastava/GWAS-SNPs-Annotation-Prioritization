# Load the required libraries
library(readxl)
library (biomaRt)

# Import and load the dataset
gwas_significant_snps <- read_excel("C:\\Users\\NAILA\\OneDrive\\Documents\\R prac\\GWAS_project\\gwas_significant_snps.xlsx")
snp_df <- gwas_significant_snps

head(snp_df)

# Filter for valid rsIDs only 
snp_df <- snp_df[grep("^rs", snp_df$SNP), ]
cat("Number of SNPs with valid rsIDs:", nrow(snp_df), "\n")

# Connect to Ensembl biomart
ensembl <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")

# Define attributes to retrieve
attributes <- c("refsnp_id", "chr_name", "chrom_start", "chrom_end",
                "consequence_type_tv", "allele", "minor_allele",
                "minor_allele_freq", "ensembl_gene_stable_id", "associated_gene")

# Query biomaRt in batches (avoid timeouts)
snp_list <- unique(snp_df$SNP)
batch_size <- 500
results_list <- list()

cat("Starting annotation...\n")

for (i in seq(1, length(snp_list), by = batch_size)) {
  batch <- snp_list[i:min(i + batch_size - 1, length(snp_list))]
  cat("Annotating SNPs", i, "to", min(i + batch_size - 1, length(snp_list)), "...\n")
  
  res <- tryCatch({
    getBM(attributes = attributes,
          filters = "snp_filter",
          values = batch,
          mart = ensembl)
  }, error = function(e) {
    warning("Batch failed:", conditionMessage(e))
    return(NULL)
  })
  
  if (!is.null(res) && nrow(res) > 0) {
    results_list[[length(results_list) + 1]] <- res
  }
}

# Combine all batches
annotated_snps <- do.call(rbind, results_list)
cat("Total annotated SNPs:", nrow(annotated_snps), "\n")

# Merge back with original df to retain Z-scores and P-values
final_annotated <- merge(snp_df, annotated_snps, by.x = "SNP", by.y = "refsnp_id", all.x = TRUE)

# Save result
write.csv(final_annotated, file = "C:\\Users\\NAILA\\OneDrive\\Documents\\R prac\\GWAS_project\\annotated_snps_biomart.csv", row.names = FALSE)
cat("Annotation saved to annotated_snps_biomart.csv\n")
