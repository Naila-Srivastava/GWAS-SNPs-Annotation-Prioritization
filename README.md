# GWAS SNPs Annotation, Prioritization and Interpretation associated with Asthma

## Overview

Genome-wide association studies (GWAS) have identified numerous single nucleotide polymorphisms (SNPs) associated with complex diseases such as asthma. However, not all SNPs are functionally relevant. This project implements a reproducible pipeline to filter, annotate, and prioritize asthma-associated SNPs using public GWAS datasets, statistical filtering, functional annotations, regulatory potential predictions, and visualization approaches.

The workflow integrates Python (Jupyter Notebooks), R (biomaRt and visualization), and online bioinformatics resources (NHGRI GWAS Catalog, FUMA and Pharos) to identify biologically meaningful variants. The ultimate goal is to highlight candidate SNPs and genes with potential roles in asthma pathogenesis, providing a foundation for downstream functional studies and personalized medicine approaches.

## Tools & Technologies

- **Languages:** Python (Jupyter Notebook) and R
- **Databases & Resources:** NHGRI GWAS Catalog, FUMA GWAS and Pharos 
- **Libraries:**
  Python: (`Pandas`, `NumPy`, `Matplotlib` and `Seaborn`)
  R: (`biomaRt`, `ggplot2`, `ggthemes` `readxl`, `tidyr`, `data.table`, `dplyr`, `stringr` and `igraph`)
- **Reproducibility:** RMarkdown for reporting and GitHub for version control

## Data Source

- NHGRI GWAS Catalog: GCST010042 (Han Y. et al.), containing asthma-associated SNPs and their metadata.

- Additional data integration from:
  - FUMA GWAS for regulatory annotation and deleteriousness prediction
  - Pharos for scoring and prioritising eQTL genes
 
## How to Run

1. Clone this repository: `git clone <repo_url>`
                          `cd asthma-gwas-snp-prioritization`
2. Install dependencies: `pip install -r requirements.txt`
3. Run analysis
4. Generate final report

## Methodology

## Features

- Automated SNP filtering by p-value and trait relevance.
- Functional annotation using Ensembl BioMart.
- Prediction of regulatory effects using FUMA.
- eQTL mapping for gene expression association.
- Scoring & prioritization of SNPs integrating multiple evidence sources.
- Multi-level visualization: Manhattan plots, graphs and networks.

## Visualizations

* Manhattan plot for genome-wide SNP significance
* Network graph showing top genes relationships
* Minor Allele Frequency Distribution graph
* SNPs Consequences barplot
* Top genes associated with asthma
* Genes associated with the most frequent KEGG Pathways

## Key Takeaways

- Successfully implemented a bioinformatics pipeline to prioritize asthma SNPs.
- Identified candidate SNPs with high regulatory potential and disease association.
- Integrated multiple datasets for a systems-level perspective of asthma genetics.
- Established a reproducible framework for GWAS SNP annotation and visualization.

## What’s Next?

- Extend prioritization to other immune-related traits (e.g., atopy, allergic rhinitis).
- Incorporate machine learning models for SNP classification.
- Expand to multi-omics integration (epigenetics, transcriptomics).
- Publish the pipeline as a ready-to-use workflow (Nextflow/Snakemake).

## References

1. NHGRI GWAS Catalog: https://www.ebi.ac.uk/gwas/ [Han Y. et al. (GCST010042)]
2. FUMA GWAS: https://fuma.ctglab.nl/
3. Pharos: https://pharos.nih.gov/

## Project Structure

```plaintext
Asthma pGWAS SNP Prioritization and Interpretation/  
│
├── README.md                       # You're reading this now  
├── .gitignore  
├── requirements.txt                # Python dependencies   
│
├── Python/                         # Jupyter notebook (Data cleaning & preprocessing)  
├── R/                              # R scripts  
├── report                          # RMarkdown report   
├── results/                        # Processed results, figures and plots  
└── data/                           # Input datasets  
