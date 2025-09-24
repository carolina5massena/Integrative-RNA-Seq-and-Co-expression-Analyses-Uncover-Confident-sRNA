# Integrative RNA-Seq and Co-expression Analyses Uncover Confident sRNA

A notebook-driven pipeline to identify and prioritize high-confidence regulatory sRNAs by integrating RNA-seq expression, genomic position categories, target-filtering logic, and co-expression network (WGCNA) analyses.

This repository accompanies the article “Integrative RNA-Seq and Co-expression Analyses Uncover Confident sRNA Networks in Staphylococcus aureus Strains of the Brazilian ST239 and USA Lineages”.

**Important scope note**  
This repository contains **three Jupyter notebooks that implement the created code** for the core integration steps. Some analyses must be performed **prior to running these notebooks** and are **not implemented here**—specifically:  
> 1) running one or more **sRNA→mRNA prediction programs**, and  
> 2) performing **differential expression (DEG)** analyses for sRNAs and mRNAs.  
> The notebooks assume you have already generated these prerequisite results and saved them as tabular files.

## Contents

### Position_Categorization.ipynb
Categorizes sRNAs by genomic context (e.g., intergenic, antisense, UTR-derived, intragenic) via overlap with genome annotations.

<img width="678" height="615" alt="image" src="https://github.com/user-attachments/assets/dbc7baf3-39d3-4336-bcd1-ba38865d0308" />

### WGCNA.ipynb
Builds co-expression networks (WGCNA), identifies trait-associated modules, and extracts hub genes/sRNAs. Exports figures and Cytoscape-ready networks.

<img width="459" height="244" alt="image" src="https://github.com/user-attachments/assets/94741eca-afff-47cd-b514-08dd81767ac3" />

### Filtros_sRNA_mRNAtarget.ipynb
Applies biologically motivated filters to select sRNA→mRNA pairs by combining prediction resources and expression evidence.

<img width="913" height="440" alt="image" src="https://github.com/user-attachments/assets/fec07b3b-7940-4e5a-b931-9b27616445b3" />


(Optional) Target predictions

data/targets/predictions.tsv with candidate sRNA→mRNA pairs and scores
