# Integrative RNA-Seq and Co-expression Analyses Uncover Confident sRNA

## Overview

This repository provides a complete computational framework to identify and prioritize biologically consistent sRNA–mRNA regulatory interactions by integrating:

- Genomic position classification  
- RNA-seq differential expression analysis  
- Cross-strain regulatory consistency  
- Co-expression network analysis (WGCNA)  
- Multi-tool interaction prediction filtering  

The workflow combines sequence-based predictions with transcriptomic and network-level evidence to produce high-confidence regulatory interactions.

<img width="913" height="440" alt="image" src="https://github.com/user-attachments/assets/fec07b3b-7940-4e5a-b931-9b27616445b3" />
---

# Repository Structure

```
.
├── 00_Test_data
├── 01_Position_Classification/
├── 02_Co-expression_analysis/
├── 03_Filters/
└── README.md
```

The modules must be executed in order.

---

# Required External Tools

## Interaction Prediction

- **IntaRNA**  
  https://github.com/BackofenLab/IntaRNA  

- **RNAplex (ViennaRNA Package)**  
  https://www.tbi.univie.ac.at/RNA/  

- **TargetRNA3**  
  https://cs.wellesley.edu/~btjaden/TargetRNA3 

- **sRNARFTarget**  
  https://github.com/BioinformaticsLabAtMUN/sRNARFTarget

---

## RNA-seq Processing

- **SRA Toolkit**  
  https://github.com/ncbi/sra-tools  

- **Nextflow**  
  https://www.nextflow.io/  

- **nf-core/rnaseq**  
  https://github.com/nf-core/rnaseq  

- **nf-core/differentialabundance**  
  https://github.com/nf-core/differentialabundance  

---

## Network and Functional Analysis

- **WGCNA (R package)**  
  https://doi.org/10.1186/1471-2105-9-559

---


# Software Requirements

## Python ≥ 3.9

```
pip install pandas numpy
```

## R ≥ 4.1

```
install.packages(c(
  "tidyverse",
  "WGCNA",
  "flashClust",
  "fastDummies"
))
```

---

# Recommended Execution Order

1. Run **01_Position_Classification**
2. Run RNA-seq processing and DESeq2
3. Run **02_Co-expression_analysis**
4. Run **03_Filters**

---

# Citation

If this repository is used in research, please cite:

- IntaRNA  
- ViennaRNA / RNAplex  
- TargetRNA3  
- sRNARFTarget  
- WGCNA  
- nf-core RNA-seq  

---

# Contact

For issues or reproducibility questions, please open a GitHub issue in this repository.
