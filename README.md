# Integrative RNA-Seq and Co-expression Analyses Uncover Confident sRNA

## Overview

This repository provides a framework to identify and prioritize biologically consistent sRNA–mRNA regulatory interactions by integrating:

1. Run **01 Position Classification**
2. Run RNA-seq processing (AUXILIARY COMMANDS) **nf-core/rnaseq**  
3. Run **02 Co-expression analysis**
4. Run DESeq2 (AUXILIARY COMMANDS) **nf-core/differentialabundance**
5. Run sRNA-mRNA Target prediction programs
6. Run **03 Filters**

The workflow combines sequence-based predictions with transcriptomic and network-level evidence to produce high-confidence regulatory interactions.
A more complete description of each of the modes and their use is available in the readme file within each module.

<img width="913" height="440" alt="image" src="https://github.com/user-attachments/assets/fec07b3b-7940-4e5a-b931-9b27616445b3" />
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


# Contact

For issues or reproducibility questions, please open a GitHub issue in this repository.
