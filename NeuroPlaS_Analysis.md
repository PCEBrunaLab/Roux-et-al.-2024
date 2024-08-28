
# Neuroblastoma Plasticity Signature (NeuroPlaS) Analysis

## Introduction

The Neuroblastoma Plasticity Signature (NeuroPlaS) was developed to identify plasticity states in neuroblastoma cells using single-cell RNA sequencing data. This method utilizes differential expression analysis to identify genes that are upregulated in INT cells compared to ADRN and MES cells. Additionally, the robustness of NeuroPlaS is validated across various datasets.
## Methods

### NeuroPlaS Identification and validation

1. **Differential Expression Analysis**: To identify genes that are upregulated in INT cells compared to ADRN and MES cells, differential expression analysis was conducted using the `FindMarkers` function from the Seurat package.
   - Genes were considered differentially expressed if the adjusted p-value was less than 0.05 and log fold change (logFC) > 1.

2. **Validation Across Datasets**: The robustness of NeuroPlaS was validated across different datasets, including:
   - Patient-derived xenografts (PDX)
   - Patient-derived organoids (PDO)
   - Cell line samples
   - Patient samples from public datasets (Great Ormand Street Hospital with 5 patients and Princess Maxima Center with 16 patients)

The expression of NeuroPlaS genes was assessed in cells annotated as ADRN, MES, and SYM using the `AddModuleScore` function from the Seurat package. This function calculates a module score for each cell based on the expression levels of NeuroPlaS genes.

3. **NuroPlaS in DTPs**: Cells are grouped into "Survivor" and "Non-Survivor" categories based on whether their barcodes match those of known Cisplatin survivors and calculates a module score for each cell based on the expression levels of NeuroPlaS genes. 

### Survival Analysis

1. **Enrichment Score Calculation**: In a large bulk RNA-sequencing dataset comprising 498 primary neuroblastoma patient samples (GSE62564/GSE49711), the enrichment score of NeuroPlaS was calculated using Gene Set Variation Analysis (GSVA) with the single-sample Gene Set Enrichment Analysis (ssGSEA).

2. **Kaplan-Meier Survival Analysis**: Based on the enrichment score, samples were categorized into “high” and “low” NeuroPlaS groups. Kaplan-Meier survival analysis was conducted, and the log-rank test was used to determine the correlation between high NeuroPlaS scores and overall survival (OS).

3. **Stratification by MYCN Status**: To assess the interaction between MYCN status and NeuroPlaS, the analysis was stratified by MYCN amplification status (amplified vs. non-amplified). Patients were grouped based on MYCN amplification status and their NeuroPlaS classification (high vs. low).



