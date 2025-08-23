# Breast Cancer Early Onset & Related Traits in UK Populations

This project investigates the **genetic architecture of breast cancer** in the UK population, leveraging **GWAS summary statistics** and **GTEx eQTL data**. The pipeline integrates data preprocessing, harmonisation, quality control, genomic control correction, and colocalisation analysis, setting a foundation for fine-mapping, Mendelian randomization, polygenic risk scoring, and functional annotation.

## Table of Contents
- [Overview & Vision](#overview--vision)
- [Project Goals](#project-goals)
- [Data Sources](#data-sources)
- [Project Workflow](#project-workflow)
- [Execution Commands](#execution-commands)
- [Project Structure](#project-structure)
- [Results Summary](#results-summary)
- [Future Insights & Extensions](#future-insights--extensions)
- [License](#license)

## Overview & Vision

This project envisions a **comprehensive computational pipeline** to explore **UK-specific breast cancer genetics** and related traits such as age at menarche and mammographic density.

Key aspects of the vision:
- Integrate **large-scale GWAS datasets** with **functional genomics**
- Harmonise genetic data with **tissue-specific eQTLs** for mechanistic insights
- Identify **candidate genes** for further study
- Create a **reproducible, scalable framework** for population-specific genomic analysis

## Project Goals

1. **Identify genetic variants** associated with breast cancer risk in British ancestry individuals
2. **Preprocess and QC GWAS summary statistics** to ensure accuracy and consistency
3. **Harmonise GWAS and eQTL data** to prepare for colocalisation analysis
4. **Apply genomic control** to correct inflation (λ)
5. **Perform colocalisation analysis** to find shared causal variants between GWAS signals and gene expression
6. Establish a **pipeline framework** for downstream analyses like fine-mapping, Mendelian Randomization (MR), and Polygenic Risk Scores (PRS)

## Data Sources

### GWAS Summary Statistics
- **Study accessions:** GCST90435602, GCST90435601, GCST90435600
- **Traits:** Breast cancer and related PheCodes (174, 174.1, 174.11)
- **Discovery samples:** 11,874-12,898 cases; 388,549 controls of British ancestry
- **Source:** [GWAS Catalog](https://www.ebi.ac.uk/gwas/)

> ⚠️ Due to large file sizes, summary statistics are **not uploaded in this repo**. Links are provided in `.gitattributes`.

### eQTL Data
- **GTEx v8 cis-eQTLs**
- Breast Mammary Tissue:
  - `Breast_Mammary_Tissue.v8.egenes.txt`
  - `Breast_Mammary_Tissue.v8.eqtl_signifpairs.txt`

> ⚠️ GTEx files are referenced in `.gitattributes` due to size.

## Project Workflow

### Stage 0: Setup & Planning
- Initialize Linux environment with **Conda or Docker** for reproducibility
- Initialize **Git repository** for version control and tracking progress

### Stage 1: Data Acquisition
- Download GWAS summary statistics from GWAS Catalog
- Download GTEx v8 eQTL data

### Stage 2: Data Preprocessing
- Clean GWAS summary statistics: harmonise columns, remove low-quality variants, ambiguous alleles, and very low MAF SNPs
- Harmonise with eQTL: match variant positions, align effect alleles
- Purpose: ensure **high-quality input** for downstream analyses

### Stage 3: Quality Control & Visualization
- Generate **Manhattan and QQ plots** to inspect GWAS signals and inflation
- Detect and correct **genomic inflation** (λ >1) via genomic control
- Visualize effect size distributions

### Stage 4: Post-GWAS Colocalisation
- Merge harmonised GWAS and eQTL datasets
- Run **coloc analysis** per gene to detect **shared causal variants**
- Outputs:
  - Gene-level summary (`coloc_summary.tsv`)
  - SNP-level coloc results (`coloc_<gene>_snp.tsv`)

### Future Stages
- **Fine-mapping:** narrow credible variant sets with susieR or FINEMAP
- **Mendelian Randomization (MR):** test causality between traits and breast cancer
- **Polygenic Risk Score (PRS):** build PRS models and evaluate risk prediction
- **Functional annotation & integrative genomics:** use tools like CADD, RegulomeDB
- **Cross-population analyses & clinical translation**

## Execution Commands

### Preprocess GWAS Summary Statistics
`python scripts/1_preprocess.py`

#### Generated Files
- `cleaned_gwas.tsv` - Created by running: `python scripts/1_preprocess.py`
> ⚠️ .tsv files are referenced in `.gitattributes` due to size.

### Harmonise GWAS and eQTL Data
`python scripts/2_harmonization.py`

### Genomic Control Correction
`Rscript scripts/4_genomic_control_correction.R`

### Generate GWAS Plots
`Rscript scripts/3_generate_plots.R`

### Bash commands to fix harmonized GC file
```bash
# Remove quotes
sed 's/"//g' results/harmonized_GC.tsv > results/harmonized_GC_nq.tsv

# Fix variant_id and filter MAF
cat results/harmonized_GC_nq.tsv \
  | sed 's/"//g' \
  | awk -F'\t' 'NR==1{print; next} {split($9,a," "); $9=a[1]; print}' OFS='\t' \
  | awk -F'\t' 'NR==1 || ($9>0 && $9<1)' \
  > results/harmonized_GC_fixed2.tsv
```

### Run Colocalisation
`Rscript scripts/colocalisation.R`



## Project Structure

```
trait-choice-breast-cancer/
│── data/                 # Raw GWAS + GTEx files (links in .gitattributes)
│── results/
│   ├── cleaned_gwas.tsv
│   ├── harmonized.tsv
│   ├── harmonized_GC.tsv
│   ├── plots/
│   └── coloc/
│── scripts/
│   ├── 1_preprocess.py
│   ├── 2_harmonization.py
│   ├── 3_generate_plots.R
│   ├── 4_genomic_control_correction.R
│   └── colocalisation.R
│── README.md
│── .gitattributes
```

## Results Summary

- After preprocessing and harmonisation, ~260 SNPs remained for analysis
- Genomic inflation λ = 1.375, corrected with genomic control
- Colocalisation analysis: top PP.H4 ≈ 0.699; suggestive, but below strict 0.7 threshold
- Plots generated for visual QC: Manhattan, QQ, and effect size distribution

## Future Insights & Extensions

- **Fine-mapping:** Identify credible causal variants using susieR / FINEMAP
- **Mendelian Randomization:** Assess causality of risk factors (e.g., age at menarche) on breast cancer
- **Polygenic Risk Scores (PRS):** Build and validate UK-specific PRS models
- **Multi-omics integration:** Include methylation, proteomics, and splicing QTLs
- **Functional annotation:** Prioritize variants with CADD, RegulomeDB, and ENCODE data
- **Cross-population comparison:** Compare UK-specific results with other ancestries
- **Clinical translation:** Develop PRS thresholds for risk stratification in healthcare
- **Dashboard & pipeline packaging:** Interactive visualization using R Shiny or Python Streamlit

## License

This project is licensed under the MIT License - see the LICENSE file for details.
