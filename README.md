## Project Title

Assessing Immune Correlates Of Protection In Young Children Following Influenza Vaccination: Application Of Information-theoretic Approach

## Overview

This repository contains the main code used for the analysis of immune correlates of protection (CoP) in young children following influenza vaccination.

The analysis applies the information-theoretic framework proposed by Alonso et al. (2024) to evaluate the surrogacy value of hemagglutination inhibition (HI) antibody titers against H1N1 and H3N2 influenza infection endpoints. The methodology adapts functions from the `Surrogate` R package and extends them to the study-specific setting.

## Project Description

The primary objective of this project is to assess whether post-vaccination HI titers can serve as valid surrogate markers for protection against influenza infection.

Analyses include:

- H1N1 primary analysis
- H1N1 secondary analysis
- H3N2 analysis

---

## Reproducibility

### Software Requirements

- R version 4.5.3
- R package environment managed using `renv`

Package versions and dependencies are captured in:

- `renv.lock`
- `renv/`
- `.Rprofile`

### Installation

1. Install the `renv` package:

```r
install.packages("renv")
```

2. Ensure the following folder and files are present in the project directory:

```text
renv/
renv.lock
.Rprofile
```

3. Restore the project environment:

```r
renv::restore()
```

---


## Local project folder structure:

```text
.
├── Code/                 # Analysis scripts and functions
│   ├── CoP_analysis_H1N1_search_SV.qmd
│   ├── CoP_analysis_H3N2_search_SV.R
│   ├── CoP_summary.qmd
│   ├── Helper_functions.R
│   ├── ICA.BinCont.Cov_inner.R
│   ├── ICA.BinCont_full_loop.R
│   └── ICA_run_job_H1N1.R
├── Plot/                 # Generated figures and visualizations
├── Raw_data/             # Input/raw datasets
├── renv/                 # R package environment
├── Results/              # Analysis outputs
|   └── Summary/              # Summary tables and reports
├── .RData                # Saved R workspace
├── .Rhistory             # R command history
├── .Rprofile             # Project-specific R configuration
├── README.md             # Project documentation
└── renv.lock             # Reproducible package snapshot
```

---

## Analysis Workflow

1. `CoP_analysis_H1N1_search_SV.qmd` - Identify starting values for H1N1 model estimation.
2. `ICA_run_job_H1N1.R` - Perform primary and secondary H1N1 analyses.
3. `CoP_analysis_H3N2_search_SV.qmd` - Identify starting values for H3N2 models estimation. 
(Full primary and secondary analyses for H3N2 were not run because no models showed adequate goodness of fit).
4. `CoP_summy.qmd` - Generate summary statistics and plots for the H1N1 analyses.


## Authors

Trang Tu
[@trangtph](https://github.com/trangtph)
trang.tph97@gmail.com


## References
1. Alonso Abad A, Ong F, Stijven F, Van der Elst W, Molenberghs G, Van Keilegom I, et al. An information-theoretic approach for the assessment of a continuous outcome as a surrogate for a binary true endpoint based on causal inference: Application to vaccine evaluation. Statistics in Medicine. 2024;43(6):1083–102
2. [Introduction to renv](https://rstudio.github.io/renv/articles/renv.html#installing-packages)

