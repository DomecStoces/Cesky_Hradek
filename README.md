# Carabid Beetle Assembly Along an Elevational Gradient: Analytical Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19444351.svg)](https://doi.org/10.5281/zenodo.19444351)

This repository contains the R scripts and datasets used to analyze carabid beetle assemblages across a continuous submontane-to-montane elevational gradient (390–820 m a.s.l.) in the Czech Republic. The workflow follows the exact methodology described in Section 2.3 of the associated manuscript.

## 0. Preliminary Data Processing
Before running the main analyses, the raw field data was processed and formatted.
* **Scripts used:** `Birch.R` and `Conversion_dataset.R`
* **Description:** These scripts were used for the preliminary processing of raw trap datasets, transforming them into the finalized `df` (dataframe) and `data_long` formats required for downstream community and trait-based analyses.

---

## Step 1: Ordination of Species Composition and Multicollinearity Check (Method 2.3.1)
To analyze shifts in carabid assemblage composition, we utilized multivariate ordination. Rare species were down-weighted to prevent artificial stretching of the chi-squared distance metric.
* **Tool used:** CANOCO 5 (for partial Detrended Canonical Correspondence Analysis - pDCCA)
* **Dataset used:** `FINAL_CANOCO.xlsx` (specifically the `env` and `sp` sheets).
* **Script used:** `Multicol.R`
* **Description:** Before performing the pDCCA in CANOCO, we validated the inclusion of environmental covariates (Temperature, Precipitation, Wind) by testing for multicollinearity. `Multicol.R` uses the `corrplot` and `dplyr` packages to compute Spearman rank correlations and visualize them. It also calculates Variance Inflation Factors (VIF) using the `car` package to ensure values remain below the threshold of 5.

---

## Step 2: Calculation of Beta-Diversity Components (Method 2.3.2)
We adopted the framework proposed by Šizling et al. (2026) to evaluate complementary, non-additive indices of community variation.
* **Script used:** `Beta_diversity.R`
* **Dataset used:** `data_long.xlsx` (aggregated by Locality, Elevation, Exposition, Year, Month, and Species).
* **Description:** 1. Computes three presence-absence dissimilarity matrices using `vegan::designdist`: Jaccard (overall dissimilarity), Simpson (species turnover), and a richness difference index (nestedness/uniformity).
  2. Matrices are square-root transformed to ensure Euclidean properties.
  3. Tests multivariate homogeneity of group dispersions using `betadisper` (PERMDISP).
  4. Runs Permutational Multivariate Analysis of Variance (PERMANOVA) using `adonis2`. Permutations (n=999) are strictly constrained within individual sampling sites (`strata = df_agg$Locality`) to account for the repeated-measures design.

---

## Step 3: Trait-based and Functional Diversity Indices (Method 2.3.3)
We evaluated how key niche traits (trophic strategy, body size, dispersal ability) respond to the elevational gradient.
* **Scripts used:** `FD_RaoQ.R` and `CWM.R` (for PCA and collinearity of traits)
* **Dataset used:** `df.csv` and `Rao_diversity.xlsx`
* **Description:** * Ordinal traits (trophic strategy, dispersal ability) were rescaled to a 0–1 interval.
  * `FD_RaoQ.R` computes the multidimensional functional diversity of the carabid assemblages using Rao's quadratic entropy (Q). Because the traits mix continuous and categorical variables, Gower's distance is calculated first (via the `cluster` package), followed by Rao's Q computation (via the `fundiversity` package).
  * `CWM.R` calculates Community-Weighted Means and visually checks for trait redundancy via PCA (`FactoMineR`) and Spearman correlations.

---

## Step 4 & 5: Phylogenetic Diversity and GAM Modelling (Methods 2.3.4 & 2.3.5)
We assessed evolutionary assembly (Mean Phylogenetic Diversity and Standardized Effect Size) and modelled community metrics against the elevational gradient.
* **Script used:** `Phylo.R`
* **Datasets used:** `PD.xlsx` (for `meanPD`) and `df.csv` (for `SESpd`).
* **Description:** * Models community metrics using Generalized Additive Models (GAMs) via the `mgcv` package. 
  * `meanPD` and `SESpd` are modelled with a Gaussian error distribution, `s(Altitude_scaled)` as a cubic regression spline, and `Locality` and `Year` as random-effect smoothers (`bs = "re"`).
  * Evaluates model fit and assumptions using `gam.check()` and concurvity diagnostics.
  * **Spatial Autocorrelation:** The script explicitly tests for residual spatial dependencies. It extracts Pearson residuals from the GAMs, calculates Moran’s I (`spdep`), and constructs empirical variograms (`gstat`) to confirm a pure nugget effect (absence of spatial autocorrelation).
  * Finally, it uses `ggplot2` to visualize the GAM trendline and significance thresholds ($\pm 1.96$) for phylogenetic clustering and overdispersion.

---

## Dependencies
Ensure you have **R (version 4.5.2+)** and the following core packages installed:
`vegan`, `mgcv`, `mgcViz`, `qgam`, `DHARMa`, `spdep`, `gstat`, `sp`, `picante`, `cluster`, `fundiversity`, `FactoMineR`, `corrplot`, `car`, `dplyr`, `tidyr`, `ggplot2`, `readxl`.
