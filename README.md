# In Vitro Germination of Polyembryonic Mango (*Mangifera indica* L.) cv. 'Chulucanas'

## Overview

This repository contains the R code used for the statistical analyses and graphical outputs of the study entitled:

**"In Vitro Germination of Polyembryonic Mango (*Mangifera indica* L.) cv. 'Chulucanas' Seeds in Liquid Medium Supplemented with Coconut Water and Lemon Juice"**

The study evaluated the effects of fruit maturity stage, coconut water, and lemon juice on the in vitro germination and establishment of polyembryonic mango seeds cv. 'Chulucanas'.

## Experimental Design

The experiment was conducted using a completely randomized design (CRD) with a **2 × 2 × 3 factorial arrangement**:

- **Maturity stage (MS):** physiological maturity (PM) and commercial maturity (CM)
- **Coconut water (CW):** 0 and 20% (v/v)
- **Lemon juice (LJ):** 0, 1.5, and 3.0% (v/v)

The factorial combination resulted in **12 treatments (T1–T12)**. Each treatment consisted of four replicates with five seeds per replicate, totaling **20 seeds per treatment and 240 seeds in the experiment**.

## Response Variables

The statistical analyses include fruit and seed characterization and the following in vitro responses:

- Developed embryos (DE)
- Developed plumules (DP)
- Leaves of developed seedlings (LDS)
- Phenolic exudation (PE)
- Microbial contamination (MC)
- In vitro establishment (IVE)

Fruit characterization variables include fruit length, fruit width, fruit weight, total soluble solids, fruit firmness, and seed weight without endocarp.

## Statistical Analysis

Statistical analyses were performed in **R**.

Developed embryos were analyzed using a **Conway–Maxwell–Poisson (COM-Poisson) generalized linear model with a log link**, whereas developed plumules were analyzed using a **Poisson model with a log link**.

Phenolic exudation, microbial contamination, and in vitro establishment were analyzed using **binomial generalized linear models with a logit link and bias-reduced estimation**.

The effects of maturity stage (MS), coconut water (CW), lemon juice (LJ), and their interactions were evaluated using **Wald chi-square tests**.

When appropriate, estimated marginal means and Tukey-adjusted multiple comparisons were obtained using the `emmeans` package.

The number of leaves was analyzed as a conditional response for seeds that developed shoots.

Multivariate analyses included:

- Spearman correlation analysis
- Principal component analysis (PCA)
- Hierarchical clustering and heatmap visualization based on standardized Z-scores

## Repository Contents

### `Supplementary_Material_S1_Statistical_Analysis.R`

Complete R script used to:

- Import and prepare the experimental data
- Generate descriptive statistics
- Fit statistical models
- Perform model diagnostics
- Conduct Wald chi-square tests
- Estimate marginal means and multiple comparisons
- Generate Tables 1–3
- Generate Figures 1–5
- Perform Spearman correlation analysis
- Perform principal component analysis (PCA)
- Generate the Z-score heatmap
- Export statistical results and graphical outputs

## Main R Packages

The analyses use the following R packages:

`tidyverse`, `googlesheets4`, `janitor`, `car`, `glmmTMB`, `DHARMa`, `emmeans`, `brglm2`, `multcomp`, `multcompView`, `Hmisc`, `FactoMineR`, `ComplexHeatmap`, `circlize`, `cowplot`, `scales`, `writexl`, and related dependencies.

## Reproducibility

To reproduce the analyses:

1. Download or clone this repository.
2. Open `Supplementary_Material_S1_Statistical_Analysis.R` in R or RStudio.
3. Install the required R packages if they are not already installed.
4. Run the script sequentially from beginning to end.
5. Statistical tables, diagnostic outputs, and figures will be generated automatically.

A random seed is specified in the script to improve computational reproducibility.

## Output

The R script generates the statistical results associated with the manuscript, including:

- **Table 1:** Fruit and seed characterization according to maturity stage
- **Table 2:** Effects of maturity stage, coconut water, lemon juice, and their interactions on in vitro responses
- **Table 3:** Model-estimated probabilities for phenolic exudation, microbial contamination, and in vitro establishment
- **Figure 1:** Phenolic exudation, microbial contamination, and in vitro establishment
- **Figure 2:** Developed embryos, developed plumules, and leaves of developed seedlings
- **Figure 3:** Spearman correlation matrices by fruit maturity stage
- **Figure 4:** Principal component analysis (PCA)
- **Figure 5:** Heatmap based on standardized Z-scores

## Authors

**Gabriela Cárdenas-Huamán, Julieta Llihua-Quispe, Max Ramírez-Rojas, and Henry Morocho-Romero**

## Data Availability

The R code required to reproduce the statistical analyses is publicly available in this GitHub repository.

## Citation

If you use the code or data from this repository, please cite the associated scientific article. Full bibliographic information will be added after publication.

## Contact

For questions regarding the statistical analyses or repository, please contact the corresponding authors through the contact information provided in the associated publication.
