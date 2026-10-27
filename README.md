# Interactive Effects of Fruit Maturity Conditions and Organic Supplements on Early *In Vitro* Establishment of Polyembryonic Mango (*Mangifera indica* L.) cv. ‘Chulucanas’

## Overview

This repository provides the experimental dataset and R scripts associated with the revised manuscript submitted to *Plants* (MDPI; previous manuscript ID: plants-4614048).

The study investigates the interactive effects of fruit maturity/postharvest conditions, coconut water (CW), and key lime juice (LJ) on early *in vitro* establishment and related developmental responses of polyembryonic mango cv. ‘Chulucanas’. The revised statistical workflow was developed following peer review. The earlier analysis script is retained for transparency and version tracking.

## Experimental Design

A completely randomized design (CRD) was used with a **2 × 2 × 3 factorial arrangement**:

- **Fruit maturity/postharvest condition (MS):** physiological maturity (PM; seeds extracted immediately after harvest) and commercial maturity (CM; fruits ripened for 10 days under ambient conditions before seed extraction).
- **Coconut water (CW):** 0 and 20% (v/v).
- **Key lime juice (LJ):** 0, 1.5, and 3.0% (v/v).

The design comprised **12 treatment combinations (T1–T12)**, with **four independently prepared medium batches per treatment** and **five seeds per batch**. Each seed was cultured in a separate vessel. Thus, the experiment included **48 independent experimental replicates and 240 seed-level observations**. The independently prepared medium batch, not the individual seed, was the experimental unit for statistical inference.

All fruits originated from a single mother tree. Because the CM condition included a 10-day postharvest ripening period, maturity condition and postharvest handling cannot be disentangled. The findings are limited to the cultivar, experimental material, tested concentrations, and four-week culture period.

## Response Variables

Fruit and seed characterization included fruit length, fruit width, fruit weight, fruit firmness, total soluble solids, and seed weight without the endocarp. These characteristics were summarized descriptively; no inferential comparisons between maturity conditions were performed.

The *in vitro* responses were:

- **DE:** number of developed embryos per seed (two weeks).
- **PE:** visually assessed brown-to-dark tissue or medium discoloration, recorded as visible exudation (two weeks).
- **MC:** visible microbial contamination (two weeks).
- **DP:** number of developed plumules per seed (four weeks).
- **LDS:** number of leaves per developed seedling (four weeks), assessed conditionally where plumules developed.
- **IVE:** early *in vitro* establishment, assessed by visible continued embryo growth at four weeks, independently of plumule or leaf formation.

Visible exudation was not chemically identified or quantified as phenolic compounds. The four-week endpoint does not establish subsequent plantlet survival or successful acclimatization.

## Statistical Analysis

Analyses were performed in **R**, using the **48 replicate-level observations** to avoid pseudoreplication. The five seed-level observations within each replicate were aggregated as appropriate for each response.

- **DE:** Conway–Maxwell–Poisson (COM-Poisson) regression with a log link.
- **DP:** Poisson regression with a log link.
- **PE, MC, and IVE:** grouped binomial generalized linear models with a logit link and bias-reduced estimation.
- **LDS:** conditional analysis restricted to replicates with developed plumules (**23 of 48**). Main effects were assessed without estimating factorial interactions because of the limited evaluable sample; LDS was not estimable for T9.

The effects of MS, CW, LJ, and their interactions were assessed using Wald chi-square tests for the count and grouped-binomial responses. Where appropriate, estimated marginal means and Tukey-adjusted pairwise comparisons were calculated using `emmeans`. Model adequacy was assessed using appropriate residual and dispersion diagnostics, including DHARMa simulation-based diagnostics for count models.

Exploratory multivariate analyses included Spearman rank correlations, principal component analysis (PCA), and hierarchical clustering with a heatmap of standardized treatment-level means (Z-scores). These analyses were used to describe patterns and associations, not to establish causality.

## Repository Contents

### `Mango_Chulucanas_Statistical_Analysis_Revised.R`

Revised statistical analysis script associated with the resubmitted manuscript. It contains the updated model specifications, replicate-level analyses, diagnostics, and multivariate workflows. **Use the CSV-compatible version of the revised script if the script above does not directly import the supplied CSV file.**

### `Mango_Chulucanas_Statistical_Analysis_Revised_CSV.R`

CSV-compatible version of the revised analysis workflow. It reads `mango_chulucanas_raw_data.csv` and exports statistical results and graphical outputs. This file should be used as the primary reproducibility script **if it is included in the repository**.

### `Supplementary_Material_S1_Statistical_Analysis.R`

Original statistical analysis script retained for transparency and version tracking. It does not supersede the revised workflow.

### `mango_chulucanas_raw_data.csv`

Experimental dataset comprising **240 seed-level observations**, grouped into **48 independently prepared medium-batch replicates** across the 12 treatment combinations.

## Main R Packages

The workflow uses R packages including `tidyverse`, `janitor`, `car`, `glmmTMB`, `DHARMa`, `emmeans`, `brglm2`, `multcomp`, `multcompView`, `Hmisc`, `FactoMineR`, `ComplexHeatmap`, `circlize`, `cowplot`, `scales`, and `writexl`.

## Reproducibility

1. Download or clone the repository.
2. Place `mango_chulucanas_raw_data.csv` in the same working directory as the CSV-compatible revised R script.
3. Open `Mango_Chulucanas_Statistical_Analysis_Revised_CSV.R` in R or RStudio, if that file is present in the repository.
4. Install the required R packages.
5. Set the working directory to the repository root and run the script sequentially.
6. Inspect the generated tables, diagnostics, and figures, and compare them with the final manuscript.

The revised CSV-compatible script specifies `set.seed(123)` for reproducible simulation-based procedures. It creates a `Results_FINAL` directory for exported outputs.

**Important:** The current CSV-compatible script retains an earlier numbering scheme for its exported figure filenames. These filenames must not be interpreted as the final figure numbers in the revised manuscript. See the mapping below.

## Tables and Figures in the Revised Manuscript

### Tables

- **Table 1:** Physical and physicochemical characteristics of mango fruits and seeds under the two maturity/postharvest conditions.
- **Table 2:** Effects of MS, CW, LJ, and their interactions on *in vitro* responses; includes a separate conditional analysis for LDS.
- **Table 3:** Model-estimated probabilities of visible exudation, microbial contamination, and early *in vitro* establishment across treatment combinations.
- **Table 4:** Treatment combinations in the factorial arrangement and mean initial culture-medium pH before adjustment and sterilization.
- **Supplementary Table S1:** Number of evaluable replicates for the conditional LDS analysis.

### Figures

- **Figure 1:** Representative photographs of morphological responses during *in vitro* culture, including visible discoloration, microbial contamination, plumule development, and leaf development.
- **Figure 2:** Observed means for (A) developed embryos, (B) developed plumules, and (C) leaves per developed seedling.
- **Figure 3:** Model-estimated probabilities and 95% confidence intervals for (A) visible exudation, (B) microbial contamination, and (C) early *in vitro* establishment.
- **Figure 4:** Spearman rank correlation matrices for fruit and seed characteristics and *in vitro* responses.
- **Figure 5:** PCA biplot of fruit and seed characteristics and *in vitro* responses.
- **Figure 6:** Hierarchical clustering heatmap of standardized treatment-level means (Z-scores).
- **Supplementary Figure S1:** DHARMa residual diagnostics for DE and DP count-response models.

### Existing Script Export Names Versus Manuscript Figure Numbers

| Current CSV-compatible script output | Corresponding revised manuscript figure |
| --- | --- |
| Figure 1: `Figure1_PE_MC_IVE.jpg` | **Figure 3** |
| Figure 2: `Figure2_DE_DP_LDS.jpg` | **Figure 2** |
| Figure 3: `Figure3_Spearman.jpg` | **Figure 4** |
| Figure 4: `Figure4_PCA.jpg` | **Figure 5** |
| Figure 5: `Figure5_Heatmap.jpg` | **Figure 6** |
| `Figure_S1_DHARMa.jpg` | **Supplementary Figure S1** |

**Figure 1 in the revised manuscript consists of representative experimental photographs and is not generated by the statistical script.** The script's exported filenames should be renamed or its export commands updated to match the manuscript before final archiving. Tables 4 and S1 should likewise be checked against the final manuscript and supplementary files; they are not necessarily generated as identically numbered worksheets by the script.

## Authors

Gabriela Cárdenas-Huamán, Julieta Llihua-Quispe, Max Ramírez-Rojas, and Henry Morocho-Romero.

## Data Availability

The experimental dataset and the revised R code are provided in this repository to support transparency and reproducibility. The original analysis script is retained for version tracking. The definitive dataset and executable script filenames should be confirmed against the files uploaded to GitHub before manuscript submission.

## Citation

If you use these data or scripts, please cite the associated scientific article. Complete bibliographic details will be added following publication.

## Contact

For questions regarding the dataset, analyses, or repository, please contact the corresponding author using the contact details provided in the associated manuscript or publication.
