# MAL-ED Anemia and Nutritional Biomarker Analysis

R code for the anemia and nutritional biomarker analysis using the MAL-ED
(Malnutrition and Enteric Disease) longitudinal cohort dataset. The script
classifies anemia using WHO 2024 age/sex-specific hemoglobin cut-offs,
compares these to the previous WHO cut-offs, and analyzes associations
between hemoglobin status and iron, vitamin B12, and lead biomarkers across
childhood.

## What this script does

- Classifies hemoglobin values as anemic/non-anemic and by severity (Mild /
  Moderate / Severe) using WHO 2024 age- and sex-specific cut-offs
- Compares anemia prevalence under the previous (2011) vs. current WHO
  cut-offs at 7 and 15 months
- Fits GEE models (exchangeable correlation) to estimate the linear and
  non-linear trend in anemia prevalence with age
- Classifies iron sufficiency (body iron), vitamin B12 status, and lead
  exposure, and cross-tabulates these against anemia status
- Computes Spearman correlations between hemoglobin and each biomarker by
  age, with Benjamini-Hochberg-adjusted p-values
- Assesses data completeness and attrition (retained vs. dropped
  participants) across follow-up visits
- Compares height-for-age z-scores (HAZ) by anemia status
- Produces publication-ready tables (as R tibbles/data frames, intended for
  export via flextable/officer) and figures (ggplot2/patchwork, exported as
  TIFF in the original analysis pipeline)

## Input data

The script expects an Excel file named `maled_data.xlsx` in the working
directory, with one row per participant. Required columns (raw name ->
renamed in script):

| Raw column | Renamed to | Description |
|---|---|---|
| `si_no` | `pid` | Participant ID |
| `sex` | `CAFSEX` | Sex (coded 1 = Boys, 2 = Girls) |
| `Hb_7`, `Hb_15`, `Hb_24`, `Hb_36`, `Hb_60`, `Hb_84`, `Hb_108`, `Hb_144` | `Hb7`...`Hb144` | Hemoglobin (g/dL) at each visit (age in months) |
| `Ferritin_7`, `Ferritin_15`, `Ferritin_24`, `Ferritin_108`, `Ferritin_144` | `Ferritin7`... | Ferritin |
| `Tfr_7`, `Tfr_15`, `Tfr_24`, `Tfr_108`, `Tfr_144` | `tfr7`... | Transferrin receptor |
| `Bodyiron_7`, `Bodyiron_15`, `Bodyiron_24`, `Bodyiron_108`, `Bodyiron_144` | `bodyiron7`... | Body iron (mg/kg) |
| `B12_108`, `B12_144` | `B12_9y`, `B12_12y` | Vitamin B12 (pg/mL) at 9y and 12y |
| `Lead_15`, `Lead_24`, `Lead_144` | `Lead_15`, `Lead_24`, `Lead_12y` | Blood lead (µg/dL) |
| `HAZ_7`, `HAZ_15`, `HAZ_24`, `HAZ_108`, `HAZ_144` | (unchanged) | Height-for-age z-scores |

**Note:** The MAL-ED dataset itself is not included in this repository, as
it is subject to a data use agreement. It is available upon request from
the MAL-ED Network / study data repository. This repository shares only
the analysis code.

## Requirements

R (version used for analysis: see `sessionInfo.txt` in this repo) with the
following packages:

- dplyr
- tidyr
- ggplot2
- ggrepel
- patchwork
- readxl
- geepack
- splines (base R)

Install with:

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "ggrepel", "patchwork", "readxl", "geepack"))
```

## Usage

1. Place `maled_data.xlsx` (with the columns listed above) in the working
   directory.
2. Run the script (e.g. `anemia_biomarker_analysis.R`) in R or RStudio.
3. Outputs are printed to the console as tables (Table 1, GEE trend table,
   Table 2, Supplementary Tables 1a/1b/2/3/4) and plotted as figures
   (Figure 2: anemia prevalence by age; Figure 3: hemoglobin-biomarker
   correlation heatmap).

## Output overview

| Output | Description |
|---|---|
| `table1` | Anemia prevalence and severity by age and sex, with Wilson 95% CI |
| `gee_table` | GEE-estimated trend in anemia prevalence with age |
| `fig2` | Anemia prevalence by age, overall and by sex |
| `table2` | Iron/B12/lead prevalence by anemia status, with Wilson 95% CI |
| `supp_table1a` | Data completeness by visit |
| `supp_table1b` | Attrition sensitivity check (retained vs. dropped) |
| `supp_table2` | Iron/B12/lead prevalence by age and sex |
| `supp_table3` | Combined iron/B12/lead/anemia status at 144 months |
| `fig3` | Spearman correlation heatmap: hemoglobin vs. biomarkers by age |
| `supp_table4` | Height-for-age z-score by anemia status |

## License

MIT License — see `LICENSE`.

## Citation

If you use this code, please cite the associated manuscript (details to be
added upon publication).
