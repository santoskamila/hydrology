# Monthly precipitation gap filling — method comparison

R workflow to fill gaps in monthly rainfall series from a network of rain gauges and to compare the performance of nine gap-filling methods through cross-validation with artificial gaps. Originally developed in MATLAB and ported to R.

## Methods

| Code | Method | File |
|---|---|---|
| `RW` | Regional weighting (ratio of long-term annual means) | `R/fill_regional_weighting.R` |
| `SLR_Var` | Simple linear regression with variable donor (most correlated gauge) | `R/fill_slr_variable.R` |
| `RV` | Regional vector (regression on deviations in anomalous months) | `R/fill_regional_vector.R` |
| `PCR_NIPALS` | Principal component regression with NIPALS (handles missing data) | `R/fill_pcr_nipals.R` |
| `KNN_Analog` | k-nearest analogue months, inverse-distance weighted | `R/fill_knn_analog.R` |
| `KNN_Simple` | Mean of the k closest gauges (logic of R's `impute.knn`) | `R/fill_knn_simple.R` |
| `MLR_Fixed` | Multiple linear regression with the 5 most correlated donors | `R/fill_mlr_fixed.R` |
| `Arithmetic_Mean` | Arithmetic mean of the neighbouring gauges | `R/fill_arithmetic_mean.R` |
| `Iterative_PCA` | One-component iterative PCA on the whole matrix | `R/validation_helpers.R` |

Also included but not used in the default validation: `fill_linear_regression.R`, `fill_pcr_v2.R` (PCR with classical PCA) and `fill_multiple_imputation.R` (stochastic regression in log space).

Metrics: RMSE, MAE, R², Nash–Sutcliffe efficiency (NSE), Spearman correlation and Willmott's index of agreement (d).

## Structure

```
precipitation-gap-filling/
├── 01_prepare_data.R          # daily Excel -> monthly totals + long-term means
├── 02_validate_methods.R      # cross-validation with artificial gaps
├── 03_plot_results.R          # heatmaps, scatter plots, Taylor diagrams, network
├── 04_summary_tables.R        # summary spreadsheet (monthly/annual totals, statistics)
├── 05_fill_gaps_pcr_nipals.R  # fills the real gaps with PCR-NIPALS
├── R/                         # functions (sourced by the scripts)
├── data/                      # input data (not versioned, see below)
├── results/                   # outputs (created when the scripts run)
└── figures/                   # plots (created when the scripts run)
```

## How to run

1. Put the daily data file in `data/` (default name: `data/PLUVS_JOAO_LEITE.xlsx`; change `input_file` at the top of `01_prepare_data.R`).
2. Open R/RStudio with **this folder** as the working directory (`setwd(".../hydrology/precipitation-gap-filling")`).
3. Install the packages once:

```r
install.packages(c("readxl", "writexl", "ggplot2", "gridExtra", "igraph"))
```

4. Run the scripts in order:

```r
source("01_prepare_data.R")
source("02_validate_methods.R")
source("03_plot_results.R")
source("04_summary_tables.R")        # optional
source("05_fill_gaps_pcr_nipals.R")  # optional
```

## Input data format

One Excel file with **one sheet per rain gauge** (the sheet name is used as the gauge code). The first four columns must be, in this order: day, month, year and daily precipitation (e.g. `DIA`, `MÊS`, `ANO`, `PREC`). Precipitation values are divided by `scale_factor` (default `1e6`) to obtain mm; negative values are treated as missing.

A month containing any missing day is treated as a gap month.

## Settings

At the top of `02_validate_methods.R`:

| Parameter | Default | Effect |
|---|---|---|
| `percent_artificial_gaps` | `5` | share of monthly values removed for validation |
| `max_neighbours` | `10` | maximum number of neighbouring gauges per gap |
| `use_parallel` | `FALSE` | run the gap-by-gap methods on several CPU cores |
| `seed` | `1` | random seed for the artificial gaps |

## Notes

- Error codes: functions return negative codes (−995 to −999) when a gap cannot be filled. As in the original code, `fill_slr_variable`, `fill_linear_regression` and `fill_regional_vector` convert these codes to 0.
- The "complete" reference dataset is built with a nearest-neighbour filling (k = 1) that approximates MATLAB's `fillmissing(..., 'knn')`.
- Because the random number generators differ, the artificial gaps (and therefore the metrics) are not identical to the MATLAB version; all deterministic steps reproduce the MATLAB results.

## Author

Kamila Almeida dos Santos
