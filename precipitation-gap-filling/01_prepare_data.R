# =============================================================================
# STEP 1: DATA PREPARATION AND CLEANING
# Reads the daily Excel file (one sheet per rain gauge), flags gaps,
# aggregates to monthly totals, computes long-term means and saves
# results/prepared_data.rds
#
# Run from the project folder:  source("01_prepare_data.R")
# Packages: readxl
# =============================================================================
rm(list = ls())
library(readxl)
source("R/utils.R")

# ---- Parameters ----------------------------------------------------------------
input_file  <- "data/PLUVS_JOAO_LEITE.xlsx"   # daily data, one sheet per gauge
scale_factor <- 1E6                           # raw values / scale_factor = mm
output_file <- "results/prepared_data.rds"
dir.create("results", showWarnings = FALSE)

# ---- 1. Load and clean daily data -----------------------------------------
cat("--- (STEP 1) Loading and cleaning daily data ---\n")
daily <- tryCatch({
  out <- list()
  for (sheet in excel_sheets(input_file)) {
    d <- as.data.frame(read_excel(input_file, sheet = sheet))
    # Expected columns: day, month, year, precipitation (e.g. DIA, MÊS, ANO, PREC)
    if (ncol(d) >= 4) names(d)[1:4] <- c("Day", "Month", "Year", "Precipitation")
    if ("Precipitation" %in% names(d)) {
      P <- d$Precipitation
      is_gap <- !is.na(P) & P < 0
      P[!is_gap] <- P[!is_gap] / scale_factor   # convert to mm
      P[is_gap] <- -999                         # gap code
      d$Precipitation <- P
    }
    out[[sheet]] <- d
  }
  out
}, error = function(e) stop(sprintf("Failed to load or clean the Excel file. Error: %s", conditionMessage(e))))
cat("Daily data loaded, cleaned and scaled.\n")

# ---- 2. Monthly aggregation -------------------------------------------------
# Equivalent of MATLAB groupsummary(..., 'sum'): sum ignoring NA + GroupCount.
# Note: gap days (-999) enter the sum, so any month containing a gap day
# becomes negative and is treated as a gap month (same logic as the original).
station_names <- names(daily)
monthly <- list()
for (st in station_names) {
  d <- daily[[st]]
  g <- split(d$Precipitation, list(d$Year, d$Month), drop = TRUE)
  keys <- do.call(rbind, strsplit(names(g), ".", fixed = TRUE))
  mt <- data.frame(Year = as.numeric(keys[, 1]),
                   Month = as.numeric(keys[, 2]),
                   GroupCount = sapply(g, length),
                   MonthlyTotal = sapply(g, sum, na.rm = TRUE))
  mt <- mt[order(mt$Year, mt$Month), ]
  rownames(mt) <- NULL
  monthly[[st]] <- mt
}
cat("Data aggregated to monthly totals.\n")

# ---- 3. Long-term means -------------------------------------------------------
annual_mean <- function(mt) {
  totals <- c()
  for (y in unique(mt$Year)) {
    dy <- mt[mt$Year == y, ]
    # only complete years (12 months) without gaps
    if (nrow(dy) == 12 && all(dy$MonthlyTotal >= 0)) totals <- c(totals, sum(dy$MonthlyTotal))
  }
  m <- if (length(totals) > 0) mean(totals) else NaN
  if (is.na(m)) m <- 0
  m
}

cat("Computing long-term means...\n")
annual_means <- list()
monthly_means <- list()
for (st in station_names) {
  annual_means[[st]] <- annual_mean(monthly[[st]])
  v <- monthly[[st]][monthly[[st]]$MonthlyTotal >= 0, ]
  mm <- numeric(12)
  if (nrow(v) > 0) {
    t <- tapply(v$MonthlyTotal, v$Month, mean)
    mm[as.numeric(names(t))] <- t
  }
  monthly_means[[st]] <- mm
}
cat("Annual and monthly means computed from valid data only.\n")

# ---- 4. Save ---------------------------------------------------------------------
saveRDS(list(monthly = monthly, annual_means = annual_means,
             monthly_means = monthly_means, station_names = station_names),
        output_file)
cat(sprintf("\n>>> PREPARED DATA SAVED TO: %s\n", output_file))
cat("--- (STEP 1) Done! ---\n")
