# =============================================================================
# STEP 5: FILL THE REAL GAPS WITH PCR-NIPALS
# Each gauge is filled using all the other gauges as neighbours.
# Requires: 01_prepare_data.R (results/prepared_data.rds)
# Output:   results/filled_data_pcr_nipals.xlsx (one sheet per gauge)
# Packages: writexl
# =============================================================================
rm(list = ls())
library(writexl)
for (f in c("R/utils.R", "R/nipals_pca.R", "R/fill_pcr_nipals.R")) source(f)

dp <- readRDS("results/prepared_data.rds")
all_monthly <- dp$monthly
gauges <- names(all_monthly)
output_file <- "results/filled_data_pcr_nipals.xlsx"

sheets <- list()
for (target in gauges) {
  neighbours <- all_monthly[setdiff(gauges, target)]   # all other gauges
  T <- all_monthly[[target]]
  T_filled <- T

  # Fill the gap months (MonthlyTotal < 0); unsuccessful fills keep the gap code
  for (r in seq_len(nrow(T))) {
    if (T$MonthlyTotal[r] < 0) {
      Px <- fill_pcr_nipals(target, all_monthly, neighbours, T$Year[r], T$Month[r], 3)
      if (Px >= 0) T_filled$MonthlyTotal[r] <- Px
    }
  }
  sheets[[target]] <- T_filled
}

write_xlsx(sheets, output_file)
cat(sprintf("PCR-NIPALS gap filling finished. File saved: %s\n", output_file))
