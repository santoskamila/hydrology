# =============================================================================
# STEP 4: SUMMARY SPREADSHEET (no plots)
# Sheets:
#   - Monthly_totals:   Gauge | Year | Month | Total_mm
#   - Annual_totals:    Gauge | Year | Total_mm
#   - Monthly_climatology
#   - Annual_mean
#   - General_statistics
# Months flagged as gaps (negative, NA or > 3000 mm) are rebuilt from the
# daily file when it is available; otherwise they are left as NA.
#
# Requires: 01_prepare_data.R (results/prepared_data.rds)
# Packages: readxl, writexl
# =============================================================================
library(readxl); library(writexl)
source("R/normalize_daily_units.R")

summary_tables <- function(rds_in = "results/prepared_data.rds",
                           xlsx_out = "results/precipitation_summary.xlsx",
                           daily_xlsx = "data/PLUVS.xlsx") {   # fallback only

  stopifnot("results/prepared_data.rds not found." = file.exists(rds_in))
  S <- readRDS(rds_in)

  monthly_all <- NULL; annual_all <- NULL; clim_all <- NULL; ann_mean_all <- NULL; stats_all <- NULL

  for (st in S$station_names) {
    obj <- S$monthly[[st]]
    Year <- obj$Year; Month <- obj$Month; Total <- obj$MonthlyTotal

    # Sanity check of monthly values
    bad <- is.na(Total) | Total < 0 | Total > 3000
    if (any(bad)) {
      cat(sprintf("Warning: invalid monthly values in %s — trying to rebuild from daily data...\n", st))
      if (file.exists(daily_xlsx) && st %in% excel_sheets(daily_xlsx)) {
        dd <- coerce_daily(as.data.frame(read_excel(daily_xlsx, sheet = st)))
        P <- normalize_daily_units(dd$P)
        g <- aggregate(list(Total = P), by = list(Year = dd$Year, Month = dd$Month),
                       FUN = function(v) sum(v, na.rm = TRUE))
        g <- g[order(g$Year, g$Month), ]
        Year <- g$Year; Month <- g$Month; Total <- g$Total
        cat(sprintf("  -> Rebuilt from %s: %s.\n", daily_xlsx, st))
      } else {
        cat(sprintf("  -> %s/%s not available; invalid months kept as NA.\n", daily_xlsx, st))
        Total[bad] <- NA
      }
    }

    Tm <- data.frame(Gauge = st, Year = Year, Month = Month, Total_mm = Total)
    monthly_all <- rbind(monthly_all, Tm)

    # Annual totals (sum of months)
    Ta <- aggregate(list(Total_mm = Tm$Total_mm), by = list(Year = Tm$Year),
                    FUN = function(v) sum(v, na.rm = TRUE))
    Ta <- data.frame(Gauge = st, Ta)
    annual_all <- rbind(annual_all, Ta)

    stats_all <- rbind(stats_all, data.frame(
      Gauge = st,
      Monthly_Mean = mean(Tm$Total_mm, na.rm = TRUE),
      Monthly_Max = max(Tm$Total_mm, na.rm = TRUE),
      Monthly_Min = min(Tm$Total_mm, na.rm = TRUE),
      Annual_Mean = mean(Ta$Total_mm, na.rm = TRUE),
      Annual_Max = max(Ta$Total_mm, na.rm = TRUE),
      Annual_Min = min(Ta$Total_mm, na.rm = TRUE)))

    if (length(S$monthly_means[[st]]) >= 12) {
      clim_all <- rbind(clim_all, data.frame(Gauge = st, Month = 1:12,
                                             Climatological_Mean_mm = S$monthly_means[[st]][1:12]))
    }
    if (!is.null(S$annual_means[[st]])) {
      ann_mean_all <- rbind(ann_mean_all, data.frame(Gauge = st, Annual_Mean_mm = S$annual_means[[st]]))
    }
  }

  sort_by <- function(d, cols) if (is.null(d)) d else d[do.call(order, d[cols]), ]
  sheets <- list(Monthly_totals = sort_by(monthly_all, c("Gauge", "Year", "Month")),
                 Annual_totals = sort_by(annual_all, c("Gauge", "Year")))
  if (!is.null(clim_all)) sheets$Monthly_climatology <- sort_by(clim_all, c("Gauge", "Month"))
  if (!is.null(ann_mean_all)) sheets$Annual_mean <- sort_by(ann_mean_all, "Gauge")
  sheets$General_statistics <- stats_all

  # If the file is open/locked, save under another name
  out <- xlsx_out
  ok <- tryCatch({ write_xlsx(sheets, out); TRUE }, error = function(e) FALSE)
  if (!ok) {
    out <- sub("\\.xlsx$", format(Sys.time(), "_%Y%m%d_%H%M%S.xlsx"), xlsx_out)
    cat(sprintf("File in use. Saving as: %s\n", out))
    write_xlsx(sheets, out)
  }
  cat(sprintf("OK: file written -> %s\n", out))
  invisible(sheets)
}

# ---- Helpers ---------------------------------------------------------------------
norm_header <- function(s) {
  # lower case, trimmed, without accents
  s <- tolower(trimws(s))
  chartr("ãâáàäêéèëíìïîõôóòöúùüûç", "aaaaaeeeeiiiiooooouuuuc", s)
}

coerce_daily <- function(T) {
  # Detects day/month/year columns by name (Portuguese or English, with or
  # without accents) or by value range, and the precipitation column.
  nm <- norm_header(names(T))
  num <- lapply(T, function(x) suppressWarnings(as.numeric(x)))

  i_day <- which(nm %in% c("dia", "day"))[1]
  i_month <- which(nm %in% c("mes", "month"))[1]
  i_year <- which(nm %in% c("ano", "year"))[1]
  in_range <- function(lo, hi) which(sapply(num, function(v) { v <- v[!is.na(v)]; length(v) > 0 && all(v >= lo & v <= hi) }))
  if (is.na(i_day)) i_day <- in_range(1, 31)[1]
  if (is.na(i_month)) i_month <- setdiff(in_range(1, 12), i_day)[1]
  if (is.na(i_year)) i_year <- in_range(1900, 2100)[1]
  if (any(is.na(c(i_day, i_month, i_year)))) stop("Day/month/year columns not found in the daily file.")

  i_p <- which(nm %in% c("prec", "precipitacao", "precipitation", "chuva", "rain", "p", "mm", "valor"))[1]
  if (is.na(i_p)) {
    for (k in setdiff(seq_along(num), c(i_day, i_month, i_year))) {
      v <- num[[k]]
      if (all(is.na(v)) || length(unique(v[!is.na(v)])) <= 2) next
      i_p <- k; break
    }
    if (is.na(i_p)) stop("Precipitation column not found in the daily file.")
  }
  list(Year = num[[i_year]], Month = num[[i_month]], Day = num[[i_day]], P = num[[i_p]])
}

summary_tables()
