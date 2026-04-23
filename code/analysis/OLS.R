# ===============================================================
# OLS.R — Ordinary Least Squares GDP Validation
# Province and regency/city level, Sulawesi & Maluku
# Multiple model specifications with year fixed effects
#
# Run from its own directory: code/analysis/
# All outputs saved to: code/analysis/outputs/
#
# Licence (code): MIT — see repository LICENSE file.
# Input data uses Hansen GFC (CC BY 4.0) and Li et al. NTL (CC BY 4.0);
# outputs derived from those datasets must carry attribution.
# NASA MODIS data is CC0; GEBCO is Public Domain.
# Administrative unit definitions derived from BIG data (proprietary).
# ===============================================================

# ===============================================================
# LOAD LIBRARIES
# ===============================================================
library(modelsummary)
library(dplyr)
library(readr)

# ===============================================================
# PATHS
# ===============================================================
DATA_PROVINCE <- "../../data_processed/data_province_2002_2024.csv"
DATA_REGENCY  <- "../../data_processed/data_regency_city_2015_2024.csv"
OUTPUT_DIR    <- "outputs"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ===============================================================
# LOAD DATA
# ===============================================================
df_province <- read_csv(
  DATA_PROVINCE,
  show_col_types = FALSE,
  col_types = cols(province_id = col_character())
)

df_regency_city <- read_csv(
  DATA_REGENCY,
  show_col_types = FALSE,
  col_types = cols(regency_city_id = col_character())
)

cat("Province rows loaded:", nrow(df_province), "\n")
cat("Regency/city rows loaded:", nrow(df_regency_city), "\n")

# ===============================================================
# DATA PREPARATION — log transformations
# ===============================================================

# Province
df_province <- df_province %>%
  mutate(
    log_gdp                  = log(gdp),
    log_ndvi_sum             = log1p(ndvi_sum),
    log_forest_loss_area_ha  = log1p(forest_loss_area_ha),
    log_chlor_a_sum_B150_D15 = log1p(chlor_a_sum_B150_D15),
    log_ntl_sum_unfiltered   = log1p(ntl_sum_unfiltered)
  )

# Regency / city
df_regency_city <- df_regency_city %>%
  mutate(
    log_gdp                  = log(gdp),
    log_ndvi_sum             = log1p(ndvi_sum),
    log_forest_loss_area_ha  = log1p(forest_loss_area_ha),
    log_chlor_a_sum_B150_D15 = log1p(chlor_a_sum_B150_D15),
    log_ntl_sum_unfiltered   = log1p(ntl_sum_unfiltered)
  )

# Remove incomplete observations
cat("\nProvince rows before NA removal:", nrow(df_province), "\n")
df_province <- df_province[complete.cases(df_province), ]
cat("Province rows after NA removal:", nrow(df_province), "\n")

cat("Regency/city rows before NA removal:", nrow(df_regency_city), "\n")
df_regency_city <- df_regency_city[complete.cases(df_regency_city), ]
cat("Regency/city rows after NA removal:", nrow(df_regency_city), "\n")

# ===============================================================
# MODEL ESTIMATION — OLS with year fixed effects
#
# Model 1: GDP ~ NTL + year FE
# Model 2: GDP ~ NTL + NDVI + year FE
# Model 3: GDP ~ NTL + NDVI + Forest Loss + year FE
# Model 4: GDP ~ NTL + NDVI + Forest Loss + Chl-a + year FE
# ===============================================================

# Province
ols1_province <- lm(
  log_gdp ~ log_ntl_sum_unfiltered + factor(year),
  data = df_province
)
ols2_province <- lm(
  log_gdp ~ log_ntl_sum_unfiltered + log_ndvi_sum + factor(year),
  data = df_province
)
ols3_province <- lm(
  log_gdp ~ log_ntl_sum_unfiltered + log_ndvi_sum +
    log_forest_loss_area_ha + factor(year),
  data = df_province
)
ols4_province <- lm(
  log_gdp ~ log_ntl_sum_unfiltered + log_ndvi_sum +
    log_forest_loss_area_ha + log_chlor_a_sum_B150_D15 + factor(year),
  data = df_province
)

# Regency / city
ols1_regency_city <- lm(
  log_gdp ~ log_ntl_sum_unfiltered + factor(year),
  data = df_regency_city
)
ols2_regency_city <- lm(
  log_gdp ~ log_ntl_sum_unfiltered + log_ndvi_sum + factor(year),
  data = df_regency_city
)
ols3_regency_city <- lm(
  log_gdp ~ log_ntl_sum_unfiltered + log_ndvi_sum +
    log_forest_loss_area_ha + factor(year),
  data = df_regency_city
)
ols4_regency_city <- lm(
  log_gdp ~ log_ntl_sum_unfiltered + log_ndvi_sum +
    log_forest_loss_area_ha + log_chlor_a_sum_B150_D15 + factor(year),
  data = df_regency_city
)

# ===============================================================
# RESULTS EXPORT — regression tables with HC3 robust SEs
# ===============================================================

export_table <- function(models, file_path) {

  coef_labels <- c(
    log_ntl_sum_unfiltered   = "Log(NTL)",
    log_ndvi_sum             = "Log(NDVI)",
    log_forest_loss_area_ha  = "Log(Forest loss)",
    log_chlor_a_sum_B150_D15 = "Log(Chl-a)"
  )

  fe_row <- data.frame(
    part      = "gof",
    term      = "Year fixed effects",
    statistic = "",
    "NTL"                              = "Yes",
    "NTL + NDVI"                       = "Yes",
    "NTL + NDVI + Forest Loss"         = "Yes",
    "NTL + NDVI + Forest Loss + Chl-a" = "Yes",
    check.names = FALSE
  )

  tbl <- modelsummary(
    models,
    vcov      = "HC3",
    coef_map  = coef_labels,
    statistic = "({std.error})",
    stars     = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
    gof_omit  = "IC|AIC|BIC|Log.Lik|Std.Errors",
    add_rows  = fe_row,
    output    = "data.frame"
  )

  write.table(
    tbl,
    file      = file_path,
    sep       = ";",
    row.names = FALSE,
    quote     = FALSE
  )
}

models_province <- list(
  "NTL"                              = ols1_province,
  "NTL + NDVI"                       = ols2_province,
  "NTL + NDVI + Forest Loss"         = ols3_province,
  "NTL + NDVI + Forest Loss + Chl-a" = ols4_province
)

models_regency_city <- list(
  "NTL"                              = ols1_regency_city,
  "NTL + NDVI"                       = ols2_regency_city,
  "NTL + NDVI + Forest Loss"         = ols3_regency_city,
  "NTL + NDVI + Forest Loss + Chl-a" = ols4_regency_city
)

export_table(
  models_province,
  file.path(OUTPUT_DIR, "TABLE_Province_OLS_GDP.csv")
)
cat("Province OLS table exported.\n")

export_table(
  models_regency_city,
  file.path(OUTPUT_DIR, "TABLE_Regency_City_OLS_GDP.csv")
)
cat("Regency/city OLS table exported.\n")

cat("\n=== OLS.R complete — all outputs saved to", OUTPUT_DIR, "===\n")
