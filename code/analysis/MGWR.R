# ===============================================================
# MGWR.R — Multiscale Geographically Weighted Regression
# GDP validation using remote sensing data, Sulawesi & Maluku,
# 2015–2024 (regency/city level)
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
library(sp)
library(GWmodel)
library(ggplot2)
library(dplyr)
library(readr)
library(purrr)
library(tibble)
library(forcats)
library(showtext)
library(grid)
library(scales)

# ===============================================================
# PATHS
# ===============================================================
DATA_REGENCY <- "../../data_processed/data_regency_city_2015_2024.csv"
OUTPUT_DIR   <- "outputs"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ===============================================================
# FONT SETUP
# Falls back to system defaults if Windows fonts are unavailable
# ===============================================================
tryCatch({
  font_add(
    family     = "Arial",
    regular    = "C:/Windows/Fonts/arial.ttf",
    bold       = "C:/Windows/Fonts/arialbd.ttf",
    italic     = "C:/Windows/Fonts/ariali.ttf",
    bolditalic = "C:/Windows/Fonts/arialbi.ttf"
  )
  font_add(
    family     = "Times New Roman",
    regular    = "C:/Windows/Fonts/times.ttf",
    bold       = "C:/Windows/Fonts/timesbd.ttf",
    italic     = "C:/Windows/Fonts/timesi.ttf",
    bolditalic = "C:/Windows/Fonts/timesbi.ttf"
  )
  showtext_auto()
  cat("Fonts loaded.\n")
}, warning = function(w) message("Font load warning: ", conditionMessage(w)),
   error   = function(e) message("Fonts not found — using system defaults."))

# ---------------------------------------------------------------
# FIGURE HELPERS
# ---------------------------------------------------------------

theme_apa7 <- function(base_family = "Times New Roman", base_size = 11) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line         = element_line(color = "black"),
      axis.ticks        = element_line(color = "black"),
      axis.title        = element_text(size = 12, face = "bold"),
      axis.text         = element_text(size = 11),
      axis.text.x       = element_text(margin = margin(t = 5)),
      panel.background  = element_blank(),
      plot.background   = element_blank(),
      plot.margin       = unit(c(6, 6, 6, 6), "mm"),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.key        = element_rect(fill = NA, color = NA),
      legend.key.size   = unit(0.4, "cm"),
      legend.text       = element_text(size = 10)
    )
}

draw_key_polygon3 <- function(data, params, size) {
  grid::rectGrob(
    width  = unit(0.25, "cm"),
    height = unit(0.25, "cm"),
    gp = gpar(
      col  = data$colour,
      fill = scales::alpha(data$fill, data$alpha),
      lwd  = 1
    )
  )
}

# ===============================================================
# LOAD DATA
# ===============================================================
df <- read_csv(
  DATA_REGENCY,
  show_col_types = FALSE,
  col_types = cols(regency_city_id = col_character())
)

cat("Rows loaded:", nrow(df), "\n")
cat("NA count:", sum(is.na(df)), "\n")

# ===============================================================
# DATA PREPARATION — log transformations
# ===============================================================
df <- df %>%
  mutate(
    log_gdp                  = log(gdp),
    log_ndvi_sum             = log1p(ndvi_sum),
    log_forest_loss_area_ha  = log1p(forest_loss_area_ha),
    log_chlor_a_sum_B150_D15 = log1p(chlor_a_sum_B150_D15),
    log_ntl_sum_unfiltered   = log1p(ntl_sum_unfiltered)
  )

# ===============================================================
# SPATIAL DATA — convert to SpatialPointsDataFrame
# Coordinates are UTM-projected (mixed zones per province);
# EPSG:32750 is used as a working CRS for distance calculations.
# ===============================================================
df_sp <- df
coordinates(df_sp) <- ~ lon + lat
proj4string(df_sp) <- CRS("+init=EPSG:32750")

# ===============================================================
# MODEL DEFINITIONS
# Model A: NTL only                         (baseline)
# Model B: NTL + NDVI + Forest Loss         (rural extension)
# Model C: NTL + NDVI + Forest Loss + Chl-a (coastal extension)
# ===============================================================
model_specs <- list(
  A = list(
    label   = "Model A (NTL)",
    formula = log_gdp ~ log_ntl_sum_unfiltered,
    n_terms = 2
  ),
  B = list(
    label   = "Model B (NTL + NDVI + Forest Loss)",
    formula = log_gdp ~ log_ntl_sum_unfiltered + log_ndvi_sum + log_forest_loss_area_ha,
    n_terms = 4
  ),
  C = list(
    label   = "Model C (NTL + NDVI + Forest Loss + Chl-a)",
    formula = log_gdp ~ log_ntl_sum_unfiltered + log_ndvi_sum +
                        log_forest_loss_area_ha + log_chlor_a_sum_B150_D15,
    n_terms = 5
  )
)

years <- 2015:2024

# ===============================================================
# MODEL ESTIMATION — MGWR loop (3 models × 10 years)
# ===============================================================
mgwr_models <- list(A = list(), B = list(), C = list())

for (spec_id in names(model_specs)) {

  spec <- model_specs[[spec_id]]
  cat("\n\n########## Running MGWR:", spec$label, "##########\n")

  for (yr in years) {

    cat("\n========== Year:", yr, "==========\n")
    df_year <- subset(df_sp, year == yr)

    mgwr_models[[spec_id]][[as.character(yr)]] <- tryCatch({
      gwr.multiscale(
        formula   = spec$formula,
        data      = df_year,
        criterion = "dCVR",
        kernel    = "bisquare",
        adaptive  = TRUE,
        bws0      = rep(nrow(df_year), spec$n_terms),
        verbose   = FALSE
      )
    }, error = function(e) {
      cat("ERROR", spec$label, "year", yr, ":", conditionMessage(e), "\n")
      NULL
    })

    model_yr <- mgwr_models[[spec_id]][[as.character(yr)]]
    if (!is.null(model_yr)) {
      bws_matrix <- model_yr[[5]]
      bws_final  <- bws_matrix[nrow(bws_matrix), ]
      term_names <- all.vars(spec$formula)[-1]
      cat(spec$label, "bandwidths for", yr, ":\n")
      cat("  Intercept:", bws_final[1], "\n")
      for (j in seq_along(term_names)) {
        cat(" ", term_names[j], ":", bws_final[j + 1], "\n")
      }
    }
  }
}

# ===============================================================
# RESULTS EXPORT — Global R² (all 3 models)
# ===============================================================
r2_all <- map_dfr(names(model_specs), function(spec_id) {

  spec <- model_specs[[spec_id]]

  map_dfr(years, function(yr) {

    model   <- mgwr_models[[spec_id]][[as.character(yr)]]
    df_year <- subset(df_sp, year == yr)
    y_obs   <- df_year$log_gdp

    r2 <- if (!is.null(model)) {
      y_hat  <- model$SDF$yhat
      ss_res <- sum((y_obs - y_hat)^2)
      ss_tot <- sum((y_obs - mean(y_obs))^2)
      1 - ss_res / ss_tot
    } else NA

    tibble(Model = spec$label, Year = as.character(yr), Value = r2)
  })
})

r2_means <- r2_all %>%
  group_by(Model) %>%
  summarise(Year = "Mean", Value = mean(Value, na.rm = TRUE), .groups = "drop")

r2_all <- bind_rows(r2_all, r2_means)

write_csv(r2_all, file.path(OUTPUT_DIR, "DATA_model_all_mgwr_R2_per_year.csv"))
cat("\nGlobal R² exported.\n")
print(r2_all)

# ===============================================================
# RESULTS EXPORT — Pointwise local R² (all models; Model C only written)
#
# Local R² is computed at each predictor's own converged bandwidth
# and then averaged across predictors, respecting MGWR's multi-scale
# structure. Each predictor's local window reflects its spatial scale.
# ===============================================================
for (spec_id in names(model_specs)) {

  spec             <- model_specs[[spec_id]]
  combined_results <- NULL

  cat("\n\n########## Local R²:", spec$label, "##########\n")

  for (yr in years) {

    model <- mgwr_models[[spec_id]][[as.character(yr)]]
    if (is.null(model)) next

    df_year <- subset(df_sp, year == yr)
    coords  <- coordinates(df_year)
    y_obs   <- df_year$log_gdp
    y_hat   <- model$SDF$yhat
    n       <- nrow(coords)

    bws_matrix <- model[[5]]
    bws_final  <- bws_matrix[nrow(bws_matrix), ]
    pred_bws   <- bws_final[-1]

    cat("  Year:", yr, "| n:", n, "| predictor bandwidths:",
        paste(round(pred_bws), collapse = ", "), "\n")

    local_r2_runs <- lapply(pred_bws, function(k) {
      k          <- min(max(4, round(k)), n)
      local_r2_k <- numeric(n)
      for (i in seq_len(n)) {
        dists        <- sqrt((coords[, 1] - coords[i, 1])^2 +
                               (coords[, 2] - coords[i, 2])^2)
        nn_idx        <- order(dists)[1:k]
        y_local       <- y_obs[nn_idx]
        yhat_local    <- y_hat[nn_idx]
        ss_res_local  <- sum((y_local - yhat_local)^2)
        ss_tot_local  <- sum((y_local - mean(y_local))^2)
        local_r2_k[i] <- ifelse(ss_tot_local == 0, NA,
                                1 - ss_res_local / ss_tot_local)
      }
      local_r2_k
    })

    local_r2 <- Reduce("+", local_r2_runs) / length(local_r2_runs)

    year_results <- data.frame(
      regency_city_id = df_year$regency_city_id,
      Local_R2        = local_r2
    ) %>%
      rename(!!paste0("local_R2_", yr) := Local_R2) %>%
      left_join(
        df %>%
          filter(year == yr) %>%
          select(regency_city_id, regency_city_name, type, province_id, geo_region, lon, lat),
        by = "regency_city_id"
      )

    if (is.null(combined_results)) {
      combined_results <- year_results
    } else {
      combined_results <- combined_results %>%
        left_join(
          year_results %>% select(regency_city_id, !!paste0("local_R2_", yr)),
          by = "regency_city_id"
        )
    }
  }

  r2_cols          <- paste0("local_R2_", years)
  combined_results <- combined_results %>%
    rowwise() %>%
    mutate(mean_local_R2 = mean(c_across(all_of(r2_cols)), na.rm = TRUE)) %>%
    ungroup()

  if (spec_id == "C") {
    write_csv(combined_results,
              file.path(OUTPUT_DIR, "DATA_model_c_mgwr_local_R2_per_year.csv"))
    cat("Model C local R² exported.\n")
  }
}

# ===============================================================
# RESULTS EXPORT — Local model coefficients (Model C only)
# ===============================================================
for (spec_id in names(model_specs)) {

  mgwr_all <- bind_rows(
    lapply(years, function(yr) {

      model <- mgwr_models[[spec_id]][[as.character(yr)]]
      if (is.null(model)) return(NULL)

      res        <- as.data.frame(model$SDF@data)
      df_year_sp <- subset(df_sp, year == yr)

      res$regency_city_id <- df_year_sp$regency_city_id
      res$year            <- yr

      return(res)
    })
  )

  if (spec_id == "C") {
    write_csv(mgwr_all,
              file.path(OUTPUT_DIR, "DATA_model_c_mgwr_local_model_per_year.csv"))
    cat("Model C local coefficients exported.\n")
  }
}

# ===============================================================
# FIGURE — Global R² per year (3 models, bar chart)
# ===============================================================
r2_df <- read_csv(
  file.path(OUTPUT_DIR, "DATA_model_all_mgwr_R2_per_year.csv"),
  show_col_types = FALSE
)

r2_df$Model <- fct_relevel(
  r2_df$Model,
  "Model A (NTL)",
  "Model B (NTL + NDVI + Forest Loss)",
  "Model C (NTL + NDVI + Forest Loss + Chl-a)"
)

model_levels <- c(
  "Model A (NTL)",
  "Model B (NTL + NDVI + Forest Loss)",
  "Model C (NTL + NDVI + Forest Loss + Chl-a)"
)

r2_plot <- r2_df %>%
  mutate(Model = fct_relevel(Model, model_levels)) %>%
  mutate(Year = ifelse(Year == "Mean", "Mean", Year)) %>%
  mutate(Year = factor(
    Year,
    levels = c(sort(unique(Year[Year != "Mean"])), "Mean")
  ))

fig_r2_per_year <- ggplot(r2_plot, aes(x = Year, y = Value, fill = Model)) +
  geom_col(
    width     = 0.75,
    position  = position_dodge(width = 0.85),
    color     = "black",
    key_glyph = draw_key_polygon3
  ) +
  scale_fill_manual(
    values = c(
      "Model A (NTL)"                              = "grey85",
      "Model B (NTL + NDVI + Forest Loss)"         = "grey60",
      "Model C (NTL + NDVI + Forest Loss + Chl-a)" = "grey35"
    ),
    guide = guide_legend(reverse = FALSE)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.02)),
    limits = c(0, 0.8),
    breaks = seq(0, 0.8, 0.1),
    labels = label_number(decimal.mark = ".", accuracy = 0.01, trim = FALSE)
  ) +
  labs(x = "Year", y = expression(bold(R^2)), fill = NULL) +
  theme_apa7(base_family = "Times New Roman") +
  theme(
    axis.ticks.x         = element_blank(),
    axis.title.x         = element_blank(),
    axis.text.x          = element_text(angle = 0, vjust = 0.5),
    legend.position      = c(0.98, 0.98),
    legend.justification = c("right", "top")
  )

ggsave(
  file.path(OUTPUT_DIR, "FIGURE_mgwr_R2_per_year.pdf"),
  fig_r2_per_year,
  width = 6, height = 4, device = cairo_pdf
)
cat("Figure: MGWR global R2 per year exported.\n")

# ===============================================================
# TABLE — Local R² descriptive statistics (Model C)
# ===============================================================
r2_local_df <- read_csv(
  file.path(OUTPUT_DIR, "DATA_model_c_mgwr_local_R2_per_year.csv"),
  show_col_types = FALSE,
  col_types = cols(regency_city_id = col_character(), .default = col_guess())
)

compute_r2_stats <- function(x) {
  x <- x[!is.na(x)]
  data.frame(
    Min    = round(min(x),            4),
    Q1     = round(quantile(x, 0.25), 4),
    Median = round(median(x),         4),
    Mean   = round(mean(x),           4),
    Q3     = round(quantile(x, 0.75), 4),
    Max    = round(max(x),            4),
    SD     = round(sd(x),             4),
    Range  = round(max(x) - min(x),   4)
  )
}

r2_year_cols  <- paste0("local_R2_", 2015:2024)
r2_stats_year <- bind_rows(
  lapply(r2_year_cols, function(col) {
    cbind(Year = sub("local_R2_", "", col), compute_r2_stats(r2_local_df[[col]]))
  })
)

r2_stats_overall <- cbind(Year = "Total", compute_r2_stats(r2_local_df$mean_local_R2))
r2_stats_table   <- bind_rows(r2_stats_year, r2_stats_overall)

cat("\n=========== LOCAL R² DESCRIPTIVE STATISTICS (Model C) ===========\n")
print(r2_stats_table, row.names = FALSE)

write_csv(r2_stats_table,
          file.path(OUTPUT_DIR, "TABLE_model_c_local_R2_descriptive.csv"))
cat("Table: local R2 descriptive statistics exported.\n")

# ===============================================================
# NOTE: Coefficient maps (spatial visualisation) require polygon
# geometry data (regency/city shapefiles) which are not included
# in the CSV datasets. To produce these maps, load a shapefile,
# join on regency_city_id, and render with ggplot2 + geom_sf().
# The coefficient data is available in:
#   DATA_model_c_mgwr_local_model_per_year.csv
# ===============================================================
message("Coefficient maps skipped: spatial polygon geometry not available in CSV inputs.")

# ===============================================================
# TABLE + FIGURE — Validation metrics & scatter plot (Model C)
# ===============================================================

# Load Model C local coefficients
mgwr_coef_val <- read_csv(
  file.path(OUTPUT_DIR, "DATA_model_c_mgwr_local_model_per_year.csv"),
  col_types = cols(
    regency_city_id          = col_character(),
    year                     = col_integer(),
    Intercept                = col_double(),
    log_ntl_sum_unfiltered   = col_double(),
    log_ndvi_sum             = col_double(),
    log_forest_loss_area_ha  = col_double(),
    log_chlor_a_sum_B150_D15 = col_double()
  )
) %>%
  rename(
    intercept                            = Intercept,
    coef_ntl_sum_unfiltered              = log_ntl_sum_unfiltered,
    coef_ndvi_sum                        = log_ndvi_sum,
    coef_forest_loss_area_ha             = log_forest_loss_area_ha,
    coef_chlor_a_sum_B150_D15            = log_chlor_a_sum_B150_D15
  )

df_val <- df %>%
  mutate(
    regency_city_id = as.character(regency_city_id),
    year            = as.integer(year)
  )

# Diagnostic: flag any observations without MGWR coefficients
missing_coef <- df_val %>%
  left_join(
    mgwr_coef_val %>% select(regency_city_id, year, intercept),
    by = c("regency_city_id", "year")
  ) %>%
  filter(is.na(intercept)) %>%
  distinct(regency_city_id, year, regency_city_name)

if (nrow(missing_coef) > 0) {
  cat("WARNING: Missing MGWR coefficients for", nrow(missing_coef),
      "regency-year observations:\n")
  print(missing_coef)
} else {
  cat("No missing MGWR coefficients.\n")
}

# Merge and compute predicted values
df_final <- df_val %>%
  left_join(mgwr_coef_val, by = c("regency_city_id", "year")) %>%
  mutate(
    log_gdp_estimate =
      intercept +
      coef_ntl_sum_unfiltered   * log_ntl_sum_unfiltered   +
      coef_ndvi_sum             * log_ndvi_sum             +
      coef_forest_loss_area_ha  * log_forest_loss_area_ha  +
      coef_chlor_a_sum_B150_D15 * log_chlor_a_sum_B150_D15,
    error = log_gdp_estimate - log_gdp
  )

# Compute validation metrics
compute_metrics <- function(data, label) {
  data.frame(
    Year    = label,
    N       = sum(!is.na(data$error)),
    RMSE    = round(sqrt(mean(data$error^2,  na.rm = TRUE)), 4),
    MAE     = round(mean(abs(data$error),    na.rm = TRUE),  4),
    Bias    = round(mean(data$error,         na.rm = TRUE),  4),
    Pearson = round(cor(data$log_gdp, data$log_gdp_estimate,
                        use = "complete.obs"), 4),
    stringsAsFactors = FALSE
  )
}

metrics_year <- df_final %>%
  group_by(year) %>%
  group_map(~ compute_metrics(.x, as.character(.y$year))) %>%
  bind_rows()

metrics_overall <- compute_metrics(df_final, "Total")
metrics_table   <- bind_rows(metrics_year, metrics_overall)

cat("\n=========== VALIDATION METRICS (Model C) ===========\n")
print(metrics_table, row.names = FALSE)

write_csv(metrics_table,
          file.path(OUTPUT_DIR, "TABLE_validation_metrics.csv"))
cat("Validation metrics exported.\n")

# Scatter plot: predicted vs actual log(GDP)
axis_min <- floor(min(df_final$log_gdp, df_final$log_gdp_estimate, na.rm = TRUE))
axis_max <- ceiling(max(df_final$log_gdp, df_final$log_gdp_estimate, na.rm = TRUE))

fig_scatter <- ggplot(df_final, aes(x = log_gdp, y = log_gdp_estimate)) +
  geom_point(alpha = 0.45, size = 1.5, color = "grey45", shape = 16) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_x_continuous(
    limits = c(axis_min, axis_max),
    breaks = seq(axis_min, axis_max, by = 1),
    labels = label_number(decimal.mark = ".", accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    limits = c(axis_min, axis_max),
    breaks = seq(axis_min, axis_max, by = 1),
    labels = label_number(decimal.mark = ".", accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = "Actual log(GDP)", y = "Predicted log(GDP)") +
  coord_fixed(ratio = 1) +
  theme_apa7(base_family = "Times New Roman")

ggsave(
  file.path(OUTPUT_DIR, "FIGURE_validation_predicted_vs_actual.pdf"),
  fig_scatter,
  width = 5, height = 5, device = cairo_pdf
)
cat("Scatter plot exported.\n")

cat("\n=== MGWR.R complete — all outputs saved to", OUTPUT_DIR, "===\n")
