# Replication code - Subnational GDP estimation via remote sensing

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19705003.svg)](https://doi.org/10.5281/zenodo.19705003)

This repository contains the replication code for a study examining the relationship between subnational gross regional domestic product (GRDP) and satellite-derived environmental indicators in Sulawesi and Maluku, Indonesia (2002–2024, province level; 2015–2024, regency/city level). The analysis combines Multiscale Geographically Weighted Regression (MGWR) with ordinary least squares (OLS) models including year fixed effects to assess spatial and temporal variation in these relationships.

The exact version of the code used in the study is archived on Zenodo and available at: https://doi.org/10.5281/zenodo.19705003

---

## Repository structure

```
.
├── code/
│   ├── data_collection/    # Python — builds and populates the database
│   └── analysis/           # R — regression models and outputs
├── data_processed/         # Panel datasets used as inputs to the R analysis
└── DATABASE/               # MySQL schema and container configuration
```

---

## Data collection pipeline

The `code/data_collection/` scripts extract and harmonise data from seven sources into a MySQL database. Each script is self-contained and documents its data source, processing steps, and output schema in its header.

| Script | Data source | Processing |
|--------|------------|------------|
| `administrative.py` | BIG boundary shapefiles | Derives type from name prefix; corrects known ID and geometry errors |
| `gdp.py` | BPS Statistics Indonesia (CSV + SIMDASI API) | Aligns two chain-linked GRDP series (base years 2000 and 2010) |
| `ndvi.py` | MODIS MOD13A1 v061 (16-day, 500 m) | Land-cover masking (IGBP classes 2/12/14); Savitzky-Golay smoothing |
| `forest_loss.py` | Hansen GFC v1.12 (30 m) | Zonal pixel counts by loss year; area conversion to hectares |
| `ntl.py` | Li et al. harmonised DMSP/VIIRS | Pixel-level 5-year moving median + MAD outlier correction |
| `chl_a.py` | MODIS-Aqua L3 chl-a + GEBCO 2025 bathymetry | Coastal zone masking at three depth/distance combos; Savitzky-Golay smoothing |

Shared utilities (CRS handling, SQL formatting, database I/O) are in `utils.py`. The database schema is in `DATABASE/schema.sql`.

---

## Analysis

`code/analysis/MGWR.R` runs `gwr.multiscale()` from the GWmodel package across three model specifications and ten years, exports global and local R² values, and produces validation figures.

`code/analysis/OLS.R` estimates pooled OLS models with year fixed effects at both province and regency/city level using HC3 robust standard errors.

Both scripts read from `data_processed/` and write all outputs to `code/analysis/outputs/`.

---

## Data sources

| Dataset | Provider | Licence |
|---------|----------|---------|
| Administrative boundaries | Badan Informasi Geospasial (BIG), 2025 | Proprietary |
| GRDP | Badan Pusat Statistik (BPS) | Open API |
| NDVI | NASA LP DAAC — MOD13A1 v061 (Didan, 2021) | CC0 |
| Land cover | NASA LP DAAC — MCD12Q1 v061 (Friedl & Sulla-Menashe, 2022) | CC0 |
| Chlorophyll-a | NASA OBPG — MODIS-Aqua R2022.0 | CC0 |
| Bathymetry | GEBCO Compilation Group, 2025 | Public domain |
| Forest loss | Hansen/UMD/Google/USGS/NASA — GFC v1.12 | CC BY 4.0 |
| Night-time light | Li et al., 2020 (figshare) | CC BY 4.0 |

Full citations are in [LICENSE](LICENSE).

---

## Citation

If you use this code, please cite:

LaurenzThiele. (2026). LaurenzThiele/research-subnational-economic-output-remote-sensing: Replication code - Subnational GDP estimation via remote sensing (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.19705004

---

## License

Code: MIT. See [LICENSE](LICENSE) for data attributions.
