"""
ndvi.py — NDVI extraction and aggregation from MODIS imagery.

# =============================================================
# DATA SOURCE
# =============================================================
# Product  : MODIS/Terra Vegetation Indices 16-Day L3 Global 500m SIN Grid V061
#            (MOD13A1 v061)
# Provider : NASA EOSDIS Land Processes DAAC (LP DAAC)
# Download : https://doi.org/10.5067/MODIS/MOD13A1.061
#            via NASA APPEARS (appeears.earthdata.nasa.gov)
# Citation : Didan, K. (2021). MODIS/Terra Vegetation Indices 16-Day L3 Global
#            500m SIN Grid V061. NASA EOSDIS Land Processes DAAC.
#            https://doi.org/10.5067/MODIS/MOD13A1.061
# Band     : NDVI (scaled: raw × 0.0001 → true NDVI in [−1, 1])
# Format   : GeoTIFF tiles, one file per 16-day composite
# Filename : …doy{YYYY}{DOY}….tif  (year + day-of-year encoded in name)
# Source folders:
#   source_data/ndvi/{province_id}/         ← NDVI tiles per province
#   source_data/ndvi/land_cover/{province_id}/  ← MCD12Q1 land cover tiles
#
# Land cover product: MCD12Q1 v061 (MODIS/Terra+Aqua Land Cover Type 1, IGBP)
# Provider : NASA EOSDIS Land Processes DAAC (LP DAAC)
# Download : https://doi.org/10.5067/MODIS/MCD12Q1.061
# Citation : Friedl, M. and Sulla-Menashe, D. (2022). MODIS/Terra+Aqua Land Cover
#            Type Yearly L3 Global 500m SIN Grid V061. NASA EOSDIS Land Processes DAAC.
#            https://doi.org/10.5067/MODIS/MCD12Q1.061
# License  : CC0
# Valid land cover classes retained: 2 (Evergreen Broadleaf), 12 (Cropland),
#   14 (Cropland/Natural Vegetation Mosaic)
#
# Administrative boundaries loaded from the database (populated by administrative.py).
#
# Spatial coverage  : Sulawesi and Maluku regions
# Temporal coverage : 2002–2024 (annual aggregation from 16-day composites)
#
# =============================================================
# PROCESSING
# =============================================================
# For each entity (province or regency/city) × year:
#  1. Identify all NDVI tiles for that year (sorted chronologically by DOY).
#  2. Load the annual land cover tile.
#  3. Compute a geometry mask (pixels inside the entity polygon).
#  4. For each tile: apply NDVI validity mask (−2000–10000), land cover mask,
#     and geometry mask; store valid pixel arrays.
#  5. Build a consistent pixel mask across all timesteps (union of valid masks).
#  6. Stack valid pixels into a (time × pixels) array.
#  7. Apply Savitzky-Golay smoothing along the time axis
#     (window_length=7, polyorder=2, mode='nearest').
#  8. Compute annual per-pixel mean, then aggregate to entity statistics.
#
# MODIS scale factor: raw DN × 0.0001 → true NDVI
#
# =============================================================
# OUTPUT
# =============================================================
# output/table_ndvi_province_<DD_MM_YYYY>.sql
# output/table_ndvi_regency_city_<DD_MM_YYYY>.sql
#   INSERT INTO ndvi (level, entity_id, entity_name, year,
#                     sum, mean, median, max, amplitude)
# =============================================================
"""

# =============================================================
# IMPORTS
# =============================================================
import os
import numpy as np
import rasterio
from rasterio.features import geometry_mask
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
from scipy.signal import savgol_filter

from utils import (
    write_sql_file,
    get_db_connection,
    load_province_boundaries_from_db,
    load_regency_city_boundaries_from_db,
    TARGET_REGIONS,
)

# =============================================================
# CONFIGURATION
# =============================================================
NDVI_BASE               = "source_data/ndvi"
LAND_COVER_BASE         = "source_data/ndvi/land_cover"
OUTPUT_DIR              = "output"

YEARS               = range(2002, 2025)
VALID_NDVI_MIN      = -2000
VALID_NDVI_MAX      = 10000
NDVI_SCALE_FACTOR   = 0.0001
VALID_LC_CLASSES    = {2, 12, 14}        # IGBP: Evergreen Broadleaf, Cropland, Crop/Nat. Mosaic
SAVGOL_WINDOW       = 7
SAVGOL_POLYORDER    = 2


# =============================================================
# SHARED PROCESSING FUNCTION
# =============================================================

def _compute_ndvi_stats(entity_id: str, geom, year: int,
                         ndvi_folder: str, land_cover_folder: str) -> dict:
    """
    Compute annual NDVI summary statistics for one entity and year.

    Returns a dict with keys: sum, mean, median, max, amplitude.
    All values are 0.0 when no valid pixels are found.
    """
    all_ndvi_files = sorted(
        [f for f in os.listdir(ndvi_folder) if f.lower().endswith(".tif")],
        key=lambda x: (int(x.split("doy")[1][:4]), int(x.split("doy")[1][4:7])),
    )
    year_ndvi_files = [f for f in all_ndvi_files if int(f.split("doy")[1][:4]) == year]

    all_lc_files = [f for f in os.listdir(land_cover_folder) if f.lower().endswith(".tif")]
    year_lc_file = [f for f in all_lc_files if int(f.split("doy")[1][:4]) == year]
    if not year_lc_file:
        raise ValueError(f"No land cover file for year {year} in {land_cover_folder}")
    year_lc_file = year_lc_file[0]

    reg_geom = [geom]

    # Precompute geometry mask using the first NDVI tile's transform
    with rasterio.open(os.path.join(ndvi_folder, year_ndvi_files[0])) as src:
        geom_mask_arr = geometry_mask(
            reg_geom,
            transform=src.transform,
            invert=True,   # True = inside polygon
            out_shape=(src.height, src.width),
        )

    with rasterio.open(os.path.join(land_cover_folder, year_lc_file)) as lc_src:
        lc_arr = lc_src.read(1).astype(int)

    annual_ndvi_list = []
    masks_list       = []

    for ndvi_file in year_ndvi_files:
        with rasterio.open(os.path.join(ndvi_folder, ndvi_file)) as src:
            ndvi_arr = src.read(1).astype(float)

            # Validity mask: raw DN range
            ndvi_mask = (ndvi_arr >= VALID_NDVI_MIN) & (ndvi_arr <= VALID_NDVI_MAX)

            # Land cover mask: classes 2, 12, 14
            lc_mask = np.isin(lc_arr, list(VALID_LC_CLASSES))

            final_mask = ndvi_mask & lc_mask & geom_mask_arr
            masks_list.append(final_mask)

            # Apply MODIS scale factor
            annual_ndvi_list.append(ndvi_arr * NDVI_SCALE_FACTOR)

    # Union of valid masks across all timesteps
    consistent_mask = (
        np.logical_or.reduce(masks_list) if masks_list
        else np.zeros_like(lc_arr, dtype=bool)
    )

    if np.sum(consistent_mask) == 0:
        per_pixel_mean = np.array([])
    else:
        # Stack: (time × pixels)
        ndvi_stack = np.array([ndvi[consistent_mask] for ndvi in annual_ndvi_list])

        # Savitzky-Golay smoothing along time axis
        ndvi_smoothed = savgol_filter(
            ndvi_stack,
            window_length=SAVGOL_WINDOW,
            polyorder=SAVGOL_POLYORDER,
            axis=0,
            mode="nearest",
        )
        per_pixel_mean = np.mean(ndvi_smoothed, axis=0)

    if per_pixel_mean.size == 0:
        return dict(sum=0.0, mean=0.0, median=0.0, max=0.0, amplitude=0.0)

    annual_max = np.max(per_pixel_mean)
    annual_min = np.min(per_pixel_mean)
    return dict(
        sum       = round(float(np.sum(per_pixel_mean)),    2),
        mean      = round(float(np.mean(per_pixel_mean)),   4),
        median    = round(float(np.median(per_pixel_mean)), 4),
        max       = round(float(annual_max),                4),
        amplitude = round(float(annual_max - annual_min),   4),
    )


# =============================================================
# SECTION 1 — PROVINCE LEVEL
# =============================================================

print("Processing province-level NDVI …")

con           = get_db_connection()
province_gdf  = load_province_boundaries_from_db(con, regions=TARGET_REGIONS)
values_list   = []

for _, province in tqdm(province_gdf.iterrows(), total=len(province_gdf), desc="Provinces"):
    pid          = province["province_id"]
    ndvi_folder  = os.path.join(NDVI_BASE, pid)
    lc_folder    = os.path.join(LAND_COVER_BASE, pid)

    for year in YEARS:
        stats = _compute_ndvi_stats(pid, province["geometry"], year, ndvi_folder, lc_folder)
        values_list.append(
            f"('province', '{pid}', '{province['province_name']}', {year}, "
            f"{stats['sum']}, {stats['mean']}, {stats['median']}, "
            f"{stats['max']}, {stats['amplitude']})"
        )

sql_province = (
    "INSERT INTO ndvi (level, entity_id, entity_name, year, "
    "sum, mean, median, max, amplitude) VALUES\n"
    + ",\n".join(values_list) + ";"
)
write_sql_file(sql_province, "ndvi_province", OUTPUT_DIR)


# =============================================================
# SECTION 2 — REGENCY / CITY LEVEL
# =============================================================

print("Processing regency/city-level NDVI …")

regency_gdf = load_regency_city_boundaries_from_db(con, regions=TARGET_REGIONS)
con.close()
values_list = []

for _, regency in tqdm(regency_gdf.iterrows(), total=len(regency_gdf), desc="Regencies/cities"):
    rid         = regency["regency_city_id"]
    pid         = regency["province_id"]
    ndvi_folder = os.path.join(NDVI_BASE, pid)
    lc_folder   = os.path.join(LAND_COVER_BASE, pid)

    for year in YEARS:
        stats = _compute_ndvi_stats(rid, regency["geometry"], year, ndvi_folder, lc_folder)
        values_list.append(
            f"('regency_city', '{rid}', '{regency['regency_city_name']}', {year}, "
            f"{stats['sum']}, {stats['mean']}, {stats['median']}, "
            f"{stats['max']}, {stats['amplitude']})"
        )

sql_regency = (
    "INSERT INTO ndvi (level, entity_id, entity_name, year, "
    "sum, mean, median, max, amplitude) VALUES\n"
    + ",\n".join(values_list) + ";"
)
write_sql_file(sql_regency, "ndvi_regency_city", OUTPUT_DIR)
