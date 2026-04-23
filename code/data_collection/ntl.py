"""
ntl.py — Night-Time Light (NTL) extraction with outlier correction.

# =============================================================
# DATA SOURCE
# =============================================================
# Product  : Harmonized Global Nighttime Light Dataset (DMSP + VIIRS, 1992–2024)
# Provider : Li, X.; Zhou, Y.; Zhao, M.; Zhao, X.
# Download : https://doi.org/10.6084/m9.figshare.9828827.v10
# Citation : Li et al. (2020). Harmonization of DMSP and VIIRS Nighttime Light
#            Data from 1992–2024 at the Global Scale (Version 10). figshare.
#            https://doi.org/10.6084/m9.figshare.9828827.v10
# Format   : GeoTIFF, one file per year (global coverage)
# Filename : NTL_Global_{year}.tif
# Source   : source_data/ntl/NTL_Global_{year}.tif
#
# Raster nodata: read from file metadata; falls back to -9999 if absent.
#
# Administrative boundaries loaded from the database (populated by administrative.py).
#
# Spatial coverage  : Sulawesi and Maluku regions
# Temporal coverage : analysis years 2002–2024
#                     (rasters loaded for 2000–2024 to support 5-year windows)
#
# =============================================================
# PROCESSING
# =============================================================
# Pixel-level 5-year moving median + MAD outlier detection:
#
#  1. Load rasters for 2000–2024; clip each to entity polygon via
#     zonal_stats(raster_out=True) and crop all-NaN borders.
#
#  2. For each entity × target year (2002–2024):
#     a. Assemble a 5-year window: offsets [−2, −1, 0, +1, +2].
#        If future years are unavailable, extend the window backwards.
#     b. Crop all arrays in the window to the same minimum shape.
#     c. Compute pixel-level moving median and MAD across the window.
#     d. Flag outliers: |DN_t − moving_median| > 3 × MAD.
#     e. Replace outlier pixels with the moving median (new_DN_t).
#
#  3. Summarise per entity and year:
#     Unfiltered (all DN > 0):
#       sum_unfiltered, data_pixels, data_pixels_greater_than_0,
#       imputed_outliers_unfiltered
#     Filtered (DN > 7, removes low-radiance background):
#       sum, data_pixels_greater_than_7, imputed_outliers
#
# =============================================================
# OUTPUT
# =============================================================
# output/table_ntl_province_<DD_MM_YYYY>.sql
# output/table_ntl_regency_city_<DD_MM_YYYY>.sql
#   INSERT INTO ntl (level, entity_id, entity_name, year,
#     data_pixels, data_pixels_greater_than_7, data_pixels_greater_than_0,
#     sum, imputed_outliers,
#     sum_unfiltered, imputed_outliers_unfiltered)
# =============================================================
"""

# =============================================================
# IMPORTS
# =============================================================
import os
import numpy as np
import pandas as pd
import rasterio
import geopandas as gpd
from rasterstats import zonal_stats
from tqdm import tqdm

from utils import (
    sql_num,
    write_sql_file,
    get_db_connection,
    load_province_boundaries_from_db,
    load_regency_city_boundaries_from_db,
    TARGET_REGIONS,
)

# =============================================================
# CONFIGURATION
# =============================================================
NTL_SOURCE_DIR          = "source_data/ntl"
OUTPUT_DIR              = "output"

LOAD_YEARS    = list(range(2000, 2025))   # years loaded for window construction
TARGET_YEARS  = range(2002, 2025)         # years written to SQL output
WINDOW_SIZE   = 5                         # target moving-window width (years)
DN_FILTER     = 7                         # DN threshold for "filtered" summary


# =============================================================
# CORE PROCESSING FUNCTION
# =============================================================

def _process_ntl(entity_gdf: gpd.GeoDataFrame,
                 id_col: str,
                 name_col: str,
                 level: str) -> list:
    """
    Run pixel-level NTL outlier detection for all entities in entity_gdf.

    Returns a list of SQL value strings (one per entity × target year).
    """

    # --- Load rasters and clip to each entity polygon --------
    # rasters dict: {entity_id: {year: 2D np.ndarray}}
    entity_ids   = entity_gdf[id_col].tolist()
    rasters      = {eid: {} for eid in entity_ids}
    current_gdf  = entity_gdf.copy()

    for year in tqdm(LOAD_YEARS, desc=f"  Loading rasters ({level})"):
        tif_path = os.path.join(NTL_SOURCE_DIR, f"NTL_Global_{year}.tif")
        if not os.path.exists(tif_path):
            print(f"  Missing raster: {tif_path}")
            continue

        with rasterio.open(tif_path) as src:
            nodata_value = src.nodata if src.nodata is not None else -9999
            raster_crs   = src.crs

        if current_gdf.crs != raster_crs:
            print(f"  ⓘ CRS mismatch — reprojecting for year {year}")
            current_gdf = current_gdf.to_crs(raster_crs)

        zs = zonal_stats(
            current_gdf,
            tif_path,
            raster_out=True,
            nodata=nodata_value,
            stats=None,
        )

        for i, row in enumerate(current_gdf.itertuples()):
            arr = zs[i]["mini_raster_array"].astype(float)
            arr = arr.filled(np.nan)   # masked pixels (outside polygon/nodata) → NaN

            # Crop all-NaN border rows and columns for efficiency
            valid_rows = ~np.all(np.isnan(arr), axis=1)
            valid_cols = ~np.all(np.isnan(arr), axis=0)
            arr_clean  = arr[np.ix_(valid_rows, valid_cols)]

            rasters[getattr(row, id_col)][year] = arr_clean

    # --- Pixel-level moving median + MAD ---------------------
    records = []

    for year in tqdm(LOAD_YEARS, desc=f"  Processing years ({level})"):
        for row in current_gdf.itertuples():
            eid = getattr(row, id_col)

            # Build 5-year window (symmetric around current year)
            arrs, offsets_used = [], []
            for offset in [-2, -1, 0, 1, 2]:
                arr = rasters[eid].get(year + offset)
                if arr is not None:
                    arrs.append(arr)
                    offsets_used.append(offset)

            # Expand backwards if window is still short
            while len(arrs) < WINDOW_SIZE:
                earliest_offset = min(offsets_used)
                prev_arr = rasters[eid].get(year + earliest_offset - 1)
                if prev_arr is None:
                    break
                arrs.insert(0, prev_arr)
                offsets_used.insert(0, earliest_offset - 1)

            if not arrs:
                continue

            # Crop to minimum shared shape across window arrays
            min_rows = min(a.shape[0] for a in arrs)
            min_cols = min(a.shape[1] for a in arrs)
            arrs_cropped = [a[:min_rows, :min_cols] for a in arrs]

            stack          = np.stack(arrs_cropped, axis=0)
            valid_pos      = ~np.all(np.isnan(stack), axis=0)
            moving_median  = np.full(stack.shape[1:], np.nan)
            mad            = np.full(stack.shape[1:], np.nan)

            moving_median[valid_pos] = np.nanmedian(stack[:, valid_pos], axis=0)
            mad[valid_pos]           = np.nanmedian(
                np.abs(stack[:, valid_pos] - moving_median[valid_pos]), axis=0
            )

            # Build per-pixel records for target years only
            if year not in TARGET_YEARS:
                continue

            offset_dict_template = dict(zip(offsets_used, range(len(arrs_cropped))))

            for r in range(min_rows):
                for c in range(min_cols):
                    DN_vals     = [a[r, c] for a in arrs_cropped]
                    offset_dict = dict(zip(offsets_used, DN_vals))
                    DN_t        = offset_dict.get(0, np.nan)
                    if np.isnan(DN_t):
                        continue

                    is_outlier = abs(DN_t - moving_median[r, c]) > 3 * mad[r, c]
                    new_DN_t   = moving_median[r, c] if is_outlier else DN_t

                    records.append({
                        "year":                 year,
                        id_col:                 eid,
                        name_col:               getattr(row, name_col),
                        "DN_t-2":               offset_dict.get(-2, np.nan),
                        "DN_t-1":               offset_dict.get(-1, np.nan),
                        "DN_t":                 DN_t,
                        "DN_t+1":               offset_dict.get(1, np.nan),
                        "DN_t+2":               offset_dict.get(2, np.nan),
                        "moving_median":        moving_median[r, c],
                        "mad":                  mad[r, c],
                        "is_outlier":           is_outlier,
                        "new_DN_t":             new_DN_t,
                    })

    df_pixels = pd.DataFrame(records)

    # --- Summarise per entity × year -------------------------

    # Unfiltered: all pixels with a valid DN
    summary_unfiltered = df_pixels.groupby([id_col, name_col, "year"]).agg(
        imputed_outliers_unfiltered=("new_DN_t", lambda x: (x != df_pixels.loc[x.index, "DN_t"]).sum()),
        sum_unfiltered=("new_DN_t", "sum"),
        data_pixels=("new_DN_t", "count"),
        data_pixels_greater_than_0=("new_DN_t", lambda x: (x > 0).sum()),
    ).reset_index()

    # Filtered: only pixels with DN > 7
    df_filtered = df_pixels[df_pixels["new_DN_t"] > DN_FILTER]
    summary_filtered = df_filtered.groupby([id_col, name_col, "year"]).agg(
        imputed_outliers=("new_DN_t", lambda x: (x != df_filtered.loc[x.index, "DN_t"]).sum()),
        sum=("new_DN_t", "sum"),
        data_pixels_greater_than_7=("new_DN_t", "count"),
    ).reset_index()

    summary = pd.merge(
        summary_filtered,
        summary_unfiltered,
        on=[id_col, name_col, "year"],
        how="outer",
    )
    summary["level"] = level

    summary = summary[[
        "level", id_col, name_col, "year",
        "data_pixels", "data_pixels_greater_than_7", "data_pixels_greater_than_0",
        "sum", "imputed_outliers",
        "sum_unfiltered", "imputed_outliers_unfiltered",
    ]]

    # Build SQL value strings
    values = [
        f"('{r['level']}', '{r[id_col]}', '{r[name_col]}', {r['year']}, "
        f"{sql_num(r['data_pixels'])}, {sql_num(r['data_pixels_greater_than_7'])}, "
        f"{sql_num(r['data_pixels_greater_than_0'])}, "
        f"{sql_num(r['sum'])}, {sql_num(r['imputed_outliers'])}, "
        f"{sql_num(r['sum_unfiltered'])}, {sql_num(r['imputed_outliers_unfiltered'])})"
        for _, r in summary.iterrows()
    ]
    return values


# =============================================================
# SECTION 1 — PROVINCE LEVEL
# =============================================================

print("Processing province-level NTL …")

con          = get_db_connection()
province_gdf = load_province_boundaries_from_db(con, regions=TARGET_REGIONS)

province_values = _process_ntl(province_gdf, "province_id", "province_name", "province")

sql_province = (
    "INSERT INTO ntl (level, entity_id, entity_name, year, "
    "data_pixels, data_pixels_greater_than_7, data_pixels_greater_than_0, "
    "sum, imputed_outliers, "
    "sum_unfiltered, imputed_outliers_unfiltered) VALUES\n"
    + ",\n".join(province_values) + ";"
)
write_sql_file(sql_province, "ntl_province", OUTPUT_DIR)


# =============================================================
# SECTION 2 — REGENCY / CITY LEVEL
# =============================================================

print("Processing regency/city-level NTL …")

regency_gdf = load_regency_city_boundaries_from_db(con, regions=TARGET_REGIONS)
con.close()

regency_values = _process_ntl(regency_gdf, "regency_city_id", "regency_city_name", "regency_city")

sql_regency = (
    "INSERT INTO ntl (level, entity_id, entity_name, year, "
    "data_pixels, data_pixels_greater_than_7, data_pixels_greater_than_0, "
    "sum, imputed_outliers, "
    "sum_unfiltered, imputed_outliers_unfiltered) VALUES\n"
    + ",\n".join(regency_values) + ";"
)
write_sql_file(sql_regency, "ntl_regency_city", OUTPUT_DIR)
