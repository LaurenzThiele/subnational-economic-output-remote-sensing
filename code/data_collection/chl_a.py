"""
chl_a.py — Chlorophyll-a concentration extraction from ocean colour imagery.

# =============================================================
# DATA SOURCE
# =============================================================
# Chlorophyll-a:
#   Product  : MODIS-Aqua Level-3 Chlorophyll-a (Version R2022.0)
#   Provider : NASA Ocean Biology Processing Group (OBPG)
#   Download : https://oceandata.sci.gsfc.nasa.gov/l3/
#              (accessed 28 October 2025)
#   Citation : NASA Ocean Biology Processing Group (2022). Moderate Resolution
#              Imaging Spectroradiometer (MODIS) Aqua Level-3 Chlorophyll-a
#              Data (2002–2024), Version R2022.0. NASA Ocean Biology DAAC.
#              https://oceandata.sci.gsfc.nasa.gov/l3/
#   Variable : chlor_a (chlorophyll-a concentration, mg m⁻³)
#   Format   : NetCDF (.nc) → converted to masked GeoTIFF in Step 1 below
#   Temporal : daily / 8-day composites (one file per period per year)
#   Source   : source_data/chl_a/{year}/*.tif  (after Step 1 conversion)
#
# Bathymetry (coastal zone definition):
#   Product  : GEBCO_2025 Grid — 15 arc-second global terrain model
#   Provider : GEBCO Bathymetric Compilation Group
#   Download : https://doi.org/10.5285/37c52e96-24ea-67ce-e063-7086abc05f29
#   Citation : GEBCO Bathymetric Compilation Group (2025). The GEBCO_2025
#              Grid — A Continuous Terrain Model for Oceans and Land at 15
#              Arc-Second Intervals. NERC EDS British Oceanographic Data
#              Centre NOC. https://doi.org/10.5285/37c52e96-24ea-67ce-e063-7086abc05f29
#   File     : source_data/chl_a/bathymetry/bathymetry_2025.tif
#   Unit     : depth in metres (negative = below sea level)
#
# Administrative boundaries loaded from the database (populated by administrative.py).
#
# Spatial coverage  : Sulawesi and Maluku regions (coastal + nearshore areas)
# Temporal coverage : 2002–2024
#
# =============================================================
# PROCESSING — STEP 1: NetCDF → masked GeoTIFF
# =============================================================
# For each .nc or .tif file in source_data/chl_a/{year}/:
#  1. Open the "chlor_a" variable; rename lon/lat dims to x/y.
#  2. Load bathymetry raster and build a "valid ocean" mask
#     (all depth values, i.e. bath >= -99999).
#  3. Reproject the mask to match the chl_a raster.
#  4. Apply the mask (set land to NaN).
#  5. If the result is entirely NaN, delete the file (corrupted).
#  6. Save as {stem}_B200m.tif with LZW compression; delete original.
#
# =============================================================
# PROCESSING — STEP 2: Spatial aggregation to admin boundaries
# =============================================================
# Three depth/distance coastal zone combinations (combos):
#   B100_D5  : depth ≥ −100 m, buffer ≤ 5 km  (shallow nearshore)
#   B150_D15 : depth ≥ −150 m, buffer ≤ 15 km (mid-shelf)
#   B200_D25 : depth ≥ −200 m, buffer ≤ 25 km (deep offshore)
#
# Geometry preprocessing (cached to GeoJSON, once per entity×combo):
#  1. Buffer the entity polygon in local UTM (metre-accurate).
#  2. Clip bathymetry to the buffer extent; build a "shallow zone" mask.
#  3. Convert the shallow mask to polygons; intersect with the buffer.
#  4. Union shallow polygons with the original entity polygon.
#  5. Fill interior holes (use exterior rings only).
#  6. Retain only sub-polygons that intersect the original entity.
#  7. Store the cleaned geometry in cleaned_geom_cache and write GeoJSON.
#
# Pixel extraction (year-first loop — one raster pass per file):
#  For each chl_a file in the year:
#    For each entity × combo:
#      1. Read the minimal raster window covering the cleaned geometry bounds.
#      2. Rasterize the cleaned geometry onto that window.
#      3. Apply valid-pixel mask (chl_a ≥ 0) and geometry mask.
#      4. Accumulate valid pixel arrays.
#
#  After all files for a year:
#    For each entity × combo:
#      1. Stack pixel arrays: (time × pixels); trim to minimum length.
#      2. Apply Savitzky-Golay smoothing along time axis
#         (window_length=7, polyorder=2, mode='nearest').
#      3. Flatten to 1D and compute finalize_stats (from utils).
#
# =============================================================
# OUTPUT
# =============================================================
# Cached GeoJSON geometries per entity×combo:
#   source_data/chl_a/bathymetry/{entity_id_safe}/{entity_id_safe}_{combo}.geojson
#
# output/table_chl_a_province_<DD_MM_YYYY>.sql
# output/table_chl_a_regency_city_<DD_MM_YYYY>.sql
#   INSERT INTO chl_a (level, entity_id, entity_name, year,
#     sum_B100_D5,   mean_B100_D5,   median_B100_D5,   max_B100_D5,   amplitude_B100_D5,
#     sum_B150_D15,  mean_B150_D15,  median_B150_D15,  max_B150_D15,  amplitude_B150_D15,
#     sum_B200_D25,  mean_B200_D25,  median_B200_D25,  max_B200_D25,  amplitude_B200_D25)
# =============================================================
"""

# =============================================================
# IMPORTS
# =============================================================
import json
import os
from pathlib import Path

import numpy as np
import rasterio
import geopandas as gpd
import shapely
import xarray as xr
import rioxarray as rio
from rasterio.features import shapes, rasterize
from rasterio.mask import mask as rio_mask
from scipy.signal import savgol_filter
from shapely import Polygon
from shapely.geometry import mapping, shape
from shapely.ops import unary_union
from tqdm import tqdm

from utils import (
    finalize_stats,
    write_sql_file,
    get_db_connection,
    load_province_boundaries_from_db,
    load_regency_city_boundaries_from_db,
    TARGET_REGIONS,
)

# =============================================================
# CONFIGURATION
# =============================================================
CHL_A_BASE          = "source_data/chl_a"
BATHYMETRY_FILE     = "source_data/chl_a/bathymetry/bathymetry_2025.tif"
GEOJSON_CACHE_BASE  = "source_data/chl_a/bathymetry"
OUTPUT_DIR          = "output"

YEARS = range(2002, 2025)

# Depth/distance coastal zone combos
COMBOS = [
    {"depth": -100, "distance": 5},    # B100_D5
    {"depth": -150, "distance": 15},   # B150_D15
    {"depth": -200, "distance": 25},   # B200_D25
]

SAVGOL_WINDOW    = 7
SAVGOL_POLYORDER = 2


# =============================================================
# STEP 1 — Convert NetCDF → masked GeoTIFF
# =============================================================

def convert_nc_to_tiff(source_base: str = CHL_A_BASE,
                        bathymetry_file: str = BATHYMETRY_FILE) -> None:
    """
    Convert all .nc and unmasked .tif files in source_base/{year}/ to
    bathymetry-masked GeoTIFFs named {stem}_B200m.tif.

    Files that are entirely NaN after masking or that cannot be read
    are deleted (they will be re-downloaded on the next run).
    """
    bath = rio.open_rasterio(bathymetry_file).squeeze()
    mask = xr.where(bath >= -99999, 1.0, np.nan)
    mask.rio.write_crs(bath.rio.crs, inplace=True)
    mask.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)

    all_files = []
    for year in YEARS:
        year_folder = Path(source_base) / str(year)
        if year_folder.exists():
            for f in year_folder.glob("*"):
                if f.suffix.lower() in [".tif", ".nc"] and "_B200m" not in f.stem:
                    all_files.append(f)

    for file in tqdm(all_files, desc="Step 1: NC → masked GeoTIFF"):
        try:
            if file.suffix.lower() == ".nc":
                with xr.open_dataset(file) as ds:
                    if "chlor_a" not in ds:
                        print(f"  Corrupted (deleted): {file} — 'chlor_a' missing")
                        file.unlink()
                        continue

                    chl_a = ds["chlor_a"].squeeze()
                    chl_a = chl_a.rename({"lon": "x", "lat": "y"})
                    chl_a.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
                    chl_a.rio.write_crs("EPSG:4326", inplace=True)

                    mask_reproj = mask.rio.reproject_match(chl_a)
                    chl_a_masked = chl_a.where(~np.isnan(mask_reproj))

                    if np.all(np.isnan(chl_a_masked.values)):
                        file.unlink()
                        continue

                    out_file = file.with_name(f"{file.stem}_B200m.tif")
                    chl_a_masked.rio.to_raster(out_file, compress="LZW")

                file.unlink()

            elif file.suffix.lower() == ".tif":
                with rio.open_rasterio(file) as tif_data:
                    if "band" in tif_data.dims:
                        tif_data = tif_data.squeeze("band", drop=True)

                    mask_reproj = mask.rio.reproject_match(tif_data)
                    tif_masked  = tif_data.where(~np.isnan(mask_reproj))

                    if np.all(np.isnan(tif_masked.values)):
                        file.unlink()
                        continue

                    out_file = file.with_name(f"{file.stem}_B200m.tif")
                    tif_masked.rio.to_raster(out_file, compress="LZW")

                file.unlink()

        except Exception:
            file.unlink()


# Run Step 1 before aggregation
convert_nc_to_tiff()


# =============================================================
# STEP 2 — Spatial aggregation helper functions
# =============================================================

def _build_coastal_geometry(reg_geom, entity_id_safe: str, combo: dict,
                              combo_label: str) -> object:
    """
    Build and cache the cleaned coastal geometry for one entity × combo.

    Buffers the entity polygon into UTM space, clips bathymetry to define
    the shallow coastal zone, unions with the original polygon, and fills
    interior holes to produce a contiguous coastal analysis area.

    Returns the cleaned Shapely geometry (or None on failure).
    """
    max_depth       = combo["depth"]
    max_distance_km = combo["distance"]

    try:
        utm_crs       = gpd.GeoSeries([reg_geom], crs="EPSG:4326").estimate_utm_crs()
        reg_geom_utm  = gpd.GeoSeries([reg_geom], crs="EPSG:4326").to_crs(utm_crs).iloc[0]
        buffer_utm    = reg_geom_utm.buffer(max_distance_km * 1000)
        buffer_geom   = gpd.GeoSeries([buffer_utm], crs=utm_crs).to_crs("EPSG:4326").iloc[0]

        with rasterio.open(BATHYMETRY_FILE) as src_bathy:
            out_image, out_transform = rio_mask(src_bathy, [mapping(buffer_geom)], crop=True)
            nodata = src_bathy.nodata

        shallow_mask = (out_image[0] < 0) & (out_image[0] >= max_depth)
        if nodata is not None:
            shallow_mask &= (out_image[0] != nodata)
        shallow_mask = shallow_mask.astype(np.uint8)

        shallow_polys = [
            shape(s)
            for s, v in shapes(shallow_mask, mask=shallow_mask, transform=out_transform)
            if v == 1
        ]

        if not shallow_polys:
            union_geom = unary_union([reg_geom])
        else:
            shallow_gdf = gpd.GeoDataFrame(
                geometry=shallow_polys, crs=src_bathy.crs
            ).to_crs("EPSG:4326")
            buffer_gdf     = gpd.GeoDataFrame(geometry=[buffer_geom], crs="EPSG:4326")
            shallow_within = gpd.overlay(shallow_gdf, buffer_gdf, how="intersection")
            union_geom     = unary_union([reg_geom] + list(shallow_within.geometry))

        if union_geom.geom_type == "Polygon":
            filled_geom = Polygon(union_geom.exterior)
        elif union_geom.geom_type == "MultiPolygon":
            filled_geom = unary_union([Polygon(p.exterior) for p in union_geom.geoms])
        else:
            filled_geom = reg_geom

        reg_multi    = shapely.MultiPolygon([reg_geom]) if reg_geom.geom_type == "Polygon" else reg_geom
        filled_multi = shapely.MultiPolygon([filled_geom]) if filled_geom.geom_type == "Polygon" else filled_geom

        valid_polys = [
            poly for poly in filled_multi.geoms
            if any(poly.intersects(ip) for ip in reg_multi.geoms)
        ]

        if len(valid_polys) == 0:
            cleaned = None
        elif len(valid_polys) == 1:
            cleaned = valid_polys[0]
        else:
            cleaned = shapely.MultiPolygon(valid_polys)

        return cleaned

    except Exception as e:
        print(f"  Geometry error for {entity_id_safe} {combo_label}: {e}")
        return None


def _precompute_geometries(entity_df, id_col: str) -> dict:
    """
    Precompute and cache cleaned coastal geometries for all entities × combos.

    Returns cleaned_geom_cache: {(entity_id_safe, combo_label): geometry or None}
    """
    cleaned_geom_cache = {}

    for _, entity in tqdm(entity_df.iterrows(), total=len(entity_df),
                          desc="  Precomputing coastal geometries"):
        entity_id      = entity[id_col]
        entity_id_safe = str(entity_id).replace(".", "_")
        reg_geom       = entity["geometry"]

        for combo in COMBOS:
            combo_label = f"B{abs(combo['depth'])}_D{combo['distance']}"

            cleaned = _build_coastal_geometry(reg_geom, entity_id_safe, combo, combo_label)
            cleaned_geom_cache[(entity_id_safe, combo_label)] = cleaned

            cache_folder = os.path.join(GEOJSON_CACHE_BASE, entity_id_safe)
            os.makedirs(cache_folder, exist_ok=True)
            geojson_path = os.path.join(cache_folder, f"{entity_id_safe}_{combo_label}.geojson")
            geojson_dict = gpd.GeoDataFrame(
                geometry=[cleaned], crs="EPSG:4326"
            ).__geo_interface__
            with open(geojson_path, "w", encoding="utf-8") as f:
                json.dump(geojson_dict, f)

    return cleaned_geom_cache


def _aggregate_chl_a(entity_df, id_col: str, name_col: str,
                      level: str, cleaned_geom_cache: dict) -> list:
    """
    Aggregate chl-a pixels to entity level for all years.

    Year-first loop: each raster file is opened once; all entities and
    combos are processed in the same pass to minimise I/O.

    Returns a list of SQL value strings (one per entity × year).
    """
    values_list = []

    for year in tqdm(YEARS, desc=f"  Processing years ({level})"):
        year_folder = os.path.join(CHL_A_BASE, str(year))
        year_files  = [
            os.path.join(year_folder, f)
            for f in os.listdir(year_folder)
            if f.lower().endswith(".tif")
        ]
        if not year_files:
            raise ValueError(f"No chl-a files for year {year}")

        pixel_store = {
            str(entity[id_col]).replace(".", "_"): {
                f"B{abs(c['depth'])}_D{c['distance']}": []
                for c in COMBOS
            }
            for _, entity in entity_df.iterrows()
        }

        for f in year_files:
            try:
                with rasterio.open(f) as src:
                    for _, entity in entity_df.iterrows():
                        entity_id_safe = str(entity[id_col]).replace(".", "_")
                        for combo in COMBOS:
                            combo_label  = f"B{abs(combo['depth'])}_D{combo['distance']}"
                            cleaned_geom = cleaned_geom_cache.get((entity_id_safe, combo_label))
                            if cleaned_geom is None:
                                continue

                            try:
                                geom_bounds = shapely.bounds(cleaned_geom)
                            except Exception:
                                continue
                            if not all(np.isfinite(geom_bounds)):
                                continue

                            try:
                                window = src.window(*geom_bounds).round_offsets().round_lengths()
                                if window.width <= 0 or window.height <= 0:
                                    continue
                            except Exception:
                                window = None

                            if window is not None:
                                arr       = src.read(1, window=window).astype(float)
                                transform = src.window_transform(window)
                            else:
                                arr       = src.read(1).astype(float)
                                transform = src.transform

                            try:
                                mask_arr = rasterize(
                                    [(mapping(cleaned_geom), 1)],
                                    out_shape=arr.shape,
                                    transform=transform,
                                    fill=0,
                                    dtype=np.uint8,
                                ).astype(bool)
                            except Exception as e:
                                print(f"  Rasterize error {entity_id_safe} {combo_label} {os.path.basename(f)}: {e}")
                                continue

                            valid_mask     = arr >= 0
                            final_mask     = valid_mask & mask_arr
                            masked_chl_a   = arr[final_mask]

                            if masked_chl_a.size > 0:
                                pixel_store[entity_id_safe][combo_label].append(masked_chl_a)

            except Exception as e:
                print(f"  Error opening {f}: {e}")
                continue

        for _, entity in entity_df.iterrows():
            entity_id      = entity[id_col]
            entity_id_safe = str(entity_id).replace(".", "_")
            entity_name    = entity[name_col]

            row = {"level": level, "entity_id": entity_id, "entity_name": entity_name, "year": year}

            for combo in COMBOS:
                combo_label = f"B{abs(combo['depth'])}_D{combo['distance']}"
                pixel_lists = pixel_store[entity_id_safe].get(combo_label, [])

                if not pixel_lists:
                    s, m, md, mx, amp = finalize_stats([])
                else:
                    try:
                        min_len              = min(a.size for a in pixel_lists)
                        pixel_stack_trimmed  = np.array([a[:min_len] for a in pixel_lists])
                        pixel_stack_smoothed = savgol_filter(
                            pixel_stack_trimmed,
                            window_length=SAVGOL_WINDOW,
                            polyorder=SAVGOL_POLYORDER,
                            axis=0,
                            mode="nearest",
                        )
                        all_pixels = pixel_stack_smoothed.flatten()
                    except Exception:
                        all_pixels = np.concatenate(pixel_lists)

                    s, m, md, mx, amp = finalize_stats([all_pixels])

                row[f"sum_{combo_label}"]       = s
                row[f"mean_{combo_label}"]      = m
                row[f"median_{combo_label}"]    = md
                row[f"max_{combo_label}"]       = mx
                row[f"amplitude_{combo_label}"] = amp

            values_list.append(row)

    columns    = list(values_list[0].keys())
    values_sql = []
    for row in values_list:
        vals = []
        for col in columns:
            v = row[col]
            if v is None or (isinstance(v, float) and np.isnan(v)):
                vals.append("0")
            elif isinstance(v, str):
                vals.append("'" + v.replace("'", "''") + "'")
            else:
                vals.append(str(v))
        values_sql.append(f"({', '.join(vals)})")

    return columns, values_sql


# =============================================================
# SECTION 1 — PROVINCE LEVEL
# =============================================================

print("Processing province-level chl-a …")

con          = get_db_connection()
province_gdf = load_province_boundaries_from_db(con, regions=TARGET_REGIONS)

province_cache = _precompute_geometries(province_gdf, "province_id")
columns, province_values = _aggregate_chl_a(
    province_gdf, "province_id", "province_name", "province", province_cache
)

sql_province = (
    f"INSERT INTO chl_a ({', '.join(columns)}) VALUES\n"
    + ",\n".join(province_values) + ";"
)
write_sql_file(sql_province, "chl_a_province", OUTPUT_DIR)


# =============================================================
# SECTION 2 — REGENCY / CITY LEVEL
# =============================================================

print("Processing regency/city-level chl-a …")

regency_gdf   = load_regency_city_boundaries_from_db(con, regions=TARGET_REGIONS)
con.close()

regency_cache = _precompute_geometries(regency_gdf, "regency_city_id")
columns, regency_values = _aggregate_chl_a(
    regency_gdf, "regency_city_id", "regency_city_name", "regency_city", regency_cache
)

sql_regency = (
    f"INSERT INTO chl_a ({', '.join(columns)}) VALUES\n"
    + ",\n".join(regency_values) + ";"
)
write_sql_file(sql_regency, "chl_a_regency_city", OUTPUT_DIR)
