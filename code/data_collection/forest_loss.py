"""
forest_loss.py — Forest loss extraction from Global Forest Watch raster.

# =============================================================
# DATA SOURCE
# =============================================================
# Product  : Hansen/UMD/Google/USGS/NASA Global Forest Change v1.12
#            (Global Forest Watch — forest loss year layer)
# Provider : University of Maryland / Global Forest Watch
# Download : https://storage.googleapis.com/earthenginepartners-hansen/
#            GFC-2024-v1.12/download.html
# Citation : Hansen, M. C., et al. (2013). High-Resolution Global Maps of
#            21st-Century Forest Cover Change. Science, 342, 850–853.
#            https://glad.earthengine.app/view/global-forest-change
# License  : Creative Commons Attribution 4.0 International (CC BY 4.0)
#            Credit: Source: Hansen/UMD/Google/USGS/NASA
# Band     : "lossyear" — pixel value encodes the year of first loss
#            Value 1 = 2001, 2 = 2002, …, 24 = 2024 (0 = no loss)
# Spatial resolution: 30 m
# Format   : GeoTIFF (merged across tiles into single file)
# File     : source_data/forest_loss/forest_loss_merged_2000_2024.tif
#
# Administrative boundaries loaded from the database (populated by administrative.py).
#
# Spatial coverage  : Sulawesi and Maluku regions
# Temporal coverage : 2002–2024 (year values 2–24 in the raster)
#
# =============================================================
# PROCESSING
# =============================================================
# For each entity polygon:
#  1. Run zonal_stats with categorical=True to count pixels per year value.
#  2. For each year value y (2–24):
#       pixel_count = count of pixels with value y inside the polygon
#       area_ha     = pixel_count × (30²) / 10000
#                     (30 m pixel → 900 m² → 0.09 ha)
#
# CRS handling: if the vector CRS differs from the raster CRS, the
# vector layer is reprojected to match the raster before zonal_stats.
#
# =============================================================
# OUTPUT
# =============================================================
# output/table_forest_loss_province_<DD_MM_YYYY>.sql
# output/table_forest_loss_regency_city_<DD_MM_YYYY>.sql
#   INSERT INTO forest_loss (level, entity_id, entity_name, year,
#                            pixel_count, area_ha)
# =============================================================
"""

# =============================================================
# IMPORTS
# =============================================================
import rasterio
import geopandas as gpd
from rasterstats import zonal_stats
from tqdm import tqdm

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
RASTER_PATH             = "source_data/forest_loss/forest_loss_merged_2000_2024.tif"
OUTPUT_DIR              = "output"

# Raster year encoding: value 1=2001, 2=2002, …, 24=2024
# Analysis range: 2002–2024 → raster values 2–24
RASTER_YEARS   = list(range(2, 25))
PIXEL_SIZE_M   = 30             # 30 m spatial resolution
PIXEL_AREA_HA  = (PIXEL_SIZE_M ** 2) / 10_000   # 0.09 ha per pixel


def _ensure_crs(gdf: gpd.GeoDataFrame, raster_path: str) -> gpd.GeoDataFrame:
    """Reproject vector layer to match raster CRS if needed."""
    with rasterio.open(raster_path) as src:
        raster_crs = src.crs
    if gdf.crs != raster_crs:
        print(f"  ⓘ CRS mismatch — reprojecting from {gdf.crs} to {raster_crs}")
        gdf = gdf.to_crs(raster_crs)
    return gdf


# =============================================================
# SECTION 1 — PROVINCE LEVEL
# =============================================================

print("Processing province-level forest loss …")

con          = get_db_connection()
province_gdf = load_province_boundaries_from_db(con, regions=TARGET_REGIONS)
province_gdf = _ensure_crs(province_gdf, RASTER_PATH)

values_list = []

for _, row in tqdm(province_gdf.iterrows(), total=len(province_gdf), desc="Provinces"):
    zs = zonal_stats(
        [row.geometry],
        RASTER_PATH,
        stats=None,
        categorical=True,
        all_touched=True,
        nodata=0,
    )
    if not zs:
        continue

    cat_counts = zs[0]
    for y in RASTER_YEARS:
        count   = cat_counts.get(y, 0)
        area_ha = count * PIXEL_AREA_HA
        values_list.append(
            f"('province', '{row['province_id']}', '{row['province_name']}', "
            f"{2000 + y}, {count}, {area_ha:.2f})"
        )

sql_province = (
    "INSERT INTO forest_loss (level, entity_id, entity_name, year, "
    "pixel_count, area_ha) VALUES\n"
    + ",\n".join(values_list) + ";"
)
write_sql_file(sql_province, "forest_loss_province", OUTPUT_DIR)


# =============================================================
# SECTION 2 — REGENCY / CITY LEVEL
# =============================================================

print("Processing regency/city-level forest loss …")

regency_gdf = load_regency_city_boundaries_from_db(con, regions=TARGET_REGIONS)
con.close()
regency_gdf = _ensure_crs(regency_gdf, RASTER_PATH)

values_list = []

for _, row in tqdm(regency_gdf.iterrows(), total=len(regency_gdf), desc="Regencies/cities"):
    zs = zonal_stats(
        [row.geometry],
        RASTER_PATH,
        stats=None,
        categorical=True,
        all_touched=True,
        nodata=0,
    )
    if not zs:
        continue

    cat_counts = zs[0]
    for y in RASTER_YEARS:
        count   = cat_counts.get(y, 0)
        area_ha = count * PIXEL_AREA_HA
        values_list.append(
            f"('regency_city', '{row['regency_city_id']}', '{row['regency_city_name']}', "
            f"{2000 + y}, {count}, {area_ha:.2f})"
        )

sql_regency = (
    "INSERT INTO forest_loss (level, entity_id, entity_name, year, "
    "pixel_count, area_ha) VALUES\n"
    + ",\n".join(values_list) + ";"
)
write_sql_file(sql_regency, "forest_loss_regency_city", OUTPUT_DIR)
