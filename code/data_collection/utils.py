"""
utils.py — Shared utilities for the data collection pipeline.

Provides:
  - SQL value formatting and INSERT statement construction
  - Dated SQL file writer
  - Geometry helpers: remove Z coordinate, swap X/Y for MySQL
  - Administrative boundary loaders (shapefile and database)
  - Database connection factory (MySQL via DATABASE/.env)
  - Chl-a pixel statistics aggregator

All scripts in this module import from here to avoid duplication.
Domain-specific logic (NDVI computation, BPS parsing, NTL outlier
detection, etc.) is NOT placed here — it stays in the domain scripts.
"""

import os
import re
import datetime
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.ops import transform
from shapely.wkt import loads as wkt_loads


# =============================================================
# REGION CONFIGURATION
# Mapping from geographic region name to province ID codes
# (Kemendagri format, matching shapefile KDPPUM column)
# =============================================================

GEO_REGIONS: dict[str, list[str]] = {
    "Sulawesi": ["71", "72", "73", "74", "75", "76"],
    "Maluku":   ["81", "82"],
}

TARGET_REGIONS: list[str] = ["Sulawesi", "Maluku"]


def region_province_ids(regions: list[str]) -> list[str]:
    """Return all province IDs belonging to the given region names."""
    return [pid for r in regions if r in GEO_REGIONS for pid in GEO_REGIONS[r]]


# =============================================================
# SQL FORMATTING
# =============================================================

def sql_str(val) -> str:
    """
    Format a scalar value for use inside a SQL INSERT statement.

      None / NaN / pd.NA           → NULL
      str                           → 'escaped string'
      float with no fractional part → integer string  (e.g. 5.0 → 5)
      other numeric                 → str(val)
    """
    if val is None:
        return "NULL"
    if isinstance(val, float) and np.isnan(val):
        return "NULL"
    try:
        if pd.isna(val):
            return "NULL"
    except (TypeError, ValueError):
        pass
    if isinstance(val, str):
        return "'" + val.replace("'", "''") + "'"
    if isinstance(val, (np.floating, float)) and float(val).is_integer():
        return str(int(val))
    return str(val)


def sql_num(val) -> str:
    """
    Like sql_str but maps None/NaN to '0' instead of NULL.
    Used for numeric columns where 0 is the correct missing sentinel
    (e.g. NTL pixel counts, forest loss area).
    """
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return "0"
    if isinstance(val, (np.floating, float)) and float(val).is_integer():
        return str(int(val))
    return str(val)


def build_insert(table: str, columns: list, rows: list) -> str:
    """
    Build a single SQL INSERT … VALUES statement.

    Parameters
    ----------
    table   : target table name
    columns : list of column name strings
    rows    : list of lists, each inner list is one row's values
              (values are formatted with sql_str)

    Returns the complete INSERT statement as a string.
    """
    col_str   = ", ".join(columns)
    val_lines = [
        "(" + ", ".join(sql_str(v) for v in row) + ")"
        for row in rows
    ]
    return f"INSERT INTO {table} ({col_str}) VALUES\n" + ",\n".join(val_lines) + ";"


def write_sql_file(sql: str, table_name: str, output_dir: str = "output") -> str:
    """
    Write a SQL string to a date-stamped file.

    File is written to:  <output_dir>/table_<table_name>_<DD_MM_YYYY>.sql

    Returns the absolute path of the written file.
    """
    os.makedirs(output_dir, exist_ok=True)
    today = datetime.date.today().strftime("%d_%m_%Y")
    path  = os.path.join(output_dir, f"table_{table_name}_{today}.sql")
    with open(path, "w", encoding="utf-8") as f:
        f.write(sql)
    print(f"Saved: {path}")
    return path


# =============================================================
# GEOMETRY HELPERS
# =============================================================

def remove_z(geom):
    """
    Drop the Z coordinate from any Shapely geometry type.
    Required before writing to MySQL ST_GeomFromText (2D only).
    """
    if geom is None or geom.is_empty:
        return geom
    return transform(lambda x, y, z=None: (x, y), geom)


def swap_xy(geom):
    """
    Swap X (longitude) and Y (latitude) coordinates.

    MySQL's ST_GeomFromText with SRID 4326 expects (lat, lon) order,
    while shapefiles and GeoPandas use (lon, lat) / (x, y) order.
    """
    if geom is None or geom.is_empty:
        return geom
    return transform(lambda x, y: (y, x), geom)


# =============================================================
# PROVINCE NAME CLEANING
# Used by administrative.py and gdp.py to normalise province names
# before string-matching across data sources.
# =============================================================

def clean_province_name(name: str) -> str:
    """
    Normalise a province name string for cross-source matching.

    Removes administrative prefixes (DKI, Istimewa, etc.), digits,
    punctuation, and whitespace, then lowercases.
    """
    s = str(name).lower()
    s = re.sub(r"\b(dki|istimewa|khusus|di |daerah|ibukota)\b", "", s)
    s = re.sub(r"\d+", "", s)
    s = re.sub(r"[^\w]", "", s)
    s = re.sub(r"\s+",  "", s)
    return s


# =============================================================
# ADMINISTRATIVE BOUNDARY LOADERS
# Replace DB geometry queries with direct shapefile reads.
# Shapefiles are the same source used by administrative.py.
# =============================================================

def load_province_boundaries(
    shapefile_path: str,
    regions: list = None,
) -> gpd.GeoDataFrame:
    """
    Load province boundaries from the administrative shapefile.

    Parameters
    ----------
    shapefile_path : relative path to the province .shp file
    regions        : list of geo_region names to keep (see GEO_REGIONS);
                     pass None to load all provinces

    Returns
    -------
    GeoDataFrame with columns:
      province_id   (str, Kemendagri code, e.g. "71")
      province_name (str)
      geometry      (EPSG:4326)
    """
    gdf = gpd.read_file(shapefile_path)
    gdf = gdf.rename(columns={"KDPPUM": "province_id", "WADMPR": "province_name"})
    gdf["province_id"] = gdf["province_id"].astype(str).str.strip()

    if regions is not None:
        target_ids = region_province_ids(regions)
        gdf = gdf[gdf["province_id"].isin(target_ids)]

    return (
        gdf[["province_id", "province_name", "geometry"]]
        .copy()
        .reset_index(drop=True)
    )


def load_regency_city_boundaries(
    shapefile_path: str,
    regions: list = None,
) -> gpd.GeoDataFrame:
    """
    Load regency/city boundaries from the administrative shapefile.

    Parameters
    ----------
    shapefile_path : relative path to the regency/city .shp file
    regions        : list of geo_region names to keep (see GEO_REGIONS);
                     pass None to load all regencies/cities

    Returns
    -------
    GeoDataFrame with columns:
      regency_city_id   (str, Kemendagri code, e.g. "71.01")
      regency_city_name (str)
      province_id       (str)
      geometry          (EPSG:4326)
    """
    gdf = gpd.read_file(shapefile_path)
    gdf = gdf.rename(columns={
        "KDPKAB": "regency_city_id",
        "WADMKK": "regency_city_name",
        "KDPPUM": "province_id",
    })
    gdf["regency_city_id"] = gdf["regency_city_id"].astype(str).str.strip()
    gdf["province_id"]     = gdf["province_id"].astype(str).str.strip()

    # Remove non-administrative placeholder rows
    gdf = gdf[gdf["regency_city_id"].notna() & (gdf["regency_city_id"] != "None")]
    gdf = gdf[~gdf["regency_city_id"].str.match(r"^\d{2}\.\-\-$")]

    if regions is not None:
        target_ids = region_province_ids(regions)
        gdf = gdf[gdf["province_id"].isin(target_ids)]

    return (
        gdf[["regency_city_id", "regency_city_name", "province_id", "geometry"]]
        .copy()
        .reset_index(drop=True)
    )


# =============================================================
# DATABASE CONNECTION
# Reads credentials from DATABASE/.env (two levels above this file's
# location: code/data_collection/ → code/ → project_root/ → DATABASE/).
# =============================================================

def get_db_connection():
    """
    Return an open mysql.connector connection using DATABASE/.env credentials.

    The .env file is resolved relative to this file's location:
      <project_root>/DATABASE/.env
    where <project_root> is four directories above this file.
    """
    try:
        from dotenv import load_dotenv
        import mysql.connector
    except ImportError as e:
        raise ImportError(
            "mysql-connector-python and python-dotenv are required for DB access. "
            f"Install with: pip install mysql-connector-python python-dotenv\n{e}"
        )

    env_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "..", "..", "DATABASE", ".env",
    )
    load_dotenv(os.path.normpath(env_path))

    import mysql.connector
    return mysql.connector.connect(
        host=os.getenv("LOCAL_HOST", "127.0.0.1"),
        port=int(os.getenv("LOCAL_PORT", "3306")),
        user=os.getenv("LOCAL_USER"),
        password=os.getenv("LOCAL_PASSWORD"),
        database=os.getenv("LOCAL_DATABASE"),
    )


def load_province_boundaries_from_db(con, regions: list = None) -> gpd.GeoDataFrame:
    """
    Load province boundaries from the database (populated by administrative.py).

    Geometries are stored in MySQL SRID 4326 with swapped axis order (lat, lon).
    They are swapped back to (lon, lat) for shapely/geopandas compatibility.

    Parameters
    ----------
    con     : open mysql.connector connection
    regions : list of geo_region names to filter by (see GEO_REGIONS);
              pass None to load all provinces

    Returns
    -------
    GeoDataFrame with columns: province_id, province_name, geometry (EPSG:4326)
    """
    cursor = con.cursor(dictionary=True)

    if regions is not None:
        target_ids = region_province_ids(regions)
        placeholders = ", ".join(["%s"] * len(target_ids))
        cursor.execute(
            f"SELECT id AS province_id, name AS province_name, "
            f"ST_AsText(geometry) AS geom_wkt "
            f"FROM province WHERE id IN ({placeholders})",
            target_ids,
        )
    else:
        cursor.execute(
            "SELECT id AS province_id, name AS province_name, "
            "ST_AsText(geometry) AS geom_wkt FROM province"
        )

    rows = cursor.fetchall()
    cursor.close()

    records = []
    for row in rows:
        geom = swap_xy(wkt_loads(row["geom_wkt"])) if row["geom_wkt"] else None
        records.append({
            "province_id":   row["province_id"],
            "province_name": row["province_name"],
            "geometry":      geom,
        })

    return gpd.GeoDataFrame(records, crs="EPSG:4326").reset_index(drop=True)


def load_regency_city_boundaries_from_db(con, regions: list = None) -> gpd.GeoDataFrame:
    """
    Load regency/city boundaries from the database (populated by administrative.py).

    Parameters
    ----------
    con     : open mysql.connector connection
    regions : list of geo_region names to filter by (see GEO_REGIONS);
              pass None to load all regencies/cities

    Returns
    -------
    GeoDataFrame with columns:
      regency_city_id, regency_city_name, province_id, geometry (EPSG:4326)
    """
    cursor = con.cursor(dictionary=True)

    if regions is not None:
        target_ids = region_province_ids(regions)
        placeholders = ", ".join(["%s"] * len(target_ids))
        cursor.execute(
            f"SELECT id AS regency_city_id, name AS regency_city_name, province_id, "
            f"ST_AsText(geometry) AS geom_wkt "
            f"FROM regency_city WHERE province_id IN ({placeholders})",
            target_ids,
        )
    else:
        cursor.execute(
            "SELECT id AS regency_city_id, name AS regency_city_name, province_id, "
            "ST_AsText(geometry) AS geom_wkt FROM regency_city"
        )

    rows = cursor.fetchall()
    cursor.close()

    records = []
    for row in rows:
        geom = swap_xy(wkt_loads(row["geom_wkt"])) if row["geom_wkt"] else None
        records.append({
            "regency_city_id":   row["regency_city_id"],
            "regency_city_name": row["regency_city_name"],
            "province_id":       row["province_id"],
            "geometry":          geom,
        })

    return gpd.GeoDataFrame(records, crs="EPSG:4326").reset_index(drop=True)


# =============================================================
# CHL-A STATISTICS HELPER
# Shared between chl_a province and regency/city processing.
# =============================================================

def finalize_stats(pixel_lists: list) -> tuple:
    """
    Aggregate a list of 1-D pixel arrays into summary statistics.

    Parameters
    ----------
    pixel_lists : list of numpy arrays (one array per input file/timestep)

    Returns
    -------
    (sum, mean, median, max, amplitude) — all rounded; zeros on empty input
    """
    if not pixel_lists:
        return 0.0, np.nan, np.nan, np.nan, np.nan
    allpix = np.concatenate(pixel_lists)
    if allpix.size == 0:
        return 0.0, np.nan, np.nan, np.nan, np.nan
    allpix = allpix[~np.isnan(allpix)]
    if allpix.size == 0:
        return 0.0, np.nan, np.nan, np.nan, np.nan
    s   = round(float(np.sum(allpix)),             2)
    m   = round(float(np.mean(allpix)),            4)
    md  = round(float(np.median(allpix)),          4)
    mx  = round(float(np.max(allpix)),             4)
    amp = round(float(np.max(allpix) - np.min(allpix)), 4)
    return s, m, md, mx, amp
