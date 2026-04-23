"""
administrative.py — Administrative boundary processing.

# =============================================================
# DATA SOURCE
# =============================================================
# Province boundaries:
#   Source   : Badan Informasi Geospasial (BIG) — BATAS_WILAYAH MapServer
#   URL      : https://geoservices.big.go.id/rbi/rest/services/BATASWILAYAH/
#              BATAS_WILAYAH/MapServer (accessed 29 October 2025)
#   Format   : ESRI Shapefile (.shp), EPSG:4326
#   File     : source_data/administrative_province_border_2024.shp
#   Key cols : KDPPUM (province code), WADMPR (province name),
#              IBUKOTA (capital city), GEO_REGION (geographic region)
#
# Regency/city boundaries:
#   Source   : Badan Informasi Geospasial (BIG) — BATAS_WILAYAH MapServer
#   URL      : https://geoservices.big.go.id/rbi/rest/services/BATASWILAYAH/
#              BATAS_WILAYAH/MapServer (accessed 29 October 2025)
#   Format   : ESRI Shapefile (.shp), EPSG:4326
#   File     : source_data/administrative_regency_city_border_2024.shp
#   Key cols : KDPKAB (regency/city code), WADMKK (full name incl. type prefix),
#              KDPPUM (province code)
#              Type (Kabupaten/Kota) is derived from the WADMKK prefix.
#
# Spatial coverage : all Indonesian provinces and regencies/cities
# Temporal coverage: administrative boundaries as of 2024
#
# =============================================================
# OUTPUT
# =============================================================
# output/table_province_<DD_MM_YYYY>.sql
#   INSERT INTO province (id, name, name_clean, capital_city, geo_region, geometry)
#
# output/table_regency_city_<DD_MM_YYYY>.sql
#   INSERT INTO regency_city (id, province_id, type, name, name_clean, geometry)
#
# Geometries are stored as MySQL ST_GeomFromText WKT with SRID 4326.
# Z coordinates are dropped; X/Y are swapped to match MySQL's (lat, lon)
# convention for SRID 4326.
# =============================================================
"""

# =============================================================
# IMPORTS
# =============================================================
import re
import pandas as pd
import geopandas as gpd

from utils import (
    remove_z,
    swap_xy,
    clean_province_name,
    write_sql_file,
)

# =============================================================
# CONFIGURATION
# =============================================================
OUTPUT_DIR             = "output"
PROVINCE_SHAPEFILE     = "source_data/administrative_province_border_2024.shp"
REGENCY_CITY_SHAPEFILE = "source_data/administrative_regency_city_border_2024.shp"

ID_PATTERN = re.compile(r"^\d{2}\.\d{2}$")

# =============================================================
# SHARED HELPERS
# =============================================================

def _prep_geometry(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Drop Z coordinates and swap X/Y for MySQL ST_GeomFromText."""
    gdf = gdf.copy()
    gdf["geometry"] = gdf["geometry"].apply(remove_z)
    gdf["geometry"] = gdf["geometry"].apply(swap_xy)
    return gdf


def _clean_regency_name(s: str) -> str:
    """Normalise regency/city name for the name_clean column."""
    s = str(s).lower()
    s = re.sub(r"\b(kabupaten|kota|adm|administrasi|kep|kepulauan)\b", "", s)
    s = re.sub(r"[^\w]", "", s)
    s = re.sub(r"\s+",   "", s)
    return s


def _q(v) -> str:
    """Format a value as a SQL string literal or NULL."""
    if v is None:
        return "NULL"
    try:
        if pd.isna(v):
            return "NULL"
    except (TypeError, ValueError):
        pass
    return f"'{str(v).replace(chr(39), chr(39)*2)}'"


def _merge_and_drop(gdf: gpd.GeoDataFrame, wrong_mask, correct_mask, label: str) -> gpd.GeoDataFrame:
    """Merge geometries from wrong rows into correct row, then drop wrong rows."""
    if not wrong_mask.any():
        raise ValueError(f"No 'wrong' rows found: {label}")
    if not correct_mask.any():
        raise ValueError(f"No 'correct' row found: {label}")
    correct_idx = gdf[correct_mask].index[0]
    merged_geom = gdf.loc[correct_mask, "geometry"].iloc[0].union(
        gdf.loc[wrong_mask, "geometry"].union_all()
    )
    gdf.at[correct_idx, "geometry"] = merged_geom
    gdf = gdf[~wrong_mask]
    print(f"  ⓘ Fixed: {label}")
    return gdf


# =============================================================
# SECTION 1 — PROVINCE
# =============================================================

print("Processing province boundaries …")

gdf_province = gpd.read_file(PROVINCE_SHAPEFILE)

dupes = gdf_province[gdf_province.duplicated("KDPPUM", keep=False)]
print(f"  Shapefile duplicate KDPPUM: {len(dupes)}")

# Fix Papua Barat Daya (wrong code 92 → 96)
mask = gdf_province["WADMPR"].str.contains("Papua Barat Daya", case=False, na=False)
gdf_province.loc[mask, "KDPPUM"] = "96"
print("  ⓘ Papua Barat Daya KDPPUM corrected: 92 → 96")

dupes = gdf_province[gdf_province.duplicated("KDPPUM", keep=False)]
if not dupes.empty:
    raise ValueError(f"Unresolved duplicate KDPPUM:\n{dupes}")

gdf_province = gdf_province.sort_values("KDPPUM")
gdf_province["KDPPUM"]     = gdf_province["KDPPUM"].astype(str)
gdf_province["name_clean"] = gdf_province["WADMPR"].apply(clean_province_name)

print(f"  Provinces: {len(gdf_province)}")

gdf_province = _prep_geometry(gdf_province)

province_values = []
for _, row in gdf_province.iterrows():
    geom_wkt = row["geometry"].wkt if row["geometry"] is not None else None
    geom_sql = f"ST_GeomFromText('{geom_wkt}', 4326)" if geom_wkt else "NULL"
    province_values.append(
        f"('{row['KDPPUM']}', "
        f"'{row['WADMPR'].replace(chr(39), chr(39)*2)}', "
        f"'{row['name_clean']}', "
        + _q(row.get("IBUKOTA")) + ", "
        + _q(row.get("GEO_REGION")) + ", "
        + geom_sql + ")"
    )

sql_province = (
    "INSERT INTO province (id, name, name_clean, capital_city, geo_region, geometry) VALUES\n"
    + ",\n".join(province_values) + ";"
)
write_sql_file(sql_province, "province", OUTPUT_DIR)


# =============================================================
# SECTION 2 — REGENCY / CITY
# =============================================================

print("Processing regency/city boundaries …")

gdf = gpd.read_file(REGENCY_CITY_SHAPEFILE)

# Remove non-administrative rows
before = len(gdf)
gdf = gdf[gdf["KDPKAB"].notna() & (gdf["KDPKAB"] != "None")]
gdf = gdf[~gdf["KDPKAB"].str.match(r"^\d{2}\.\-\-$")]
print(f"  ⓘ Removed {before - len(gdf)} non-administrative island rows")

# Fix known duplicate and misallocated geometries
gdf = _merge_and_drop(
    gdf,
    wrong_mask=gdf["WADMPR"].str.contains("Papua Barat$", case=False, na=False) & (gdf["KDPKAB"] == "92.01"),
    correct_mask=gdf["WADMPR"].str.contains("Papua Barat Daya", case=False, na=False) & (gdf["KDPKAB"] == "92.01"),
    label="Sorong — Papua Barat entry merged into Papua Barat Daya",
)
gdf = _merge_and_drop(
    gdf,
    wrong_mask=gdf["WADMKK"].str.contains("Pahuwato", case=False, na=False),
    correct_mask=gdf["WADMKK"].str.contains("Pohuwato", case=False, na=False),
    label="Pohuwato — misspelled 'Pahuwato' merged into correct entry",
)
gdf = _merge_and_drop(
    gdf,
    wrong_mask=gdf["KDPKAB"] == "91.04",
    correct_mask=gdf["KDPKAB"] == "94.01",
    label="Nabire — wrong code 91.04 merged into 94.01",
)
gdf = _merge_and_drop(
    gdf,
    wrong_mask=gdf["KDPKAB"] == "91.08",
    correct_mask=gdf["KDPKAB"] == "94.03",
    label="Paniai — wrong code 91.08 merged into 94.03",
)
gdf = _merge_and_drop(
    gdf,
    wrong_mask=gdf["KDPKAB"] == "91.28",
    correct_mask=gdf["KDPKAB"] == "94.08",
    label="Deiyai — wrong code 91.28 merged into 94.08",
)
gdf = _merge_and_drop(
    gdf,
    wrong_mask=gdf["KDPKAB"] == "91.2",
    correct_mask=gdf["KDPKAB"] == "91.20",
    label="Mamberamo Raya — format 91.2 merged into 91.20",
)

# Fix invalid KDPKAB values from misallocated island polygons
gdf["KDPKAB"] = gdf["KDPKAB"].astype(str).str.strip()
invalid_codes = gdf[~gdf["KDPKAB"].str.match(ID_PATTERN)]
if not invalid_codes.empty:
    print(f"  ⓘ Invalid KDPKAB rows: {len(invalid_codes)}")

gdf = _merge_and_drop(
    gdf,
    wrong_mask=gdf["KDPKAB"].str.contains("52.0452.07", na=False),
    correct_mask=gdf["KDPKAB"] == "52.07",
    label="Gili Kalong — merged into Sumbawa Barat (52.07)",
)
gdf = _merge_and_drop(
    gdf,
    wrong_mask=gdf["KDPKAB"].str.contains("71.05/71.10", na=False),
    correct_mask=gdf["KDPKAB"] == "71.10",
    label="Danau Moat — merged into Bolaang Mongondow Timur (71.10)",
)

# Final validation
final_dupes = gdf[gdf.duplicated("KDPKAB", keep=False)]
if not final_dupes.empty:
    raise ValueError(
        f"Unresolved KDPKAB duplicates:\n"
        f"{final_dupes[['KDPKAB','WADMKK','WADMPR']].to_string(index=False)}"
    )

final_dupes = gdf[gdf.duplicated("WADMKK", keep=False)]
if not final_dupes.empty:
    raise ValueError(
        f"Unresolved WADMKK duplicates:\n"
        f"{final_dupes[['KDPKAB','WADMKK','WADMPR']].to_string(index=False)}"
    )

invalid_after = gdf[~gdf["KDPKAB"].str.match(ID_PATTERN)]
if not invalid_after.empty:
    raise ValueError(
        f"Still invalid KDPKAB codes:\n"
        f"{invalid_after[['KDPKAB','WADMKK']].to_string(index=False)}"
    )

# Derive type and clean name from WADMKK (e.g. "Kabupaten Gowa" → type="Kabupaten", name="Gowa")
gdf["type"]       = gdf["WADMKK"].str.extract(r"^(Kabupaten|Kota)")[0].fillna("Kabupaten")
gdf["name"]       = gdf["WADMKK"].str.replace(r"^(Kabupaten|Kota)\s+", "", regex=True)
gdf["name_clean"] = gdf["name"].apply(_clean_regency_name)

# Fix province IDs for Papua Barat Daya regencies
PAPUA_BARAT_DAYA_IDS = ["92.01", "92.04", "92.05", "92.09", "92.10", "92.71"]
gdf["KDPPUM"] = gdf["KDPPUM"].astype(str).str.strip()
gdf.loc[gdf["KDPKAB"].isin(PAPUA_BARAT_DAYA_IDS), "KDPPUM"] = "96"
print(f"  ⓘ Province IDs corrected to 96 for: {PAPUA_BARAT_DAYA_IDS}")

gdf = gdf.sort_values(["KDPPUM", "KDPKAB"])
print(f"  Regencies/cities: {len(gdf)}")

gdf = _prep_geometry(gdf)

regency_values = []
for _, row in gdf.iterrows():
    geom_wkt = row["geometry"].wkt if row["geometry"] is not None else None
    geom_sql = f"ST_GeomFromText('{geom_wkt}', 4326)" if geom_wkt else "NULL"
    regency_values.append(
        f"('{row['KDPKAB']}', "
        f"'{row['KDPPUM']}', "
        f"'{row['type'].replace(chr(39), chr(39)*2)}', "
        f"'{row['name'].replace(chr(39), chr(39)*2)}', "
        f"'{row['name_clean']}', "
        + geom_sql + ")"
    )

sql_regency = (
    "INSERT INTO regency_city (id, province_id, type, name, name_clean, geometry) VALUES\n"
    + ",\n".join(regency_values) + ";"
)
write_sql_file(sql_regency, "regency_city", OUTPUT_DIR)
