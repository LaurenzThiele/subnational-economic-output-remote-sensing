"""
gdp.py — Gross Regional Domestic Product (GRDP) processing.

# =============================================================
# DATA SOURCE
# =============================================================
# Province level (two base-year series):
#   Base year 2000 (2002–2009):
#     Source : BPS Indonesia — Statistics Table
#     URL    : https://www.bps.go.id/en/statistics-table/1/MTYyNyMx/...
#     Format : CSV, header on row 2, years as columns
#     File   : source_data/gdp/province_level_2000_2013_base_year_2000.csv
#
#   Base year 2010 (2010–2024):
#     Source : BPS Indonesia — Statistics Table
#     URL    : https://www.bps.go.id/id/statistics-table/2/MjI2NyMy/...
#     Format : CSV per year, filename contains the year, header on row 1
#     Files  : source_data/gdp/*base_year_2010*.csv
#
# Regency/city level (base year 2010, 2015–2024):
#   Source : BPS Interoperability API (SIMDASI)
#   URL    : https://webapi.bps.go.id/v1/api/interoperabilitas/datasource/simdasi/
#             id/25/tahun/{year}/id_tabel/UklLSnFZZnMzMlJiSWpMOExJODIrQT09/
#             wilayah/{regency_city_id}/key/{api_key}
#   Auth   : API key passed via BPS_API_KEY environment variable
#   Format : JSON response with nested sector labels and values
#
# Administrative reference:
#   Province and regency/city IDs and names are loaded from the database
#   (populated by administrative.py).
#
# Spatial coverage : Sulawesi and Maluku regions
# Temporal coverage:
#   Province    : 2002–2024 (two overlapping base-year series)
#   Regency/city: 2015–2024
#
# =============================================================
# OUTPUT
# =============================================================
# output/table_gdp_province_<DD_MM_YYYY>.sql
#   INSERT INTO gdp (level, entity_id, entity_name, year, base_year,
#     agriculture_forestry_fisheries, mining, manufacturing,
#     electricity_gas_supply, water_waste_recycling, construction,
#     wholesale_retail_trade_repair, transportation_warehousing,
#     accommodation_food_services, information_communication,
#     finance_insurance, real_estate, business_services,
#     government_defense_social_security, education_services,
#     health_social_services, other_services, gross_domestic_product)
#
# output/table_gdp_regency_city_<DD_MM_YYYY>.sql
#   Same columns as above (base_year = 2010 for all regency rows)
#
# Intermediate CSVs cached per year:
#   source_data/gdp/regency_city_{year}.csv
# =============================================================
"""

# =============================================================
# IMPORTS
# =============================================================
import os
import re
import sys
import glob
import time
import datetime
import numpy as np
import pandas as pd
import requests
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

from utils import (
    clean_province_name,
    write_sql_file,
    get_db_connection,
    load_province_boundaries_from_db,
    load_regency_city_boundaries_from_db,
    TARGET_REGIONS,
)

# =============================================================
# CONFIGURATION
# =============================================================
GDP_SOURCE_DIR          = "source_data/gdp"
OUTPUT_DIR              = "output"

BPS_API_KEY     = os.getenv("BPS_API_KEY")   # set in environment; not stored in code
REGENCY_YEARS   = list(range(2015, 2025))
PROVINCE_YEARS  = list(range(2002, 2025))

# BPS API URL template for regency/city GRDP
BPS_URL_TEMPLATE = (
    "https://webapi.bps.go.id/v1/api/interoperabilitas/datasource/simdasi/"
    "id/25/tahun/{year}/id_tabel/UklLSnFZZnMzMlJiSWpMOExJODIrQT09/"
    "wilayah/{regency_city_id}/key/{api_key}"
)

# GDP table column order (matches INSERT statement)
GDP_COLUMNS = [
    "level", "entity_id", "entity_name", "year", "base_year",
    "agriculture_forestry_fisheries", "mining", "manufacturing",
    "electricity_gas_supply", "water_waste_recycling", "construction",
    "wholesale_retail_trade_repair", "transportation_warehousing",
    "accommodation_food_services", "information_communication",
    "finance_insurance", "real_estate", "business_services",
    "government_defense_social_security", "education_services",
    "health_social_services", "other_services", "gross_domestic_product",
]


def _sql_val(val) -> str:
    """Format a GDP sector value: NULL for NaN, else bare numeric string."""
    if pd.isna(val):
        return "NULL"
    return str(val)


# =============================================================
# SECTION 1 — PROVINCE LEVEL
# =============================================================

print("Processing province-level GDP …")

# Load province reference from database (ID + name for join)
con          = get_db_connection()
province_ref = load_province_boundaries_from_db(con, regions=TARGET_REGIONS)
province_ref["province_name_clean"] = province_ref["province_name"].apply(clean_province_name)

# --- 1.1  Base year 2010 CSVs --------------------------------
# One CSV per year; filename encodes the year.
# Row structure: rows 0–2 are headers; data starts at row 3.
# "Tahunan" (annual) columns are identified from row 2.

df_2010_list = []

for file in glob.glob(f"{GDP_SOURCE_DIR}/*base_year_2010*.csv"):
    match = re.search(r"(\d{4})", file)
    year  = int(match.group(1)) if match else None

    df = pd.read_csv(file, header=1)

    # Identify "Tahunan" (annual total) column indices from row 2
    tahunan_cols = [i for i in range(df.shape[1]) if "Tahunan" in str(df.iloc[2, i])]

    # Walk left from each Tahunan column to find its sector label (row 0)
    sector_labels = []
    for col in tahunan_cols:
        j = col
        while j >= 0 and pd.isna(df.iloc[0, j]):
            j -= 1
        sector_labels.append(df.iloc[0, j])

    keep_cols  = [0] + tahunan_cols
    tahunan_df = df.iloc[3:, keep_cols].reset_index(drop=True)
    tahunan_df.columns = ["Province"] + sector_labels

    tahunan_df = tahunan_df.dropna(how="all")
    tahunan_df["province_name_clean"] = tahunan_df["Province"].apply(clean_province_name)
    tahunan_df = tahunan_df.replace("--", np.nan)

    merged = province_ref.merge(tahunan_df, on="province_name_clean", how="left")
    merged = merged.drop(columns=["Province", "geometry"], errors="ignore")
    merged["year"]      = year
    merged["base_year"] = 2010
    df_2010_list.append(merged)

df_2010_all = pd.concat(df_2010_list, ignore_index=True)

# --- 1.2  Base year 2000 CSV ---------------------------------
# Single CSV covering 2000–2013; years are column names.

file_2000 = f"{GDP_SOURCE_DIR}/province_level_2000_2013_base_year_2000.csv"

df_2000 = pd.read_csv(file_2000, header=2).dropna(how="all")
df_2000["province_name_clean"] = df_2000["Province"].apply(clean_province_name)
df_2000 = df_2000.replace("--", np.nan)

merged_2000 = province_ref.merge(df_2000, on="province_name_clean", how="left")
merged_2000 = merged_2000.drop(columns=["Province", "geometry"], errors="ignore")
merged_2000["base_year"] = 2000

df_2000_long = merged_2000.melt(
    id_vars=["province_id", "province_name", "province_name_clean", "base_year"],
    var_name="year",
    value_name="Produk Domestik Regional Bruto",
)
df_2000_long["year"] = df_2000_long["year"].astype(int)

# --- 1.3  Combine and build SQL ------------------------------

df_gdp_province = pd.concat([df_2010_all, df_2000_long], ignore_index=True)

province_sector_cols_2010 = [
    "A Pertanian, Kehutanan dan Perikanan",
    "B Pertambangan dan Penggalian",
    "C Industri Pengolahan",
    "D Pengadaan Listrik dan Gas",
    "E Pengadaan Air, Pengelolaan Sampah, Limbah dan Daur Ulang",
    "F Konstruksi",
    "G Perdagangan Besar dan Eceran, Reparasi Mobil dan Sepeda Motor",
    "H Transportasi dan Pergudangan",
    "I Penyediaan Akomodasi dan Makan Minum",
    "J Informasi dan Komunikasi",
    "K Jasa Keuangan dan Asuransi",
    "L Real Estate",
    "M,N Jasa Perusahaan",
    "O Administrasi Pemerintahan, Pertahanan dan Jaminan Sosial Wajib",
    "P Jasa Pendidikan",
    "Q Jasa Kesehatan dan Kegiatan Sosial",
    "R,S,T,U Jasa Lainnya",
    "Produk Domestik Regional Bruto",
]

values_province = []
for _, row in df_gdp_province.iterrows():
    values_province.append(
        f"('province', '{row['province_id']}', '{row['province_name']}', "
        f"{row['year']}, {row['base_year']}, "
        + ", ".join(_sql_val(row.get(c)) for c in province_sector_cols_2010)
        + ")"
    )

sql_province = (
    "INSERT INTO gdp (level, entity_id, entity_name, year, base_year, "
    "agriculture_forestry_fisheries, mining, manufacturing, "
    "electricity_gas_supply, water_waste_recycling, construction, "
    "wholesale_retail_trade_repair, transportation_warehousing, "
    "accommodation_food_services, information_communication, "
    "finance_insurance, real_estate, business_services, "
    "government_defense_social_security, education_services, "
    "health_social_services, other_services, gross_domestic_product) VALUES\n"
    + ",\n".join(values_province) + ";"
)

write_sql_file(sql_province, "gdp_province", OUTPUT_DIR)


# =============================================================
# SECTION 2 — REGENCY / CITY LEVEL
# =============================================================

print("Processing regency/city-level GDP …")

# Load regency/city IDs from database
regency_ref = load_regency_city_boundaries_from_db(con, regions=TARGET_REGIONS)
con.close()


def _format_bps_id(db_id: str) -> int:
    """
    Convert Kemendagri code to BPS numeric ID format.

    BPS and Kemendagri codes differ for some regencies; known mismatches
    are corrected before conversion.

    Conversion: "73.24" → int("73" + "24000") = 7324000
    """
    if db_id == "73.24":   # Luwu Timur — BPS uses 73.25
        db_id = "73.25"
    if isinstance(db_id, str) and "." in db_id:
        parts = db_id.split(".")
        return int(parts[0] + parts[1].ljust(5, "0"))
    return int(db_id)


regency_ref["bps_id"] = regency_ref["regency_city_id"].apply(_format_bps_id)


def _fetch_bps_gdp(bps_id: int, db_id: str, name: str, year: int,
                   retries: int = 3, wait: int = 30):
    """
    Fetch GRDP data for one regency/city and year from the BPS API.

    Returns (DataFrame, label_order) on success, (None, None) on failure.
    """
    url = BPS_URL_TEMPLATE.format(year=year, regency_city_id=bps_id, api_key=BPS_API_KEY)

    for attempt in range(retries):
        try:
            resp = requests.get(url, timeout=20)
            if resp.status_code != 200:
                return None, None
            break
        except (requests.exceptions.ConnectTimeout, requests.exceptions.ReadTimeout):
            if attempt < retries - 1:
                time.sleep(wait)
                continue
            return None, None
        except requests.exceptions.RequestException as e:
            print(f"Request failed for {url}: {e}")
            return None, None

    data = resp.json()
    if data.get("status") != "OK":
        return None, None

    rows = data.get("data", [])
    if any(isinstance(r, dict) and r.get("message") == "Tabel tidak ditemukan" for r in rows):
        return None, None
    if any(isinstance(r, dict) and r.get("message") == "Wilayah tidak ditemukan" for r in rows):
        print(f"Region not found: {url}")
        sys.exit(1)
    if not rows or len(rows) < 2:
        return None, None

    table      = rows[1]
    kolom_key  = list(table["kolom"].keys())[0]
    records    = []
    labels     = []

    for entry in table.get("data", []):
        val_raw = entry["variables"][kolom_key]["value_raw"]
        value   = float(val_raw.replace(".", "").replace(",", ".")) if val_raw is not None else 0.0
        records.append({
            "regency_city_id":   db_id,
            "regency_city_name": name,
            "label":             entry.get("label_raw"),
            "value":             value,
        })
        labels.append(entry.get("label_raw"))

    df = pd.DataFrame(records)
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    return df, labels


# Parallel fetch per year (with CSV caching)
os.makedirs(GDP_SOURCE_DIR, exist_ok=True)
data_by_year = []

for year in REGENCY_YEARS:
    cache_path = f"{GDP_SOURCE_DIR}/regency_city_{year}.csv"

    if os.path.exists(cache_path):
        data_by_year.append(
            pd.read_csv(cache_path, dtype={"regency_city_id": str})
        )
        continue

    all_dfs = []
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = {
            executor.submit(
                _fetch_bps_gdp, bps_id, db_id, name, year
            ): (bps_id, db_id, name)
            for bps_id, db_id, name in regency_ref[
                ["bps_id", "regency_city_id", "regency_city_name"]
            ].values
        }
        for f in tqdm(as_completed(futures), total=len(futures), desc=f"Fetching {year}"):
            df, _ = f.result()
            if df is not None:
                all_dfs.append(df)

    if not all_dfs:
        print(f"No data for year {year}")
        continue

    year_df = pd.concat(all_dfs, ignore_index=True)
    wide_df = year_df.pivot_table(
        index=["regency_city_id", "regency_city_name"],
        columns="label",
        values="value",
        aggfunc="first",
    ).reset_index()
    cols    = ["regency_city_id", "regency_city_name"] + [
        c for c in wide_df.columns if c not in ["regency_city_id", "regency_city_name"]
    ]
    wide_df = wide_df[cols]
    data_by_year.append(wide_df)
    wide_df.to_csv(cache_path, index=False)
    print(f"Saved {cache_path} ({len(wide_df)} rows)")

# Build SQL INSERT
regency_sector_cols = [
    "A Pertanian, Kehutanan, dan Perikanan",
    "B Pertambangan dan Penggalian",
    "C Industri Pengolahan",
    "D Pengadaan Listrik dan Gas",
    "E Pengadaan Air; Pengelolaan Sampah, Limbah, dan Daur Ulang",
    "F Konstruksi",
    "G Perdagangan Besar dan Eceran; Reparasi Mobil dan Sepeda Motor",
    "H Transportasi dan Pergudangan",
    "I Penyediaan Akomodasi dan Makan Minum",
    "J Informasi dan Komunikasi",
    "K Jasa Keuangan dan Asuransi",
    "L Real Estat",
    "M,N Jasa Perusahaan",
    "O Administrasi Pemerintahan, Pertahanan, dan Jaminan Sosial Wajib",
    "P Jasa Pendidikan",
    "Q Jasa Kesehatan dan Kegiatan Sosial",
    "R,S,T,U Jasa Lainnya",
    "Produk Domestik Bruto",
]

values_regency = []
for year, year_df in zip(REGENCY_YEARS, data_by_year):
    for _, row in year_df.iterrows():
        values_regency.append(
            f"('regency_city', '{row['regency_city_id']}', '{row['regency_city_name']}', "
            f"{year}, 2010, "
            + ", ".join(f"'{row[c]}'" for c in regency_sector_cols)
            + ")"
        )

sql_regency = (
    "INSERT INTO gdp (level, entity_id, entity_name, year, base_year, "
    "agriculture_forestry_fisheries, mining, manufacturing, "
    "electricity_gas_supply, water_waste_recycling, construction, "
    "wholesale_retail_trade_repair, transportation_warehousing, "
    "accommodation_food_services, information_communication, "
    "finance_insurance, real_estate, business_services, "
    "government_defense_social_security, education_services, "
    "health_social_services, other_services, gross_domestic_product) VALUES\n"
    + ",\n".join(values_regency) + ";"
)

write_sql_file(sql_regency, "gdp_regency_city", OUTPUT_DIR)
