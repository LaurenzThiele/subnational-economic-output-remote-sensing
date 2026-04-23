-- =============================================================
-- Database schema for the subnational economic output project
-- MySQL 8.0  |  SRID 4326 (WGS 84)
--
-- Execution order:
--   1. This file (schema.sql)
--   2. output/table_province_*.sql
--   3. output/table_regency_city_*.sql
--   4. output/table_gdp_*.sql
--   5. output/table_ndvi_*.sql
--   6. output/table_forest_loss_*.sql
--   7. output/table_ntl_*.sql
--   8. output/table_chl_a_*.sql
-- =============================================================

-- -------------------------------------------------------------
-- 1. province
--    Populated by administrative.py
-- -------------------------------------------------------------
CREATE TABLE IF NOT EXISTS province (
    id           VARCHAR(10)  NOT NULL,
    name         VARCHAR(120) NOT NULL,
    name_clean   VARCHAR(120),
    capital_city VARCHAR(120),
    geo_region   VARCHAR(60),
    geometry     GEOMETRY     NOT NULL SRID 4326,
    PRIMARY KEY (id),
    SPATIAL INDEX (geometry)
);

-- -------------------------------------------------------------
-- 2. regency_city
--    Populated by administrative.py
-- -------------------------------------------------------------
CREATE TABLE IF NOT EXISTS regency_city (
    id          VARCHAR(10)  NOT NULL,
    province_id VARCHAR(10)  NOT NULL,
    type        VARCHAR(20),
    name        VARCHAR(120) NOT NULL,
    name_clean  VARCHAR(120),
    geometry    GEOMETRY     NOT NULL SRID 4326,
    PRIMARY KEY (id),
    SPATIAL INDEX (geometry),
    FOREIGN KEY (province_id) REFERENCES province(id)
);

-- -------------------------------------------------------------
-- 3. gdp
--    Gross Regional Domestic Product (GRDP) by sector.
--    All monetary values in billion IDR (current prices).
--    base_year: 2000 or 2010 (chain-linked series).
--    Populated by gdp.py
-- -------------------------------------------------------------
CREATE TABLE IF NOT EXISTS gdp (
    level                              VARCHAR(20)    NOT NULL,
    entity_id                          VARCHAR(10)    NOT NULL,
    entity_name                        VARCHAR(120)   NOT NULL,
    year                               SMALLINT       NOT NULL,
    base_year                          SMALLINT       NOT NULL,
    agriculture_forestry_fisheries     DOUBLE,
    mining                             DOUBLE,
    manufacturing                      DOUBLE,
    electricity_gas_supply             DOUBLE,
    water_waste_recycling              DOUBLE,
    construction                       DOUBLE,
    wholesale_retail_trade_repair      DOUBLE,
    transportation_warehousing         DOUBLE,
    accommodation_food_services        DOUBLE,
    information_communication          DOUBLE,
    finance_insurance                  DOUBLE,
    real_estate                        DOUBLE,
    business_services                  DOUBLE,
    government_defense_social_security DOUBLE,
    education_services                 DOUBLE,
    health_social_services             DOUBLE,
    other_services                     DOUBLE,
    gross_domestic_product             DOUBLE,
    PRIMARY KEY (level, entity_id, year, base_year)
);

-- -------------------------------------------------------------
-- 4. ndvi
--    Annual NDVI statistics from MODIS MOD13A1.
--    Only Evergreen Broadleaf + Cropland pixels (IGBP 2/12/14).
--    Savitzky-Golay smoothed before aggregation.
--    Populated by ndvi.py
-- -------------------------------------------------------------
CREATE TABLE IF NOT EXISTS ndvi (
    level       VARCHAR(20)   NOT NULL,
    entity_id   VARCHAR(10)   NOT NULL,
    entity_name VARCHAR(120)  NOT NULL,
    year        SMALLINT      NOT NULL,
    sum         DOUBLE,
    mean        DECIMAL(10,4),
    median      DECIMAL(10,4),
    max         DECIMAL(10,4),
    amplitude   DECIMAL(10,4),
    PRIMARY KEY (level, entity_id, year)
);

-- -------------------------------------------------------------
-- 5. forest_loss
--    Annual forest cover loss area from Hansen GFC v1.12.
--    30 m resolution; pixel_count × 0.09 ha per pixel.
--    Populated by forest_loss.py
-- -------------------------------------------------------------
CREATE TABLE IF NOT EXISTS forest_loss (
    level       VARCHAR(20)   NOT NULL,
    entity_id   VARCHAR(10)   NOT NULL,
    entity_name VARCHAR(120)  NOT NULL,
    year        SMALLINT      NOT NULL,
    pixel_count INT           NOT NULL DEFAULT 0,
    area_ha     DECIMAL(14,2) NOT NULL DEFAULT 0.00,
    PRIMARY KEY (level, entity_id, year)
);

-- -------------------------------------------------------------
-- 6. ntl
--    Annual Night-Time Light (NTL) statistics.
--    Outlier-corrected via 5-year pixel-level moving median + MAD.
--    DN threshold for filtered series: DN > 7.
--    Populated by ntl.py
-- -------------------------------------------------------------
CREATE TABLE IF NOT EXISTS ntl (
    level                        VARCHAR(20)  NOT NULL,
    entity_id                    VARCHAR(10)  NOT NULL,
    entity_name                  VARCHAR(120) NOT NULL,
    year                         SMALLINT     NOT NULL,
    data_pixels                  INT          DEFAULT 0,
    data_pixels_greater_than_7   INT          DEFAULT 0,
    data_pixels_greater_than_0   INT          DEFAULT 0,
    sum                          DOUBLE       DEFAULT 0,
    imputed_outliers             INT          DEFAULT 0,
    sum_unfiltered               DOUBLE       DEFAULT 0,
    imputed_outliers_unfiltered  INT          DEFAULT 0,
    PRIMARY KEY (level, entity_id, year)
);

-- -------------------------------------------------------------
-- 7. chl_a
--    Annual chlorophyll-a statistics for three coastal zone combos:
--      B100_D5  : depth ≥ −100 m, buffer ≤  5 km
--      B150_D15 : depth ≥ −150 m, buffer ≤ 15 km
--      B200_D25 : depth ≥ −200 m, buffer ≤ 25 km
--    Savitzky-Golay smoothed along time axis before aggregation.
--    Populated by chl_a.py
-- -------------------------------------------------------------
CREATE TABLE IF NOT EXISTS chl_a (
    level                VARCHAR(20)   NOT NULL,
    entity_id            VARCHAR(10)   NOT NULL,
    entity_name          VARCHAR(120)  NOT NULL,
    year                 SMALLINT      NOT NULL,
    -- B100_D5
    sum_B100_D5          DOUBLE,
    mean_B100_D5         DECIMAL(12,4),
    median_B100_D5       DECIMAL(12,4),
    max_B100_D5          DECIMAL(12,4),
    amplitude_B100_D5    DECIMAL(12,4),
    -- B150_D15
    sum_B150_D15         DOUBLE,
    mean_B150_D15        DECIMAL(12,4),
    median_B150_D15      DECIMAL(12,4),
    max_B150_D15         DECIMAL(12,4),
    amplitude_B150_D15   DECIMAL(12,4),
    -- B200_D25
    sum_B200_D25         DOUBLE,
    mean_B200_D25        DECIMAL(12,4),
    median_B200_D25      DECIMAL(12,4),
    max_B200_D25         DECIMAL(12,4),
    amplitude_B200_D25   DECIMAL(12,4),
    PRIMARY KEY (level, entity_id, year)
);
