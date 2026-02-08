import os

# Input directories (raw data)
ERA5_IN = os.getenv("ERA5_IN")
CMIP6_IN = os.getenv("CMIP6_IN")

# Processed data directories
ERA5_PROCESSED = os.getenv("ERA5_PROCESSED")
CMIP6_PROCESSED = os.getenv("CMIP6_PROCESSED")

# Output directory for bias correction and CFFDRS
OUT_DIR = os.getenv("OUT_DIR")

# Optional spatial mask
SHP_MASK = os.getenv("SHP_MASK")

# Clip hursmin to valid range [0.0, 100.0]
# set to FALSE for QC & comparison with older versions of the dataset
CLIP_HURSMIN = os.getenv("CLIP_HURSMIN", "TRUE").upper() == "TRUE"

# Legacy mode - skips time alignment and dimension transpose to match old pipeline behavior
# Set to TRUE to reproduce original pipeline results (for testing/comparison purposes only)
LEGACY_MODE = os.getenv("LEGACY_MODE", "FALSE").upper() == "TRUE"

era5_years = (1979, 2020)
cmip6_years = (1979, 2099)
hist_years = (1980, 2009)
sim_periods = ((2010, 2039), (2040, 2069), (2070, 2099))

geog_bbox = (-177.1875, 35.0, -35.0, 79.375)

gcm_list = ["EC-Earth3-Veg", "MPI-ESM1-2-HR", "MRI-ESM2-0", "CNRM-CM6-1-HR"]
metvars = ["tasmax", "pr", "sfcWind", "hursmin"]

qdm = {
    "oper": {
        "tasmax": "+",
        "pr": "*",
        "sfcWind": "*",
        "hursmin": "*",
    },
    "min_thresh": {
        "tasmax": None,
        "pr": 0.5,
        "sfcWind": None,
        "hursmin": None,
    },
    "quantile_vals": [
        0.005,
        0.01,
        0.02,
        0.03,
        0.04,
        0.05,
        0.06,
        0.07,
        0.08,
        0.09,
        0.1,
        0.11,
        0.12,
        0.13,
        0.14,
        0.15,
        0.16,
        0.17,
        0.18,
        0.19,
        0.2,
        0.21,
        0.22,
        0.23,
        0.24,
        0.25,
        0.26,
        0.27,
        0.28,
        0.29,
        0.3,
        0.31,
        0.32,
        0.33,
        0.34,
        0.35,
        0.36,
        0.37,
        0.38,
        0.39,
        0.4,
        0.41,
        0.42,
        0.43,
        0.44,
        0.45,
        0.46,
        0.47,
        0.48,
        0.49,
        0.5,
        0.51,
        0.52,
        0.53,
        0.54,
        0.55,
        0.56,
        0.57,
        0.58,
        0.59,
        0.6,
        0.61,
        0.62,
        0.63,
        0.64,
        0.65,
        0.66,
        0.67,
        0.68,
        0.69,
        0.7,
        0.71,
        0.72,
        0.73,
        0.74,
        0.75,
        0.76,
        0.77,
        0.78,
        0.79,
        0.8,
        0.81,
        0.82,
        0.83,
        0.84,
        0.85,
        0.86,
        0.87,
        0.88,
        0.89,
        0.9,
        0.91,
        0.92,
        0.93,
        0.94,
        0.95,
        0.96,
        0.97,
        0.98,
        0.99,
        0.995,
    ],
}
