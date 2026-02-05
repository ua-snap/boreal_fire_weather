# Bias Correction of CMIP6 GCMs and Calculation of CFFDRS Fire Weather Indices

Run the following pipeline to bias correct CMIP6 data variables using ERA5, then calculate CFFDRS fire weather indices from the results. 

This README assumes you have downloaded copies of the CMIP6 and ERA5 source data. Contact SNAP for information about obtaining source data for this processing pipeline.

### Setup

Create a micromamba environment using the `environment.yml`
```bash
micromamba env create -f environment.yml
micromamba activate cffdrs
```

Set environment variables
```bash
# Input directories (raw data)
export ERA5_IN=/path/to/era5/input/data
export CMIP6_IN=/path/to/cmip6/input/data

# Processed data directories (output of 01_*.py and 02_*.py scripts)
export ERA5_PROCESSED=/path/to/era5/processed
export CMIP6_PROCESSED=/path/to/cmip6/processed

# Output directory for bias correction and CFFDRS
export OUT_DIR=/path/for/output/files
```

Optionally, include a shapefile to subset the data during the bias correction step
```bash
export SHP_MASK=/path/to/shapefile/mask.shp
```

**Note:** If your data is already processed (daily summaries with correct variable names and calendar), you can skip scripts `01_*.py` and `02_*.py` and set `ERA5_PROCESSED` and `CMIP6_PROCESSED` to point directly to your existing processed data.

### Directory Structure

The pipeline uses the following directory structure:
```
ERA5_IN/                          # Raw ERA5 input data (hourly)
├── t2m_1979.nc                   # Temperature files
├── tp_1979.nc                    # Precipitation files
├── u10_1979.nc                   # Wind u-component files
├── v10_1979.nc                   # Wind v-component files
├── d2m_1979.nc                   # Dewpoint files
└── ...

CMIP6_IN/                         # Raw CMIP6 input data (daily)
├── CNRM-CM6-1-HR/
├── EC-Earth3-Veg/
├── MPI-ESM1-2-HR/
└── MRI-ESM2-0/

ERA5_PROCESSED/                   # Processed ERA5 (daily summaries)
├── tasmax_era5_1979.nc
├── pr_era5_1979.nc
├── sfcWind_era5_1979.nc
├── hursmin_era5_1979.nc
└── ...

CMIP6_PROCESSED/                  # Processed CMIP6 (standardized)
├── CNRM-CM6-1-HR/
│   ├── tasmax_CNRM-CM6-1-HR_1979.nc
│   └── ...
├── EC-Earth3-Veg/
├── MPI-ESM1-2-HR/
└── MRI-ESM2-0/

OUT_DIR/                          # Final outputs
├── bias_corrected/               # Outputs from 03_*.py
│   ├── CNRM-CM6-1-HR/
│   ├── EC-Earth3-Veg/
│   ├── MPI-ESM1-2-HR/
│   └── MRI-ESM2-0/
└── cffdrs/                       # Outputs from 04_*.py
    ├── era5/
    ├── CNRM-CM6-1-HR/
    ├── EC-Earth3-Veg/
    ├── MPI-ESM1-2-HR/
    └── MRI-ESM2-0/
```

### Pipeline Scripts

#### 01_process_era5.py
Process ERA5 hourly data to daily summaries:
- Calculate daily maximum temperature (`tasmax`) from hourly temperature
- Calculate daily total precipitation (`pr`) from hourly precipitation
- Calculate daily mean wind speed (`sfcWind`) from `u10`/`v10` components
- Calculate daily minimum relative humidity (`hursmin`) from temperature and dewpoint
- Modify geographic extent to region of interest
- Standardize calendar to 'noleap' (365-day years)
- Input: `ERA5_IN/`
- Output: `ERA5_PROCESSED/`

#### 02_process_cmip6.py
Process CMIP6 daily data:
- Modify geographic extent to match ERA5
- Standardize calendar to 'noleap'
- Convert units (e.g., K → degC, m/s → km/hour)
- Calculate wind speed from u/v components where needed
- Input: `CMIP6_IN/{gcm}/`
- Output: `CMIP6_PROCESSED/{gcm}/`

#### 03_bias_correct_gcms.py
Apply Quantile Delta Mapping (QDM) bias correction to CMIP6 data using ERA5 as reference:
- Train QDM model using historical overlap period (1980-2009)
- Apply bias correction to historical GCM data (1980-2009)
- Apply bias correction to future projection periods (2010-2039, 2040-2069, 2070-2099)
- Optional: Apply spatial mask from shapefile
- Handle precipitation thresholds with jitter method
- Output: `OUT_DIR/bias_corrected/{gcm}/`

#### 04_calculate_cffdrs.py
Calculate Canadian Forest Fire Danger Rating System (CFFDRS) indices:
- Fine Fuel Moisture Code (FFMC)
- Duff Moisture Code (DMC)
- Drought Code (DC)
- Initial Spread Index (ISI)
- Build-up Index (BUI)
- Fire Weather Index (FWI)
- Daily Severity Rating (DSR)

Calculated for both ERA5 and bias-corrected CMIP6 data.
- Output: `OUT_DIR/cffdrs/era5/` and `OUT_DIR/cffdrs/{gcm}/`

### References

> Note that the set of functions performing the CFFDRS calculations are described in Van Wagner 1985 and 1987. Variable names were kept the same as in the original functions. Furthermore, in this set of function equation numbers also directly match up to those described in the original Fortran code. We also used `cffdrs` R package to help guide the creation and use of this script.

- Van Wagner, C.E. and T. L. Pickett. Equations and FORTRAN program for the Canadian Forest Fire Weather Index System. 1985. Canadian Forestry Service, Petawawa National Forestry Institute, Chalk River, Ontario. Forestry Technical Report 33. 18 p.

- Van Wagner, C.E. Development and structure of the Canadian Forest Fire Weather Index System. 1987. Canadian Forestry Service, Headquarters, Ottawa. Forestry Technical Report 35. 35 p.

- Wang, X. et al. “cffdrs: An R package for the Canadian Forest Fire Danger Rating System.” Ecological Processes, 6(1), 5. https://ecologicalprocesses.springeropen.com/articles/10.1186/s13717-017-0070-z