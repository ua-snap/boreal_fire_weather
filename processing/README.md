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

**Optional Settings:**

Shapefile to subset the data during bias correction:
```bash
export SHP_MASK=/path/to/shapefile/mask.shp
```

External bias-corrected data (for 03b script input only):
```bash
# Points 03b to existing bias-corrected files to apply hursmin correction
# If not set, 03b will use OUT_DIR/bias_corrected as both input and output
export CMIP6_BIAS_CORRECTED=/path/to/existing/bias_corrected
```

**Optional Configuration Flags:**

```bash
# CLIP_HURSMIN: Clip relative humidity to valid range [0, 100]
# Set to FALSE for comparison with original dataset (default: TRUE)
export CLIP_HURSMIN=TRUE

# LEGACY_MODE: Replicate original pipeline behavior for testing
# Skips time alignment, dimension ordering, uses buggy chunking (default: FALSE)
# WARNING: For diagnostic purposes only, not for production use
export LEGACY_MODE=FALSE
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
- **New features:**
  - Explicit time coordinate alignment to prevent temporal mismatches
  - Consistent dimension ordering: `(time, lat, lon)` across all outputs
  - Optimized Dask chunking for improved performance
  - Optional clipping of `hursmin` to physically valid range [0, 100]
  - Distributed computing with `persist()` for efficient data reuse
  - Configuration via `CLIP_HURSMIN` and `LEGACY_MODE` environment variables
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
- Input: `OUT_DIR/bias_corrected/{gcm}/` (output from 03 or 03b)
- Output: `OUT_DIR/cffdrs/era5/` and `OUT_DIR/cffdrs/{gcm}/`

#### 03b_fix_bias_corrected_gcms.py
Apply hursmin clamping correction to existing bias-corrected data:
- Reads existing bias-corrected GCM files
- Applies `CLIP_HURSMIN` routine to clamp hursmin values to [0, 100]
- Writes corrected files to `OUT_DIR/bias_corrected/{gcm}/`
- Use when working with external bias-corrected data that wasn't processed with the clipping correction
- Input: `CMIP6_BIAS_CORRECTED/` (or `OUT_DIR/bias_corrected/` if not set)
- Output: `OUT_DIR/bias_corrected/{gcm}/`

### Alternative Workflows

#### Option 1: Full Pipeline (Default)
Process everything from scratch:
```bash
python 01_process_era5.py      # Process ERA5 raw data
python 02_process_cmip6.py     # Process CMIP6 raw data
python 03_bias_correct_gcms.py # Apply bias correction with QDM
python 04_calculate_cffdrs.py  # Calculate fire weather indices
```

#### Option 2: Using Pre-processed Data
Skip steps 01-02 if you already have processed daily data:
```bash
# Point to existing processed data
export ERA5_PROCESSED=/path/to/existing/era5/processed
export CMIP6_PROCESSED=/path/to/existing/cmip6/processed

# Run bias correction and CFFDRS calculation
python 03_bias_correct_gcms.py
python 04_calculate_cffdrs.py
```

#### Option 3: Using External Bias-Corrected Data (Correction Only)
If you have existing bias-corrected data (e.g., from an earlier version of the pipeline) that needs the `hursmin` clamping correction applied:

```bash
# Point to external bias-corrected data
export CMIP6_BIAS_CORRECTED=/path/to/existing/bias_corrected

# Point to ERA5 processed data (needed for CFFDRS calculation)
export ERA5_PROCESSED=/path/to/existing/era5/processed

# Set output directory for corrected files
export OUT_DIR=/path/for/corrected/output

# Enable hursmin clipping (if not already applied)
export CLIP_HURSMIN=TRUE

# Apply hursmin clamping correction only
python 03b_fix_bias_corrected_gcms.py

# Calculate CFFDRS indices using the corrected data from OUT_DIR
python 04_calculate_cffdrs.py
```

**Note:** The `03b_fix_bias_corrected_gcms.py` script reads from `CMIP6_BIAS_CORRECTED` (if set) or `OUT_DIR/bias_corrected` (if not set), applies corrections, and writes all variables to `OUT_DIR/bias_corrected`. The `04_calculate_cffdrs.py` script always reads from `OUT_DIR/bias_corrected`.

### Quality Control

#### Comparing Outputs

This refactored pipeline includes methodological improvements over the original ([Zenodo v7783759](https://zenodo.org/records/7783759)). Due to these improvements, outputs will differ numerically from the original dataset. See `qc/pipeline_comparison_analysis.md` for detailed technical analysis.

To compare outputs with a reference dataset:

```bash
# Run pipeline (set CLIP_HURSMIN=FALSE for comparison with original dataset)
export CLIP_HURSMIN=FALSE
python 03_bias_correct_gcms.py
python 04_calculate_cffdrs.py

# Compare bias corrected outputs
cd qc/bias_corrected/
python ../compare_datasets.py /path/to/reference/bias_corrected /path/to/new/bias_corrected

# Compare CFFDRS outputs
cd ../cffdrs/
python ../compare_datasets.py /path/to/reference/cffdrs /path/to/new/cffdrs
```

See `qc/README.md` for detailed QC workflow documentation.

#### LEGACY_MODE Testing

For diagnostic purposes, a `LEGACY_MODE` flag attempts to replicate the original pipeline's behavior:

```bash
export LEGACY_MODE=TRUE
export CLIP_HURSMIN=FALSE
python 03_bias_correct_gcms.py
```

**Warning:** LEGACY_MODE intentionally includes bugs from the original code and should **not** be used for production work.

### Additional Documentation

- **`qc/pipeline_comparison_analysis.md`** - Technical analysis of differences between this pipeline and the original version
- **`qc/README.md`** - QC workflow documentation

### References

> Note that the set of functions performing the CFFDRS calculations are described in Van Wagner 1985 and 1987. Variable names were kept the same as in the original functions. Furthermore, in this set of function equation numbers also directly match up to those described in the original Fortran code. We also used `cffdrs` R package to help guide the creation and use of this script.

- Van Wagner, C.E. and T. L. Pickett. Equations and FORTRAN program for the Canadian Forest Fire Weather Index System. 1985. Canadian Forestry Service, Petawawa National Forestry Institute, Chalk River, Ontario. Forestry Technical Report 33. 18 p.

- Van Wagner, C.E. Development and structure of the Canadian Forest Fire Weather Index System. 1987. Canadian Forestry Service, Headquarters, Ottawa. Forestry Technical Report 35. 35 p.

- Wang, X. et al. “cffdrs: An R package for the Canadian Forest Fire Danger Rating System.” Ecological Processes, 6(1), 5. https://ecologicalprocesses.springeropen.com/articles/10.1186/s13717-017-0070-z