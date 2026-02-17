# Revision of bias corrected CMIP6 GCMs and Recalculation of CFFDRS Fire Weather Indices

Run the following pipeline to revise bias corrected + downscaled CMIP6 data, then recalculate CFFDRS fire weather indices from the results. 

This README assumes you have downloaded copies of the bias corrected + downscaled CMIP6 data from this USGS ScienceBase data release: https://www.sciencebase.gov/catalog/item/67ead89cd34ed02007f8357f. Contact SNAP for information about obtaining source data for this processing pipeline.

## Background

In the original processing of this dataset, some humidity values in the bias corrected + downscaled CMIP6 data (`hursmin` variable) were outside the range of 0-100. This was a data artifact of the quantile delta mapping method, and had downstream effects on the CFFDRS indices. The goal of this processing pipeline is to limit the `hursmin` variable values to 0-100 and recalculate the CFFDRS indices.

**NOTE:**  SNAP originally attempted to rerun the entire bias correction + downscaling pipeline from source data. Valid outputs were produced, but did not exactly match the original bias corrected + downscaled CMIP6 dataset. Since achieving byte-for-byte reproducibility of the original dataset was somewhat out of scope for this project, SNAP decided to focus on the revision of files in the data release. However, the entire pipeline is included here for future work; scripts `01_process_era.py`, `02_process_cmip6.py`, and `03_bias_correct_gcms.py` constitute this earlier part of the pipeline. The instructions below first describe the revision process and are limited to scripts `03b_fix_bias_corrected_gcms.py` and `04_calculate_cffdrs.py`. Instructions for the earlier part of the pipeline are included at the end of this README. 

### Setup (revision only)

Create a micromamba environment using the `processing/environment.yml`
```bash
micromamba env create -f processing/environment.yml
```

Set environment variables
```bash
# Path to bias corrected + downscaled CMIP6 files from the data release
export CMIP6_BIAS_CORRECTED=/path/to/existing/bias_corrected

# Output directory for revised bias corrected + downscaled CMIP6 files and recalculated CFFDRS indices
export OUT_DIR=/path/for/output/files
```

**Optional Settings:**

```bash
# Path to preprocessed ERA5 data (output of 01_process_era5.py) if wanting to recalculate CFFDRS indices for ERA5
# NOTE: these files are *not* included in the data release linked above!
export ERA5_PROCESSED=/path/to/era5/processed
```

**Optional Configuration Flags:**

```bash
# CLIP_HURSMIN: Clip relative humidity to valid range [0, 100]
# Set to FALSE if trying reproduce original dataset (default: TRUE)
export CLIP_HURSMIN=TRUE

# LEGACY_MODE: Replicate original pipeline behavior for testing
# Skips time alignment, dimension ordering, uses buggy chunking (default: FALSE)
# WARNING: For diagnostic purposes only, not for production use
export LEGACY_MODE=FALSE
```

### Revision Scripts

#### 03b_fix_bias_corrected_gcms.py
Apply `hursmin` clamping correction to existing bias-corrected + downscaled CMIP6 data:
- Reads existing bias-corrected + downscaled GCM files
- Applies `CLIP_HURSMIN` routine to clamp hursmin values to [0, 100]
- Writes corrected files to `OUT_DIR/bias_corrected/{gcm}/`

#### 04_calculate_cffdrs.py
Calculates Canadian Forest Fire Danger Rating System (CFFDRS) indices.
Calculated for bias-corrected + downscaled CMIP6 data and optionally for ERA5 data.
- Input from `OUT_DIR/bias_corrected/{gcm}/`
- Writes to `OUT_DIR/cffdrs/{gcm}/` and optionally to `OUT_DIR/cffdrs/era5/` 

### Directory Structure

The revision pipeline expects the following directory structure:
```
ERA5_PROCESSED/                   # Processed ERA5 (daily summaries)
├── tasmax_era5_1979.nc
├── pr_era5_1979.nc
├── sfcWind_era5_1979.nc
├── hursmin_era5_1979.nc
└── ...

CMIP6_BIAS_CORRECTED/                  # Bias corrected CMIP6 (from data release)
├── CNRM-CM6-1-HR/
│   ├── tasmax_CNRM-CM6-1-HR_1979.nc
│   └── ...
├── EC-Earth3-Veg/
├── MPI-ESM1-2-HR/
└── MRI-ESM2-0/
```

And produces the following directory structure:

```
OUT_DIR/                          # Final outputs
├── bias_corrected/               # Outputs from 03b_fix_bias_corrected_gcms.py
│   ├── CNRM-CM6-1-HR/
│   ├── EC-Earth3-Veg/
│   ├── MPI-ESM1-2-HR/
│   └── MRI-ESM2-0/
└── cffdrs/                       # Outputs from 04_calculate_cffdrs.py
    ├── era5/
    ├── CNRM-CM6-1-HR/
    ├── EC-Earth3-Veg/
    ├── MPI-ESM1-2-HR/
    └── MRI-ESM2-0/
```

### Running the Revision Pipeline

The pipeline takes between 3 and 4 hours to complete when run like so:

```bash
# Start a screen session on a compute node
srun --partition=t2small --cpus-per-task=24 --pty /bin/bash

# Navigate to the repo processing/ directory
cd ~/boreal_fire_weather/processing

# Activate environment
micromamba activate cffdrs

# Run the scripts
python 03b_fix_bias_corrected_gcms.py; python 04_calculate_cffdrs.py
```


### Quality Control

The QC scripts require the use of a separate environment, installed and activated like so:
```bash
micromamba env create -f qc/environment.yml
micromamba activate cffdrs_qc 
```

Here is a quick example comparing the revised files to the input files from the data release. This command chooses 10 random matching files and compares all variables, outputting a text summary to terminal and text file:
```bash
python qc/compare_datasets.py $CMIP6_BIAS_CORRECTED $OUT_DIR "" 10 --text-only | tee $OUT_DIR/qc/cmip6_qc.txt
python qc/compare_datasets.py /path/to/reference/cffdrs $OUT_DIR "" 10 --text-only | tee $OUT_DIR/qc/cffdrs_qc.txt
```

Alternatively, produce notebooks with a visual comparison of the files. This command chooses 5 random matching files for a single variable, and outputs a notebook comparison (with HTML version) for each file:
```bash
python qc/compare_datasets.py $CMIP6_BIAS_CORRECTED $OUT_DIR "hursmin" 5 --shapefile qc/shp/ecos.shp
python qc/compare_datasets.py /path/to/reference/cffdrs $OUT_DIR "ffmc" 5 --shapefile qc/shp/ecos.shp
```

We can also compare the NaN patterns between the revised files and the input files from the data release. A single notebook output (with HTML version) chooses a subset of NaN comparisons to plot, while the text output shows comparison results for all matching files. This command chooses 50 random matching files and compares NaN patterns for a single variable:
```bash
python qc/compare_nans.py $CMIP6_BIAS_CORRECTED $OUT_DIR $OUT_DIR/qc/nan_qc "ffmc" 50 --shapefile qc/shp/ecos.shp | tee $OUT_DIR/qc/cmip6_nan_qc.txt
python qc/compare_nans.py /path/to/reference/cffdrs $OUT_DIR $OUT_DIR/qc/nan_qc "ffmc" 50 --shapefile qc/shp/ecos.shp | tee $OUT_DIR/qc/cffdrs_nan_qc.txt
```

Alternatively, compare all files and all variables. 
```bash
python qc/compare_nans.py $CMIP6_BIAS_CORRECTED $OUT_DIR $OUT_DIR/qc/nan_qc "" --shapefile qc/shp/ecos.shp | tee $OUT_DIR/qc/cmip6_nan_qc.txt
python qc/compare_nans.py /path/to/reference/cffdrs $OUT_DIR $OUT_DIR/qc/nan_qc "" --shapefile qc/shp/ecos.shp | tee $OUT_DIR/qc/cffdrs_nan_qc.txt
```








 # ⚠️ SECTION BELOW is TBD, under construction! ⚠️
### Setup (complete bias correction pipeline)

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