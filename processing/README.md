# Revision of bias corrected CMIP6 GCMs and Recalculation of CFFDRS Fire Weather Indices

Run the following pipeline to revise bias corrected + downscaled CMIP6 data, then recalculate CFFDRS fire weather indices from the results. 

This README assumes you have downloaded copies of the bias corrected + downscaled CMIP6 data from this USGS ScienceBase data release: https://www.sciencebase.gov/catalog/item/67ead89cd34ed02007f8357f. Contact SNAP for information about obtaining additional source data for this processing pipeline.

## Background

In the original processing of this dataset, some humidity values in the bias corrected + downscaled CMIP6 data (`hursmin` variable) were outside the range of 0-100. This was a data artifact specific to the quantile delta mapping method used, and had downstream effects on the CFFDRS indices. The goal of this processing pipeline is to limit the `hursmin` variable values to a range of 0-100 and recalculate the CFFDRS indices.

**NOTE:**  SNAP originally attempted to rerun the entire bias correction + downscaling pipeline from source data. Valid outputs were produced, but did not exactly match the original bias corrected + downscaled CMIP6 dataset. Since achieving byte-for-byte reproducibility of the original dataset was somewhat out of scope for this project, SNAP decided to focus on the revision of bias corrected + downscaled CMIP6 files and corresponding CFFDRS files in the data release. However, the entire pipeline is included here for future work; scripts `01_process_era.py`, `02_process_cmip6.py`, and `03_bias_correct_gcms.py` constitute this earlier part of the pipeline, and instructions for running the full pipeline are included at the end of this README. The instructions in the "Revision-Only Pipeline" section below describe the revision process and are limited to scripts `03b_fix_bias_corrected_gcms.py` and `04_calculate_cffdrs.py`.


## Revision-Only Pipeline

### 1. Setup

Create a micromamba environment using the `processing/environment.yml`
```bash
micromamba env create -f processing/environment.yml
```

Set environment variables:
```bash
# Path to bias corrected + downscaled CMIP6 files from the data release
export CMIP6_BIAS_CORRECTED=/path/to/existing/data

# Output directory for revised bias corrected + downscaled CMIP6 files and recalculated CFFDRS indices
export OUT_DIR=/path/for/output/files
```

**Optional Settings:**

```bash
# Path to preprocessed ERA5 data (output of 01_process_era5.py) if wanting to recalculate CFFDRS indices for ERA5
# NOTE: these files are *not* included in the data release linked above!
export ERA5_PROCESSED=/path/to/era5/processed
```

### 2. Revision Scripts

#### `03b_fix_bias_corrected_gcms.py`
- Applies `CLIP_HURSMIN` routine to clamp hursmin values to [0, 100]
- Writes corrected files to `OUT_DIR/bias_corrected/{gcm}/`

#### `04_calculate_cffdrs.py`
- Calculates Canadian Forest Fire Danger Rating System (CFFDRS) indices for bias-corrected + downscaled CMIP6 data, and optionally for ERA5 data.
- Input from `OUT_DIR/bias_corrected/{gcm}/`
- Writes to `OUT_DIR/cffdrs/{gcm}/` and optionally to `OUT_DIR/cffdrs/era5/` 

### 3. Directory Structure

The revision pipeline expects the following directory structure:
```
ERA5_PROCESSED/                   # Optional processed ERA5 (daily summaries)
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

### 4. Running the Revision Pipeline

The pipeline takes between 3 and 4 hours to complete when run like so:

```bash
# Start a screen session on a compute node
srun --partition=t2small --cpus-per-task=24 --pty /bin/bash

# Navigate to the repo processing/ directory
cd ~/boreal_fire_weather/processing

# Activate environment
micromamba activate cffdrs

# Run the scripts
python 03b_fix_bias_corrected_gcms.py
python 04_calculate_cffdrs.py
```


### 5. Quality Control

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

Alternatively, produce notebooks with visual comparison of files that fail QC. This command chooses 5 random matching files for a single variable, and outputs a notebook comparison (with HTML version) for each file:
```bash
python qc/compare_datasets.py $CMIP6_BIAS_CORRECTED $OUT_DIR "hursmin" 5 

python qc/compare_datasets.py /path/to/reference/cffdrs $OUT_DIR "ffmc" 5 
```

We can also compare the NaN patterns between the revised files and the input files from the data release. A single notebook output (with HTML version) chooses a subset of failing NaN comparisons to plot, while the text output shows comparison results for all matching files. This command chooses 50 random matching files and compares NaN patterns for a single variable:
```bash
python qc/compare_nans.py $CMIP6_BIAS_CORRECTED $OUT_DIR $OUT_DIR/qc/nan_qc "ffmc" 50 | tee $OUT_DIR/qc/cmip6_nan_qc.txt

python qc/compare_nans.py /path/to/reference/cffdrs $OUT_DIR $OUT_DIR/qc/nan_qc "ffmc" 50 | tee $OUT_DIR/qc/cffdrs_nan_qc.txt
```

Alternatively, compare all files and all variables. 
```bash
python qc/compare_nans.py $CMIP6_BIAS_CORRECTED $OUT_DIR $OUT_DIR/qc/nan_qc "" | tee $OUT_DIR/qc/cmip6_nan_qc.txt

python qc/compare_nans.py /path/to/reference/cffdrs $OUT_DIR $OUT_DIR/qc/nan_qc "" | tee $OUT_DIR/qc/cffdrs_nan_qc.txt
```







## Complete Bias Correction + Revision Pipeline

### 1. Setup

⚠️ _Note that the code below was unable to identically reproduce the original dataset; there may be unidentified issues with the codebase, so please use with caution!_

Create a micromamba environment using the `processing/environment.yml`
```bash
micromamba env create -f processing/environment.yml
micromamba activate cffdrs
```

Set environment variables
```bash
# Input directories (raw source data -  not included in data release linked above!)
export ERA5_IN=/path/to/era5/input/data
export CMIP6_IN=/path/to/cmip6/input/data

# Processed data directories (output of 01_process_era5.py and 02_process_cmip6.py scripts or obtained from author -  not included in data release linked above!)
export ERA5_PROCESSED=/path/to/era5/processed
export CMIP6_PROCESSED=/path/to/cmip6/processed

# Output directory for bias corrected GCMs and CFFDRS indices
export OUT_DIR=/path/for/output/files
```

**Optional Settings:**

```bash
# Shapefile to subset the data during bias correction.
export SHP_MASK=/path/to/shapefile/mask.shp

# Enable/disable hursmin clipping to 0-100 range; default is TRUE.
export CLIP_HURSMIN=TRUE

# Legacy mode skips time alignment and dimension transposition to match old pipeline behavior; default is FALSE.
# Warning: LEGACY_MODE intentionally includes bugs from the original code and should **not** be used for production work!
export LEGACY_MODE=FALSE

# Path to existing bias-corrected CMIP6 data if skipping bias-correction step.
# If not set, default is to use OUT_DIR/bias_corrected 
export CMIP6_BIAS_CORRECTED=/path/to/existing/data
```


### 2. Pipeline Scripts

#### `01_process_era5.py`
Process ERA5 hourly data to daily summaries:
- Calculate daily maximum temperature (`tasmax`) from hourly temperature
- Calculate daily total precipitation (`pr`) from hourly precipitation
- Calculate daily mean wind speed (`sfcWind`) from `u10`/`v10` components
- Calculate daily minimum relative humidity (`hursmin`) from temperature and dewpoint
- Modify geographic extent to region of interest
- Standardize calendar to 'noleap' (365-day years)
- Input: `ERA5_IN/`
- Output: `ERA5_PROCESSED/`

#### `02_process_cmip6.py`
Process CMIP6 daily data:
- Modify geographic extent to match ERA5
- Standardize calendar to 'noleap'
- Convert units (e.g., K → degC, m/s → km/hour)
- Calculate wind speed from u/v components where needed
- Input: `CMIP6_IN/{gcm}/`
- Output: `CMIP6_PROCESSED/{gcm}/`

#### `03_bias_correct_gcms.py`
Apply Quantile Delta Mapping (QDM) bias correction to CMIP6 data using ERA5 as reference:
- Train QDM model using historical overlap period (1980-2009)
- Apply bias correction to historical GCM data (1980-2009)
- Apply bias correction to future projection periods (2010-2039, 2040-2069, 2070-2099)
- Handle precipitation thresholds with jitter method
- **⭐️ New features in this codebase:**
  - Optional spatial mask from shapefile
  - Explicit time coordinate alignment to prevent temporal mismatches
  - Consistent dimension ordering: `(time, lat, lon)` across all outputs
  - Optimized Dask chunking for improved performance
  - Optional clipping of `hursmin` to physically valid range [0, 100]
  - Distributed computing with `persist()` for efficient data reuse
  - Configuration via `CLIP_HURSMIN` and `LEGACY_MODE` environment variables
- ⚠️ New features can be disabled by setting `CLIP_HURSMIN=FALSE` and `LEGACY_MODE=TRUE` 
- Output: `OUT_DIR/bias_corrected/{gcm}/`

#### `04_calculate_cffdrs.py`
- Calculate Canadian Forest Fire Danger Rating System (CFFDRS) indices for both ERA5 and bias-corrected CMIP6 data.
- Input: `OUT_DIR/bias_corrected/{gcm}/`
- Output: `OUT_DIR/cffdrs/era5/` and `OUT_DIR/cffdrs/{gcm}/`

### 3. Directory Structure

The full pipeline expects the following directory structure:
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

CMIP6_PROCESSED/                  # Processed CMIP6
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
├── bias_corrected/               # Outputs from 03_fix_bias_corrected_gcms.py
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


### 4. Running the Full Pipeline

```bash
# Start a screen session on a compute node
srun --partition=t2small --cpus-per-task=24 --pty /bin/bash

# Navigate to the repo processing/ directory
cd ~/boreal_fire_weather/processing

# Activate environment
micromamba activate cffdrs

# Run the scripts
python 01_process_era5.py
python 02_process_cmip6.py
python 03_fix_bias_corrected_gcms.py
python 04_calculate_cffdrs.py
```

### 5. Quality Control

The QC scripts require the use of a separate environment, installed and activated like so:
```bash
micromamba env create -f qc/environment.yml
micromamba activate cffdrs_qc 
```

Here is a quick example comparing the revised files to the input files from the data release. This command chooses 10 random matching files and compares all variables, outputting a text summary to terminal and text file:
```bash
python qc/compare_datasets.py /path/to/reference/cmip6 $OUT_DIR "" 10 --text-only | tee $OUT_DIR/qc/cmip6_qc.txt

python qc/compare_datasets.py /path/to/reference/cffdrs $OUT_DIR "" 10 --text-only | tee $OUT_DIR/qc/cffdrs_qc.txt
```

Alternatively, produce notebooks with visual comparison of files that fail QC. This command chooses 5 random matching files for a single variable, and outputs a notebook comparison (with HTML version) for each file:
```bash
python qc/compare_datasets.py /path/to/reference/cmip6 $OUT_DIR "hursmin" 5 

python qc/compare_datasets.py /path/to/reference/cffdrs $OUT_DIR "ffmc" 5 
```

We can also compare the NaN patterns between the revised files and the input files from the data release. A single notebook output (with HTML version) chooses a subset of failing NaN comparisons to plot, while the text output shows comparison results for all matching files. This command chooses 50 random matching files and compares NaN patterns for a single variable:
```bash
python qc/compare_nans.py /path/to/reference/cmip6 $OUT_DIR $OUT_DIR/qc/nan_qc "ffmc" 50 | tee $OUT_DIR/qc/cmip6_nan_qc.txt

python qc/compare_nans.py /path/to/reference/cffdrs $OUT_DIR $OUT_DIR/qc/nan_qc "ffmc" 50 | tee $OUT_DIR/qc/cffdrs_nan_qc.txt
```

Alternatively, compare all files and all variables. 
```bash
python qc/compare_nans.py /path/to/reference/cmip6 $OUT_DIR $OUT_DIR/qc/nan_qc "" | tee $OUT_DIR/qc/cmip6_nan_qc.txt

python qc/compare_nans.py /path/to/reference/cffdrs $OUT_DIR $OUT_DIR/qc/nan_qc "" | tee $OUT_DIR/qc/cffdrs_nan_qc.txt
```



## References

> Note that the set of functions performing the CFFDRS calculations are described in Van Wagner 1985 and 1987. Variable names were kept the same as in the original functions. Furthermore, in this set of function equation numbers also directly match up to those described in the original Fortran code. The original author also used `cffdrs` R package to help guide the creation and use of the `04_calculate_cffdrs.py` script.

- Van Wagner, C.E. and T. L. Pickett. Equations and FORTRAN program for the Canadian Forest Fire Weather Index System. 1985. Canadian Forestry Service, Petawawa National Forestry Institute, Chalk River, Ontario. Forestry Technical Report 33. 18 p.

- Van Wagner, C.E. Development and structure of the Canadian Forest Fire Weather Index System. 1987. Canadian Forestry Service, Headquarters, Ottawa. Forestry Technical Report 35. 35 p.

- Wang, X. et al. “cffdrs: An R package for the Canadian Forest Fire Danger Rating System.” Ecological Processes, 6(1), 5. https://ecologicalprocesses.springeropen.com/articles/10.1186/s13717-017-0070-z