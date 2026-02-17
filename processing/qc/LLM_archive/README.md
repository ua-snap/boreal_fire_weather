# Quality Control (QC) for Boreal Fire Weather Pipeline

## Overview

This directory contains QC tools and documentation for validating the bias correction pipeline against the original implementation ([Zenodo v7783759](https://zenodo.org/records/7783759)).

**Key Files:**
- `compare_datasets.py` - Script for comparing NetCDF outputs between pipeline versions
- `compare_nans.py` - Script for focused NaN pattern analysis with smart grouping
- `qc_utils.py` - Utility functions for dataset comparison and visualization
- `visualize_comparison.ipynb` - Parameterized notebook template for visual QC reports
- `visualize_nan_comparison.ipynb` - Parameterized notebook template for NaN-focused reports
- `pipeline_comparison_analysis.md` - Detailed analysis of differences between pipelines (summarized below)
- `qc/bias_corrected/original_vs_fully_processed.txt` - Results from initial QC comparison

---

## Background: Why QC Failed

The refactored pipeline produces **different results** from the original, with all 80 test files showing discrepancies:

### Two Types of Differences Found

1. **Shape/Dimension Ordering** (CNRM-CM6-1-HR only)
   - Original: `(lat, lon, time)` ordering
   - New: `(time, lat, lon)` ordering
   
2. **Numerical Differences** (All other models)
   - ~45-50% of grid cells differ
   - Mean differences: 1e-3 to 4e-2 (1-4% of typical values)
   - Max differences: 0.5 to 23 units depending on variable

### Root Cause Analysis

The primary cause is **time coordinate handling**. The refactored pipeline adds:

1. **Time alignment** - Ensures ERA5 and CMIP6 timestamps match exactly
   ```python
   if not ref.time.equals(hst.time):
       hst = hst.reindex(time=ref.time, method="nearest", tolerance=pd.Timedelta("12h"))
   ref = ref.assign_coords(time=ref.time.dt.floor("D"))  # Floor to midnight
   ```

2. **Dimension standardization** - Enforces consistent `(time, lat, lon)` ordering
   ```python
   export_ds = export_ds.transpose("time", "lat", "lon")
   ```

The original pipeline lacked both, causing:
- Temporal misalignment in QDM monthly groupings
- Inconsistent output dimension ordering

**See [pipeline_comparison_analysis.md](pipeline_comparison_analysis.md) for complete technical details.**

---

## Legacy Mode: Testing the Hypothesis

To validate that time alignment causes the differences, a `LEGACY_MODE` flag was added to skip these improvements:

### Configuration Variables

Both defined in `config.py`:

```python
# Skip time alignment and transpose (reproduces original behavior)
LEGACY_MODE = os.getenv("LEGACY_MODE", "FALSE").upper() == "TRUE"

# Clip relative humidity to [0, 100] range
CLIP_HURSMIN = os.getenv("CLIP_HURSMIN", "TRUE").upper() == "TRUE"
```

### When to Use Each Mode

| Mode | Use For | Don't Use For |
|------|---------|---------------|
| **Normal** (default) | Production runs, final datasets | Backward compatibility testing |
| **Legacy** (LEGACY_MODE=TRUE, CLIP_HURSMIN=FALSE) | Testing, validation, hypothesis confirmation | Production, when robustness matters |

**⚠️ Legacy mode reproduces original pipeline issues:**
- Temporal misalignment between ERA5/CMIP6
- Inconsistent dimension ordering (especially CNRM)
- Helper function bugs (masked by exception handling)

**✅ Normal mode improvements:**
- Proper temporal alignment
- Consistent output format
- More reproducible results
- Bug fixes

---

## QC Script: `compare_datasets.py`

### Features

- **Visual HTML reports by default** for failed CFFDRS comparisons
- Compares all data variables in matching NetCDF files
- **Automatically handles dimension ordering differences**
- Optional text-only mode for traditional QC reports
- Optional filtering by variable names
- Optional random sampling to limit files checked per variable
- Detailed difference reporting (max, mean, count)
- Clear pass/fail status

### Usage

```bash
python compare_datasets.py <old_dir> <new_dir> <output_dir> [variables] [max_files] [--text-only] [-o output.txt]
```

**Arguments:**
- `old_dir` - Path to reference/old dataset directory
- `new_dir` - Path to new dataset directory to validate
- `output_dir` - Directory where QC outputs (notebooks/html) will be saved
- `variables` - (Optional) Comma-separated variable names (e.g., `tasmax,pr,hursmin`)
- `max_files` - (Optional) Maximum files to randomly sample per variable
- `--text-only` - (Optional) Generate text report only, skip visualizations
- `-o, --output` - (Optional) Save text report to file instead of stdout
- `--shapefile` - (Optional) Custom shapefile path for visualizations

**Examples:**

```bash
# Generate visual HTML reports (default)
python compare_datasets.py /data/old /data/new /output/qc

# Text report to stdout
python compare_datasets.py /data/old /data/new /output/qc --text-only

# Text report to file
python compare_datasets.py /data/old /data/new /output/qc --text-only -o qc_report.txt

# Compare specific variables with visualizations
python compare_datasets.py /data/old /data/new /output/qc tasmax,pr,hursmin

# Sample 10 files per variable with text output
python compare_datasets.py /data/old /data/new /output/qc tasmax,pr 10 --text-only

# Sample 20 files across all variables
python compare_datasets.py /data/old /data/new /output/qc "" 20
```

### Output

**Default (visualizations):**
- HTML reports in `<output_dir>/html/`
- Executed notebooks in `<output_dir>/notebooks/`
- Console summary of comparison results

**With `--text-only`:**
- Per-file pass/fail status
- Difference statistics for failed files
- Summary counts by variable
- Overall QC result (exit 0 = pass, exit 1 = fail)
- Output to stdout or file (with `-o` flag)

---

## Testing Workflow

Complete workflow to validate the time alignment hypothesis:

### 1. Run Legacy Mode (Reproduces Original)

```bash
export LEGACY_MODE=TRUE
export CLIP_HURSMIN=FALSE  # Disable clipping for pure comparison
export OUT_DIR=/data/test_legacy
python 03_bias_correct_gcms.py
```

### 2. Run Normal Mode (Refactored Pipeline)

```bash
export LEGACY_MODE=FALSE  # or omit (default)
export CLIP_HURSMIN=FALSE
export OUT_DIR=/data/test_normal
python 03_bias_correct_gcms.py
```

### 3. Compare Legacy vs Original

```bash
python qc/compare_datasets.py \
    /data/original_zenodo \
    /data/test_legacy \
    "tasmax,pr,hursmin,sfcWind" \
    20
```

**Expected result:** Close match (or identical) to original dataset

### 4. Compare Normal vs Original

```bash
python qc/compare_datasets.py \
    /data/original_zenodo \
    /data/test_normal \
    "tasmax,pr,hursmin,sfcWind" \
    20
```

**Expected result:** Differences in values due to improved time alignment

### 5. Compare Legacy vs Normal

```bash
python qc/compare_datasets.py \
    /data/test_legacy \
    /data/test_normal \
    "tasmax,pr,hursmin,sfcWind" \
    20
```

**Expected result:** Shows impact of time alignment and transpose improvements

---

## Interpretation Guide

### If Legacy Mode Matches Original

✅ **Confirms hypothesis:** Time alignment is the primary cause of differences

**Recommendation:** Adopt normal mode (with improvements) for production

### If Legacy Mode Still Differs

⚠️ **Additional factors involved:** Investigate further

**Check:**
- Package versions (especially `xclim`, `xarray`, `dask`)
- Input data preprocessing
- Random seed effects (if any)
- Floating-point precision issues

### Acceptable Differences in Normal Mode

**Small numerical differences** (~1-4% mean change) are expected and acceptable because:
- Better temporal alignment changes which data points are paired
- Monthly quantile groupings are more accurate
- Results are more reproducible and scientifically sound

**These are improvements, not errors.**

---

## Quick Reference Commands

```bash
# Set both flags for complete legacy behavior
export LEGACY_MODE=TRUE CLIP_HURSMIN=FALSE

# Run bias correction
python 03_bias_correct_gcms.py

# Visual QC check (default - generates HTML reports)
python qc/compare_datasets.py /old /new /output/qc

# Quick text QC check (sample 20 files)
python qc/compare_datasets.py /old /new /output/qc "" 20 --text-only

# Full text QC on all files to output file
python qc/compare_datasets.py /old /new /output/qc --text-only -o qc_report.txt

# QC specific variables with visualizations
python qc/compare_datasets.py /old /new /output/qc tasmax,pr,hursmin
```

---

---

## Visual QC System

The QC script **generates visual HTML reports by default** that compare datasets across all timesteps.

### Setup

1. **Create the QC conda environment:**

```bash
cd processing/qc
conda env create -f environment.yml
conda activate cffdrs_qc
```

2. **The environment includes:**
   - `papermill` - For parameterized notebook execution
   - `nbconvert` - For HTML export
   - `matplotlib`, `geopandas` - For mapping
   - `pandas` - For summary tables
   - All necessary data science packages

### How It Works

When `compare_datasets.py` detects differences, it automatically:

1. **Compares** all timesteps individually between reference and new datasets
2. **Identifies** timesteps where differences occur (value changes or NaN pattern changes)
3. **Generates** summary statistics table showing only timesteps with differences
4. **Creates** 4-panel comparison plots for:
   - Timestep with maximum value difference
   - Timestep with maximum NaN pattern changes (if different from above)
5. **Each 4-panel plot shows:**
   - Panel 1: Original dataset values
   - Panel 2: New dataset values
   - Panel 3: Difference map (new - original) with diverging colormap
   - Panel 4: NaN pattern comparison
6. **Overlays** boreal region shapefile boundary on all maps
7. **Produces** a markdown-formatted QC summary with key findings
8. **Exports** executed notebooks as HTML reports

### Usage

**Generate visualizations (default behavior):**
```bash
python compare_datasets.py /old/cffdrs /new/cffdrs /output/qc
```

**Generate text report only:**
```bash
python compare_datasets.py /old/cffdrs /new/cffdrs /output/qc --text-only
```

**Save text report to file:**
```bash
python compare_datasets.py /old/cffdrs /new/cffdrs /output/qc --text-only -o qc_report.txt
```

**Specify custom shapefile:**
```bash
python compare_datasets.py /old/cffdrs /new/cffdrs /output/qc --shapefile ../shp/custom_boundary.shp
```

### Output Structure

```
<output_dir>/
├── notebooks/          # Executed Jupyter notebooks (.ipynb)
│   ├── cffdrs_CNRM-CM6-1-HR_1980_comparison.ipynb
│   └── ...
└── html/              # HTML reports (open in browser)
    ├── cffdrs_CNRM-CM6-1-HR_1980_comparison.html
    └── ...
```

### Visualization Features

#### Timestep-by-Timestep Analysis
- **Comprehensive comparison:** All timesteps analyzed automatically
- **Smart filtering:** Summary statistics shown only for timesteps with differences
- **Targeted visualization:** Plots focus on timesteps with greatest differences

#### Summary Statistics Table
For each timestep with differences, displays:
- Date
- Mean and standard deviation (reference vs new)
- Maximum absolute difference
- Count of cells with value differences
- Count of cells with NaN pattern changes

#### Four-Panel Plots
- **Original data** (warm colormap)
- **New data** (warm colormap, same scale as original)
- **Difference** (diverging red-blue colormap, white = no change)
- **NaN patterns** (white = both valid, colors indicate NaN mismatches)

#### Spatial Context
- Boreal region boundary overlaid in black
- Map extent set to shapefile bounds
- Lat/lon axes for georeferencing

#### QC Summary
Markdown-formatted summary including:
- Total timesteps analyzed
- Count of timesteps with differences
- Maximum value difference and when it occurred
- Maximum NaN changes and when they occurred
- Overall statistics across all timesteps

### QC Utilities Module

The `qc_utils.py` module provides reusable functions:

- `compare_timestep()` - Compare single timestep between datasets
- `compare_all_timesteps()` - Compare all timesteps, return DataFrame
- `create_comparison_plot()` - Generate 4-panel visualization
- `create_summary_table()` - Format statistics table for display
- `generate_qc_summary()` - Create markdown-formatted findings

These functions can be imported and used in custom analysis notebooks or scripts.

### Customization

Edit `visualize_comparison.ipynb` to customize:
- Colormap choices
- Percentile limits for color scaling
- Figure sizes
- Additional plots or statistics
- Table formatting

The template uses papermill's `parameters` tag, so key settings can be injected at runtime.

---

## NaN Pattern Analysis: `compare_nans.py`

For datasets where **NaN pattern differences are the primary concern**, use the specialized NaN comparison workflow.

### Purpose

While `compare_datasets.py` compares all aspects of datasets (values + NaNs), `compare_nans.py` focuses exclusively on **NaN pattern differences** with intelligent grouping:

- **Groups files** by model + variable (e.g., "CNRM-CM6-1-HR_tasmax")
- **Selects representatives** from each group (file with max NaN differences)
- **Limits output** to 10 plots maximum for manageable reports
- **Analyzes all timesteps** to find worst-case NaN mismatches
- **Color-coded visualization** showing where NaNs appeared/disappeared

### Usage

```bash
python compare_nans.py <old_dir> <new_dir> <output_dir> [variables] [max_files] [--shapefile <path>]
```

**Arguments:**
- `old_dir` - Path to reference/old dataset directory
- `new_dir` - Path to new dataset directory
- `output_dir` - Directory for NaN analysis output
- `variables` - (Optional) Comma-separated variable names to analyze
- `max_files` - (Optional) Maximum files to sample per variable before grouping
- `--shapefile` - (Optional) Custom shapefile path for boundary overlay

**Examples:**

```bash
# Analyze all CFFDRS files
python compare_nans.py /data/old/cffdrs /data/new/cffdrs /output/nan_analysis

# Analyze specific variables
python compare_nans.py /data/old /data/new /output/nan_analysis tasmax,pr,hursmin

# Sample 50 files, then group and select representatives
python compare_nans.py /data/old /data/new /output/nan_analysis "" 50

# With custom shapefile
python compare_nans.py /data/old/cffdrs /data/new/cffdrs /output/nan_analysis --shapefile shp/ecos.shp
```

### How It Works

1. **Analysis Phase:**
   - Scans all matching file pairs
   - For each file, compares NaN patterns across **all timesteps**
   - Calculates total NaN differences and worst timestep

2. **Grouping Phase:**
   - Groups files by `model_variable` (e.g., "MPI-ESM1-2-HR_tasmax")
   - Files in same group have similar NaN patterns (same model/variable combination)

3. **Selection Phase:**
   - Selects **one representative file per group** (max NaN differences)
   - Limits to **10 groups maximum** for concise reports

4. **Visualization Phase:**
   - Generates HTML report with one plot per representative file
   - Each plot shows the **timestep with maximum NaN pattern difference**
   - Color codes: Red (NaN→Valid), Blue (Valid→NaN), Gray (both NaN)

### Output Structure

```
<output_dir>/
├── nan_pattern_comparison.ipynb    # Executed notebook
└── nan_pattern_comparison.html     # HTML report with plots
```

### Visualization Features

#### Representative Sampling
- **Intelligent grouping** prevents redundant plots of similar files
- **Max 10 plots** keeps reports focused and manageable
- **Group statistics** show how many files each representative covers

#### Summary Tables
- **Overview statistics:** Total groups, files shown, total NaN differences
- **Group-level table:** Model, variable, group size, representative file chosen
- **Per-file details:** Total NaN diffs, timesteps with diffs, worst timestep count

#### NaN Pattern Maps
Each plot shows:
- **White areas:** Valid data in both datasets
- **Red areas:** Original had NaN, new has valid data (NaN→Valid)
- **Blue areas:** Original had valid data, new has NaN (Valid→NaN)
- **Gray areas:** Both datasets have NaN
- **Black boundary:** Boreal region shapefile overlay

#### When to Use This vs `compare_datasets.py`

| Use `compare_datasets.py` | Use `compare_nans.py` |
|--------------------------|----------------------|
| Need value difference analysis | Only care about NaN patterns |
| Comparing small number of files | Analyzing hundreds of files |
| Want detailed per-file reports | Want high-level summary by model+variable |
| All files are important | Many files have similar patterns |

---

## Further Documentation

- **[pipeline_comparison_analysis.md](pipeline_comparison_analysis.md)** - Comprehensive technical analysis of all differences
- **[cffdrs_calculation_analysis.md](cffdrs_calculation_analysis.md)** - Math investigation of NaN propagation in CFFDRS calculations
- **qc/bias_corrected/original_vs_fully_processed.txt** - Raw QC results showing specific failures
- **Original dataset:** [Zenodo v7783759](https://zenodo.org/records/7783759)

