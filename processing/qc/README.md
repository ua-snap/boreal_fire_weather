# Quality Control (QC) for Boreal Fire Weather Pipeline

## Overview

This directory contains QC tools and documentation for validating the bias correction pipeline against the original implementation ([Zenodo v7783759](https://zenodo.org/records/7783759)).

**Key Files:**
- `compare_datasets.py` - Script for comparing NetCDF outputs between pipeline versions
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

- Compares all data variables in matching NetCDF files
- **Automatically handles dimension ordering differences**
- Optional filtering by variable names
- Optional random sampling to limit files checked per variable
- Detailed difference reporting (max, mean, count)
- Clear pass/fail status

### Usage

```bash
python compare_datasets.py <old_dir> <new_dir> [variables] [max_files]
```

**Arguments:**
- `old_dir` - Path to reference/old dataset directory
- `new_dir` - Path to new dataset directory to validate
- `variables` - (Optional) Comma-separated variable names (e.g., `tasmax,pr,hursmin`)
- `max_files` - (Optional) Maximum files to randomly sample per variable

**Examples:**

```bash
# Compare all files
python compare_datasets.py /data/old /data/new

# Compare specific variables
python compare_datasets.py /data/old /data/new tasmax,pr,hursmin

# Sample 10 files per variable
python compare_datasets.py /data/old /data/new tasmax,pr 10

# Sample 20 files across all variables
python compare_datasets.py /data/old /data/new "" 20
```

### Output

- Per-file pass/fail status
- Difference statistics for failed files
- Summary counts by variable
- Overall QC result (exit 0 = pass, exit 1 = fail)

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

# Quick QC check (sample 20 files)
python qc/compare_datasets.py /old /new "" 20

# Full QC on all files
python qc/compare_datasets.py /old /new

# QC specific variables only
python qc/compare_datasets.py /old /new tasmax,pr,hursmin
```

---

## Further Documentation

- **[pipeline_comparison_analysis.md](pipeline_comparison_analysis.md)** - Comprehensive technical analysis of all differences
- **qc/bias_corrected/original_vs_fully_processed.txt** - Raw QC results showing specific failures
- **Original dataset:** [Zenodo v7783759](https://zenodo.org/records/7783759)

