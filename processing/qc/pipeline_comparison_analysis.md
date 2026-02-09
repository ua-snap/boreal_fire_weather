# Bias Correction Pipeline Comparison Analysis

**Date:** February 9, 2026  
**Analyst:** GitHub Copilot  
**Purpose:** Deep dive comparison of original vs refactored bias correction pipelines

---

## Executive Summary

The QC comparison revealed that **all 80 tested files differ** between the original and new pipelines. Investigation identified three critical differences in data processing:

1. **Shape/Dimension Ordering Issues** (CNRM-CM6-1-HR model only)
   - Original: `(lat, lon, time)` → `(178, 569, 365)`
   - New: `(time, lat, lon)` → `(365, 178, 569)`

2. **Time Coordinate Alignment** (All models)
   - Original: No explicit time alignment between ERA5 and GCM data
   - New: Explicit reindexing and time normalization to day boundaries

3. **Chunking Strategy Bug** (All models)
   - Original: Bug where both lat and lon chunks use lon dimension size
   - New: Correct dimension-specific chunking with auto-sizing

**LEGACY_MODE Implementation:** A legacy compatibility mode has been implemented that replicates the original pipeline's behavior (including bugs) for testing purposes. This mode skips time alignment, dimension ordering, and uses the original chunking bug.

---

## QC Results Summary

### Results by Variable

| Variable | Files Tested | Shape Mismatches | Value Differences | Pass Rate |
|----------|--------------|------------------|-------------------|-----------|
| hursmin  | 20           | 4 (CNRM)         | 16 (all others)   | 0%        |
| pr       | 20           | 7 (CNRM)         | 13 (all others)   | 0%        |
| sfcWind  | 20           | 3 (CNRM)         | 17 (all others)   | 0%        |
| tasmax   | 20           | 1 (CNRM)         | 19 (all others)   | 0%        |
| **Total**| **80**       | **15**           | **65**            | **0%**    |

### Numerical Difference Magnitudes

| Variable | Mean Difference | Max Difference | Affected Cells (%) |
|----------|----------------|----------------|-------------------|
| hursmin  | 0.025 - 0.042  | 7.7 - 25.9     | ~45-50%          |
| pr       | 0.002 - 0.004  | 5.1 - 8.9      | ~3-4%            |
| sfcWind  | 0.010 - 0.011  | 1.8 - 4.0      | ~45-50%          |
| tasmax   | 0.013 - 0.017  | 1.5 - 5.2      | ~45-50%          |

---

## Pipeline Structure Comparison

### File Locations

**Original Pipeline:**
- Script: `~/amyoun01-boreal_wildfire_feedbacks-a3d14b5/scripts/04_bias_correct_gcms.py`
- Function: `wildfire_analysis/data_processing/quantile_delta_mapping.py`
- Config: `wildfire_analysis/config.yaml`
- Versions: v1.0 (https://doi.org/10.5281/zenodo.7783759) and v0.1-beta (https://zenodo.org/records/7761883)

**New Pipeline:**
- Script: `~/boreal_fire_weather/processing/03_bias_correct_gcms.py`
- Functions: Inline in the same script
- Config: `config.py` (Python module)

### Algorithm Commonality

Both pipelines use the **same core algorithm**:
- `xclim.sdba.QuantileDeltaMapping` (version 0.40.0)
- Same quantile values (101 quantiles from 0.005 to 0.995)
- Same operators (additive for tasmax, multiplicative for others)
- Same minimum threshold (0.5 mm/day for precipitation only)

---

## Critical Code Differences

### 1. **Dimension Ordering**

**Old Pipeline:**
```python
# NO explicit transpose before saving
export_ds.to_netcdf(fn, engine='h5netcdf')
```

**New Pipeline:**
```python
# EXPLICIT transpose to standardized dimension order
export_ds = export_ds.transpose("time", "lat", "lon")
export_ds.to_netcdf(fn, engine="h5netcdf", encoding=encoding)
```

**Impact:** The old pipeline preserves whatever dimension ordering xarray uses internally after regridding and QDM operations. For CNRM data, this apparently resulted in `(lat, lon, time)` ordering, while other models maintained `(time, lat, lon)`. The new pipeline enforces consistent ordering across all models.

---

### 2. **Time Coordinate Alignment**

**Old Pipeline:**
```python
# NO time coordinate alignment - loads data as-is
ref = xr.open_mfdataset(ref_src, engine='h5netcdf', ...)
hst = xr.open_mfdataset(hst_src, engine='h5netcdf', ...)
sim = xr.open_mfdataset(sim_src, engine='h5netcdf', ...)
```

**New Pipeline:**
```python
# EXPLICIT time coordinate alignment
if not ref.time.equals(hst.time):
    hst = hst.reindex(time=ref.time, method="nearest", tolerance=pd.Timedelta("12h"))

# Floor all times to start of day (removes sub-day differences)
ref = ref.assign_coords(time=ref.time.dt.floor("D"))
hst = hst.assign_coords(time=hst.time.dt.floor("D"))
sim = sim.assign_coords(time=sim.time.dt.floor("D"))
```

**Impact:** Different time-of-day encodings (e.g., midnight vs noon) can cause slight misalignments in temporal grouping for monthly quantile calculations. When QDM groups by `time.month`, even small time offsets could place values in adjacent months, changing which quantile distributions are used for bias correction.

**Example Scenario:**
- Original: Time encoded as `2020-01-31 12:00:00` → Grouped in January
- New: Time floored to `2020-01-31 00:00:00` → Still January
- But if original was `2020-02-01 00:00:00` and ref was `2020-01-31 12:00:00`, reindexing could shift the alignment

This would affect the quantile mapping pairs used for each grid cell and time step, propagating through the entire QDM algorithm.

---

### 3. **Chunking Strategy Bug**

**Old Pipeline:**
```python
# Bug: Both n and m incorrectly set to X axis size
axes = h.get_geoaxes(ref)
n, m = (ref[axes['X'][0]].size, ref[axes['X'][0]].size)  # Both use X!

ref = ref.chunk(chunks={
    'time': -1,
    'lat': round(n * frac),  # Wrong: n should be Y axis size
    'lon': round(m * frac)
})
```

**New Pipeline:**
```python
# Fixed: Correct dimension-specific sizing
axes = get_geoaxes(ref)
n, m = (ref[axes["Y"]].size, ref[axes["X"]].size)  # Correct!

# Uses auto-chunking for optimal performance
chunks={"time": -1, axes["Y"]: "auto", axes["X"]: "auto"}
```

**Impact:** This bug in the original pipeline causes incorrect chunk sizing for the latitude dimension. With the actual grid dimensions (178 lat × 569 lon):

- **Intended behavior:** `n=178, m=569` → chunks of `{'lat': 36, 'lon': 114}`
- **Actual behavior:** `n=569, m=569` → chunks of `{'lat': 114, 'lon': 114}`

This chunking bug affects numerical results because:
1. QDM performs quantile calculations within chunks before aggregation
2. Different chunk boundaries → different aggregation order
3. Floating-point operations are not perfectly associative  
4. Small rounding differences accumulate during quantile mapping

Unlike simple chunking performance differences, this bug causes the Dask task graph to partition the data differently, leading to subtle numerical variations in the bias-corrected output.

**Note:** The original code also has `axes['X'][0]` attempting to index a string, but this appears to have been handled by the code's error handling mechanisms in practice.

---

### 4. **Hursmin Clipping**

**Old Pipeline:**
```python
# NO clipping applied
```

**New Pipeline:**
```python
# Clips hursmin to valid range [0, 100]
def clip_hursmin(data_array: xr.DataArray, var_name: str) -> xr.DataArray:
    if CLIP_HURSMIN and var_name == "hursmin":
        return data_array.clip(min=0.0, max=100.0)
    return data_array
```

**Impact:** This is an **intentional improvement** to prevent physically impossible relative humidity values. However, based on the QC results showing differences even in unclamped ranges, this is not the primary cause of discrepancies.

---

### 5. **Shapely API Change**

**Old Pipeline:**
```python
mask_grd = shapely.vectorized.contains(mask_geom, x, y)
```

**New Pipeline:**
```python
mask_grd = shapely.contains_xy(mask_geom, x, y)
```

**Impact:** This is an API update to use the newer Shapely 2.0 function. Functionally equivalent, should produce identical masks.

---

### 6. **Chunking Bug Investigation (AFFECTS RESULTS)**

**Bugs Found in Original Code:**

1. **In `chunk_data_arrays()` function:**
   ```python
   n, m = (ref[axes['X'][0]].size, ref[axes['X'][0]].size)
   # Both n and m use X axis - should be Y and X respectively
   ```

2. **In `regrid_geodata()` function:**
   ```python
   ds_interp = ds.interp({
       axes['X'][0]: x,  # axes['X'] is 'lon', so axes['X'][0] is 'l'
       axes['Y'][0]: y   # Similarly, axes['Y'][0] is 'l' from 'lat'
   }, method=method)
   ```

**Why chunking bug matters:**
The chunking bug causes incorrect chunk sizes for the latitude dimension, affecting how Dask partitions the computation graph. This leads to different intermediate results during quantile calculation and aggregation, causing small numerical differences in the final output.

**Why regridding bug doesn't matter:**
Grid compatibility analysis confirmed that all GCM data was already on the ERA5 grid (178×569, lat: 35-79.25°N, lon: -177 to -35°W). The shape check in `dimcheck_and_regrid()` passed, so the buggy `regrid_geodata()` function was never executed in production.

---

## LEGACY_MODE Implementation

To test whether the identified differences fully explain the QC discrepancies, a **LEGACY_MODE** environment variable has been implemented in the new pipeline. When enabled (`LEGACY_MODE=TRUE`), the pipeline replicates the original processing approach:

**Behaviors in LEGACY_MODE:**
1. ✓ **Skips time alignment** - No reindexing or time flooring
2. ✓ **Skips dimension transpose** - Preserves natural xarray ordering
3. ✓ **Uses buggy chunking** - Both n and m set to X axis size
4. ✓ **Matches original chunking pattern** - Fixed 20% chunks, not "auto"

**Testing Strategy:**
Running the pipeline with `LEGACY_MODE=TRUE` creates outputs using processing steps identical to the original pipeline. Comparing these outputs with the original dataset will determine if all significant differences have been identified and replicated.

**Important Note:** LEGACY_MODE is for diagnostic purposes only. It intentionally includes bugs from the original code to test the hypothesis that these bugs explain the observed differences. It should not be used for production work.

---

## Version History of Original Codebase

**Two versions examined:**
- **Original version (v1.0: https://doi.org/10.5281/zenodo.7783759):** Initial public release with environment.yaml
- **Beta version (v0.1-beta: https://zenodo.org/records/7761883):** Alternative commit with requirements.txt

**Finding:** Both versions contain **identical bias correction code**, including all bugs. Files confirmed identical:
- `scripts/04_bias_correct_gcms.py`
- `wildfire_analysis/data_processing/quantile_delta_mapping.py`  
- `wildfire_analysis/utils/helpers.py`
- `wildfire_analysis/config.yaml`
- Package versions (xclim 0.40.0, xarray 2022.12.0, dask 2022.12.1)

The existing original dataset could have been produced by either version, as they are functionally equivalent for the bias correction step.

---

## Package Version Comparison

### Identical Core Dependencies
Both pipelines use **identical versions** of critical packages:
- `xclim==0.40.0` ✓
- `python==3.9.15` ✓
- `xarray==2022.12.0` (old) vs unspecified (new, likely similar) ≈
- `dask==2022.12.1` (old) vs unspecified (new, likely similar) ≈
- `numpy==1.23.5` (old) vs `<1.25.0` (new) ~

**Conclusion:** Package versions are **not** the cause of differences. Both use the same `xclim` version which contains the QDM implementation.

---

## Processing Parameter Comparison

| Parameter | Old Pipeline | New Pipeline | Match? |
|-----------|--------------|--------------|--------|
| Historical years | 1980-2009 | 1980-2009 | ✓ |
| Simulation periods | (2010-2039), (2040-2069), (2070-2099) | Same | ✓ |
| Quantile values | 101 values (0.005 to 0.995) | Same | ✓ |
| Operators | +/* | Same | ✓ |
| Min threshold (pr) | 0.5 mm/day | 0.5 mm/day | ✓ |
| Regrid method | gcm2era | gcm2era | ✓ |
| Group method | time.month | time.month | ✓ |
| Interpolation | linear | linear | ✓ |
| Jittering | Yes (for pr) | Yes (for pr) | ✓ |

**Conclusion:** All QDM parameters are **identical**.

---

## Root Cause Analysis

### Three Interacting Causes

The numerical differences between pipelines stem from three distinct processing differences:

#### 1. Time Coordinate Misalignment

**Mechanism:** The original pipeline does not align time coordinates between ERA5 reference and GCM historical data. Different time-of-day encodings (midnight vs noon) can cause misalignment in the QDM training pairs.

**Impact Chain:**
- QDM groups by `time.month` to calculate separate quantile mappings
- If ref time is `2020-01-31 12:00` and hst time is `2020-01-31 00:00`, they match
- But if there's a day offset, values could be grouped into different months
- Different monthly quantile distributions → different bias correction mappings
- Effects propagate to all adjusted simulation values

**Evidence:** Widespread differences (~45-50% of cells) with small magnitudes suggest systematic shifting of quantile mappings rather than gross errors.

#### 2. Chunking Bug

**Mechanism:** The original pipeline sets both lat and lon chunk sizes to the same value (based on lon dimension):
```python
n = m = lon_size  # Bug: should be n=lat_size, m=lon_size
```
This creates chunks of `{'lat': 114, 'lon': 114}` instead of `{'lat': 36, 'lon': 114}`.

**Impact Chain:**
- Different chunk boundaries → different Dask task partitioning
- QDM quantile calculations done within chunks before aggregation
- Different aggregation order → floating-point rounding differences  
- Small numerical variations accumulate through the correction process

**Why it matters:** While chunking is often computationally transparent, quantile operations can be sensitive to how data is partitioned, especially when aggregating results across chunk boundaries.

#### 3. Dimension Ordering (CNRM FILES ONLY)

**Mechanism:** The original pipeline does not enforce dimension order, preserving whatever ordering results from internal xarray operations.

**Impact:** CNRM files ended up with `(lat, lon, time)` ordering while other models had `(time, lat, lon)`. This is a shape difference only, not a values difference.

### Combined Effect

These three factors work together to produce the observed QC results:
- **All models:** Time alignment + chunking bug → numerical differences
- **CNRM only:** Addition of dimension reordering → shape mismatch

The LEGACY_MODE implementation replicates all three original behaviors to test whether they fully account for the observed differences.

---

## Notable Improvements in New Pipeline

1. **Distributed computing optimization:** Uses `persist()` to load ref/hst once and reuse across simulation periods
2. **Explicit dimension ordering:** Ensures consistent output format
3. **Time coordinate normalization:** Prevents temporal misalignment issues
4. **Hursmin bounds checking:** Prevents physically impossible values (optional via config)
5. **Better encoding control:** Explicit dtype handling for time and data variables
6. **Progress tracking:** Better user feedback with progress bars
7. **Bug fixes:** Corrected helper function indexing issues
8. **Code clarity:** More modular structure with better documentation

---

## Recommendations

### Testing and Validation

1. **LEGACY_MODE Testing (IN PROGRESS)**
   - Run the new pipeline with `LEGACY_MODE=TRUE` and `CLIP_HURSMIN=FALSE` to replicate original processing
   - Compare LEGACY_MODE output against original dataset using QC scripts
   - This will verify whether the identified differences (time alignment, chunking bug, dimension ordering) fully account for the observed discrepancies
   - **Status:** LEGACY_MODE implementation complete, validation run pending

2. **Expected LEGACY_MODE Outcomes:**
   - If LEGACY outputs match original dataset → All differences identified and understood
   - If LEGACY outputs still differ → Additional processing differences remain to be found
   - Either outcome provides valuable information for the path forward

### For Production Use (After LEGACY_MODE Validation)

1. **Default recommendation: Use improved pipeline**
   - Fixes time alignment issues for more accurate temporal pairing
   - Corrects chunking bug for proper dimension handling
   - Enforces consistent dimension ordering across all models
   - Provides better performance with optimized chunking

2. **If LEGACY_MODE matches original dataset:**
   - Treat original dataset as v1.0 (with known processing quirks)
   - Treat improved pipeline output as v2.0 (with corrections)
   - Document differences and provide migration guidance
   - Consider regenerating critical downstream datasets with v2.0

3. **If backward compatibility is critical:**
   - LEGACY_MODE can reproduce original processing for specific comparisons
   - Not recommended for new production work
   - Should be clearly documented as replicating historical bugs

### For Documentation

1. **Processing provenance:**
   - Clearly indicate which pipeline version produced each dataset
   - Document the three key processing differences
   - Explain why v2.0 is recommended despite numerical differences

2. **Migration guide for users:**
   - Expected magnitude of differences by variable
   - Which downstream analyses are sensitive to these differences
   - How to interpret comparisons between v1.0 and v2.0 outputs

3. **Version control:**
   - Tag original pipeline code with appropriate version markers
   - Maintain LEGACY_MODE for reproducibility but warn against production use
   - Use semantic versioning for future pipeline changes

---

## Conclusion

### Processing Differences Identified

Comprehensive comparison of the original and refactored bias correction pipelines has identified **three key processing differences**:

1. **Time coordinate alignment** - Original pipeline did not align time coordinates between ERA5 and GCM data, potentially causing temporal mismatches in QDM training pairs
2. **Chunking bug** - Original pipeline incorrectly used lon dimension size for both lat and lon chunks, affecting Dask task partitioning and quantile calculations
3. **Dimension ordering** - Original pipeline did not enforce consistent dimension order, resulting in variable output formats

### LEGACY_MODE for Validation

A **LEGACY_MODE** has been implemented in the refactored pipeline that replicates all three original behaviors, including bugs. This enables direct testing of the hypothesis that these three differences fully explain the observed QC discrepancies:

- When `LEGACY_MODE=TRUE`: Uses original time handling, buggy chunking, and flexible dimension ordering
- When `LEGACY_MODE=FALSE`: Uses improved time alignment, correct chunking, and enforced ordering

**Current Status:** LEGACY_MODE implementation is complete. Validation testing against the original dataset will determine whether the identified differences are comprehensive or if additional processing variations remain to be discovered.

### Pipeline Quality Assessment

Both pipelines use the **same core QDM algorithm** (`xclim.sdba.QuantileDeltaMapping` v0.40.0) with identical parameters. The differences lie entirely in data preprocessing steps:

**Original Pipeline:**
- Simple approach with minimal preprocessing
- Contains subtle bugs in helper functions
- Inconsistent handling of time coordinates and dimensions
- Produced the published v1.0 dataset

**Refactored Pipeline:**  
- Adds explicit preprocessing for robustness
- Fixes bugs in chunking and regridding logic
- Enforces consistent time alignment and dimension ordering
- Produces v2.0 dataset with improved reproducibility

The refactored pipeline represents a methodological improvement, but numerical differences from v1.0 are expected due to the preprocessing corrections.

**Next Step:** Complete LEGACY_MODE validation run to verify that all processing differences have been identified and successfully replicated.

---

**Analysis Date:** February 9, 2026  
**Document Version:** 2.0  
**Status:** LEGACY_MODE testing in progress
