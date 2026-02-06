# Quality Control (QC) for Boreal Fire Weather Pipeline

## Purpose

This QC process validates that the current processing pipeline produces results identical to the reference version ([Zenodo v7783759](https://zenodo.org/records/7783759)). It is designed to be run after executing the pipeline with `CLIP_HURSMIN=FALSE`.

## Overview

The `compare_datasets.py` script compares NetCDF files with identical names between two directories and checks that their values are identical. It works on both bias-corrected GCM outputs and CFFDRS outputs, regardless of directory structure (as long as there are no duplicate filenames within each directory).

## Features

- Compares all data variables in matching NetCDF files
- Optional filtering by variable names
- Optional random sampling to limit number of files checked per variable
- Detailed difference reporting (max, mean, count of different values)
- Clear pass/fail status for each file and overall summary

## Usage

### Basic syntax

```bash
python compare_datasets.py <old_dir> <new_dir> [variables] [max_files]
```

**Arguments:**
- `old_dir` - Path to reference/old dataset directory
- `new_dir` - Path to new dataset directory to validate
- `variables` - (Optional) Comma-separated variable names (e.g., `tasmax,pr,hursmin`)
- `max_files` - (Optional) Maximum files to randomly sample per variable

### Examples

```bash
# Compare all files in both directories
python compare_datasets.py /data/old /data/new

# Compare only specific variables
python compare_datasets.py /data/old /data/new tasmax,pr,hursmin

# Compare 20 randomly selected files per variable
python compare_datasets.py /data/old /data/new tasmax,pr 20

# Compare 10 files across all variables
python compare_datasets.py /data/old /data/new "" 10
```

## Output

Results are printed to the terminal with:
- Per-file pass/fail status
- Detailed difference statistics for failed comparisons
- Summary counts by variable
- Overall QC pass/fail result

The script exits with code 0 if all files match, or code 1 if any differences are found.
