#!/usr/bin/env python
"""
NaN Pattern Comparison Tool for Boreal Fire Weather Data

This script identifies and visualizes NaN pattern differences between datasets.
It groups files by model+variable and selects representative examples to plot.

Usage:
    python compare_nans.py <old_dir> <new_dir> <output_dir> [variables] [max_files] [--shapefile path]

Arguments:
    old_dir     : Path to directory containing reference/old dataset files
    new_dir     : Path to directory containing new dataset files
    output_dir  : Directory where output HTML report will be saved
    variables   : (Optional) Comma-separated list of variables (e.g., "ffmc,dmc,dc")
    max_files   : (Optional) Maximum number of files to randomly sample per variable
    --shapefile : (Optional) Custom shapefile path for boundary overlay

Examples:
    # Analyze all matching files
    python compare_nans.py /data/old /data/new /output/nan_qc

    # Analyze specific variables
    python compare_nans.py /data/old /data/new /output/nan_qc "ffmc,dmc,dc"

    # Sample 50 files per variable
    python compare_nans.py /data/old /data/new /output/nan_qc "" 50

    # Custom shapefile
    python compare_nans.py /data/old /data/new /output/nan_qc --shapefile ../shp/custom.shp
"""

import sys
import argparse
from pathlib import Path
import xarray as xr
import numpy as np
from collections import defaultdict
import random
import subprocess
import json
import pandas as pd


def find_matching_files(dir1: Path, dir2: Path) -> dict:
    """Find files with identical names in two directories."""
    files1 = {f.name: f for f in dir1.rglob("*.nc")}
    files2 = {f.name: f for f in dir2.rglob("*.nc")}

    common_names = set(files1.keys()) & set(files2.keys())

    if not common_names:
        print(f"ERROR: No matching files found between directories")
        print(f"  Directory 1: {dir1}")
        print(f"  Directory 2: {dir2}")
        sys.exit(1)

    return {name: (files1[name], files2[name]) for name in common_names}


def group_files_by_variable(file_pairs: dict) -> dict:
    """Group file pairs by variable name."""
    var_groups = defaultdict(list)

    for filename, (path1, path2) in file_pairs.items():
        # Extract variable name (assumes format: var_model_year.nc or similar)
        parts = filename.split("_")
        if len(parts) >= 1:
            var_name = parts[0]
            var_groups[var_name].append((filename, path1, path2))

    return dict(var_groups)


def extract_model_from_filename(filename: str) -> str:
    """Extract model name from filename."""
    # Common patterns: var_MODEL_year.nc, cffdrs_MODEL_year.nc, etc.
    parts = filename.replace(".nc", "").split("_")
    if len(parts) >= 2:
        return parts[1]  # Assumes second part is model name
    return "unknown"


def analyze_nan_patterns(ref_path: Path, new_path: Path) -> dict:
    """
    Analyze NaN pattern differences across all timesteps.

    Returns dict with:
        - total_nan_diff: Total count of cells with different NaN patterns across all timesteps
        - worst_timestep: Index of timestep with most NaN differences
        - worst_timestep_count: NaN diff count at worst timestep
        - n_timesteps: Total number of timesteps
        - variables: List of variables in the dataset
    """
    try:
        ds_ref = xr.open_dataset(ref_path)
        ds_new = xr.open_dataset(new_path)

        # Get all data variables
        data_vars = [v for v in ds_ref.data_vars if v in ds_new.data_vars]

        if not data_vars:
            ds_ref.close()
            ds_new.close()
            return None

        n_times = len(ds_ref.time) if "time" in ds_ref.dims else 1

        # Analyze across all timesteps and variables
        total_nan_diff = 0
        worst_timestep = 0
        worst_timestep_count = 0
        timestep_counts = []

        for t_idx in range(n_times):
            timestep_diff = 0

            for var in data_vars:
                if n_times > 1:
                    data_ref = ds_ref[var].isel(time=t_idx).values
                    data_new = ds_new[var].isel(time=t_idx).values
                else:
                    data_ref = ds_ref[var].values
                    data_new = ds_new[var].values

                nan_ref = np.isnan(data_ref)
                nan_new = np.isnan(data_new)
                nan_diff = np.sum((nan_ref & ~nan_new) | (~nan_ref & nan_new))

                timestep_diff += nan_diff

            timestep_counts.append(timestep_diff)
            total_nan_diff += timestep_diff

            if timestep_diff > worst_timestep_count:
                worst_timestep_count = timestep_diff
                worst_timestep = t_idx

        result = {
            "total_nan_diff": int(total_nan_diff),
            "worst_timestep": int(worst_timestep),
            "worst_timestep_count": int(worst_timestep_count),
            "n_timesteps": int(n_times),
            "variables": data_vars,
            "mean_nan_diff_per_timestep": (
                float(total_nan_diff / n_times) if n_times > 0 else 0
            ),
            "timesteps_with_diffs": int(sum(1 for c in timestep_counts if c > 0)),
        }

        ds_ref.close()
        ds_new.close()

        return result

    except Exception as e:
        print(f"    ERROR analyzing {ref_path.name}: {str(e)}")
        return None


def select_representative_files(analysis_results: dict, max_plots: int = 10) -> list:
    """
    Select representative files to visualize based on model+variable grouping.

    Returns list of tuples: (filename, ref_path, new_path, analysis_info, group_label)
    """
    # Group by model + variable
    groups = defaultdict(list)

    for filename, info in analysis_results.items():
        if info["analysis"]["total_nan_diff"] == 0:
            continue  # Skip files with no differences

        model = extract_model_from_filename(filename)
        # Use first variable as representative (most files have one primary variable)
        var = (
            info["analysis"]["variables"][0]
            if info["analysis"]["variables"]
            else "unknown"
        )
        group_key = f"{model}_{var}"

        groups[group_key].append((filename, info))

    # Select one file per group with maximum NaN difference
    selected = []
    for group_key, files in sorted(groups.items()):
        # Sort by total NaN diff and take the one with most differences
        files_sorted = sorted(
            files, key=lambda x: x[1]["analysis"]["total_nan_diff"], reverse=True
        )
        filename, info = files_sorted[0]

        selected.append(
            {
                "filename": filename,
                "ref_path": info["ref_path"],
                "new_path": info["new_path"],
                "analysis": info["analysis"],
                "group": group_key,
                "group_size": len(files),
            }
        )

    # Sort by total NaN difference and take top max_plots
    selected_sorted = sorted(
        selected, key=lambda x: x["analysis"]["total_nan_diff"], reverse=True
    )

    return selected_sorted[:max_plots]


def generate_visualization(
    selected_files: list, output_dir: Path, shapefile_path: Path, qc_utils_path: str
):
    """Generate HTML visualization report using papermill."""
    if not selected_files:
        print("\nNo files with NaN differences to visualize.")
        return

    print("\n" + "=" * 80)
    print("GENERATING NaN PATTERN VISUALIZATION")
    print("=" * 80)
    print(f"Creating report for {len(selected_files)} representative file(s)...\n")

    # Path to notebook template
    template_nb = Path(__file__).parent / "visualize_nan_comparison.ipynb"

    if not template_nb.exists():
        print(f"ERROR: Notebook template not found at {template_nb}")
        return

    # Create output directories
    output_dir.mkdir(parents=True, exist_ok=True)
    html_output_dir = output_dir / "html"
    html_output_dir.mkdir(parents=True, exist_ok=True)

    # Output paths
    output_nb = output_dir / "nan_pattern_comparison.ipynb"
    output_html = html_output_dir / "nan_pattern_comparison.html"

    try:
        # Execute notebook with papermill
        result = subprocess.run(
            [
                "papermill",
                str(template_nb),
                str(output_nb),
                "-p",
                "qc_utils_path",
                qc_utils_path,
                "-p",
                "shapefile_path",
                str(shapefile_path),
                "-y",
                f"selected_files: {json.dumps(selected_files)}",
            ],
            capture_output=True,
            text=True,
            timeout=600,
        )

        if result.returncode != 0:
            print(f"ERROR executing notebook: {result.stderr}")
            return

        # Convert to HTML
        html_filename = "nan_pattern_comparison.html"
        result_html = subprocess.run(
            [
                "jupyter",
                "nbconvert",
                "--to",
                "html",
                "--output-dir",
                str(html_output_dir),
                "--output",
                html_filename,
                str(output_nb),
            ],
            capture_output=True,
            text=True,
            timeout=120,
        )

        if result_html.returncode != 0:
            print(f"WARNING: Failed to convert to HTML: {result_html.stderr}")
            print(f"Notebook saved at: {output_nb}")
        else:
            print(f"\n✓ HTML report generated: {output_html}")

        print("=" * 80)

    except subprocess.TimeoutExpired:
        print("ERROR: Visualization timed out")
    except Exception as e:
        print(f"ERROR: {str(e)}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare NaN patterns between two dataset directories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument("old_dir", help="Reference/old dataset directory")
    parser.add_argument("new_dir", help="New dataset directory")
    parser.add_argument("output_dir", help="Output directory for reports")
    parser.add_argument(
        "variables",
        nargs="?",
        default="",
        help="Comma-separated variable names (optional)",
    )
    parser.add_argument(
        "max_files",
        nargs="?",
        type=int,
        default=None,
        help="Max files to sample per variable (optional)",
    )
    parser.add_argument(
        "--shapefile",
        default="shp/ecos.shp",
        help="Path to shapefile for boundary overlay",
    )

    args = parser.parse_args()

    # Convert to Path objects
    old_dir = Path(args.old_dir).resolve()
    new_dir = Path(args.new_dir).resolve()
    output_dir = Path(args.output_dir)
    shapefile_path = Path(args.shapefile)
    qc_utils_path = str(Path(__file__).parent.resolve())

    # Validate input directories
    if not old_dir.exists():
        print(f"ERROR: Reference directory not found: {old_dir}")
        sys.exit(1)
    if not new_dir.exists():
        print(f"ERROR: New dataset directory not found: {new_dir}")
        sys.exit(1)

    print(f"Reference directory: {old_dir}")
    print(f"New directory: {new_dir}")
    print(f"Output directory: {output_dir}")

    # Find matching files
    print("\nFinding matching files...")
    file_pairs = find_matching_files(old_dir, new_dir)
    print(f"Found {len(file_pairs)} matching file pairs")

    # Filter by variables if specified
    if args.variables:
        var_list = [v.strip() for v in args.variables.split(",")]
        var_groups = group_files_by_variable(file_pairs)

        filtered_pairs = {}
        for var in var_list:
            if var in var_groups:
                for filename, path1, path2 in var_groups[var]:
                    filtered_pairs[filename] = (path1, path2)

        file_pairs = filtered_pairs
        print(f"Filtered to {len(file_pairs)} files for variables: {var_list}")

    # Sample files if max_files specified
    if args.max_files and len(file_pairs) > args.max_files:
        var_groups = group_files_by_variable(file_pairs)
        sampled_pairs = {}

        for var, files in var_groups.items():
            sample_size = min(args.max_files, len(files))
            sampled = random.sample(files, sample_size)
            for filename, path1, path2 in sampled:
                sampled_pairs[filename] = (path1, path2)

        file_pairs = sampled_pairs
        print(f"Sampled {len(file_pairs)} files (max {args.max_files} per variable)")

    # Analyze NaN patterns
    print("\nAnalyzing NaN patterns...")
    analysis_results = {}

    for i, (filename, (ref_path, new_path)) in enumerate(file_pairs.items(), 1):
        print(f"  [{i}/{len(file_pairs)}] {filename}...", end=" ")

        analysis = analyze_nan_patterns(ref_path, new_path)

        if analysis:
            analysis_results[filename] = {
                "ref_path": str(ref_path),
                "new_path": str(new_path),
                "analysis": analysis,
            }

            if analysis["total_nan_diff"] > 0:
                print(f"✗ {analysis['total_nan_diff']:,} NaN differences")
            else:
                print("✓ Identical")
        else:
            print("ERROR")

    # Summary
    files_with_diffs = sum(
        1
        for info in analysis_results.values()
        if info["analysis"]["total_nan_diff"] > 0
    )

    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print(f"Total files analyzed: {len(analysis_results)}")
    print(f"Files with NaN differences: {files_with_diffs}")
    print(f"Files identical: {len(analysis_results) - files_with_diffs}")

    if files_with_diffs > 0:
        # Select representative files for visualization
        selected = select_representative_files(analysis_results, max_plots=10)

        print(f"\nSelected {len(selected)} representative files for visualization")
        print(f"  (1 per model+variable group, max 10 total)")

        # Generate visualization
        generate_visualization(selected, output_dir, shapefile_path, qc_utils_path)
    else:
        print("\nNo NaN pattern differences found. No visualization generated.")

    print(f"{'='*80}\n")


if __name__ == "__main__":
    main()
