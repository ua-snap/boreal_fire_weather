#!/usr/bin/env python
"""
Quality Control Script for Boreal Fire Weather Processing Pipeline

This script compares NetCDF datasets from two different processing pipeline runs.
By default, it generates visual HTML reports with maps showing differences.
Use --text-only for traditional text-based QC reports.

Usage:
    python compare_datasets.py <old_dir> <new_dir> <output_dir> [variables] [max_files] [--text-only] [-o output.txt]

Arguments:
    old_dir     : Path to directory containing reference/old dataset files
    new_dir     : Path to directory containing new dataset files to validate
    output_dir  : Directory where QC outputs (notebooks/html) will be saved
    variables   : (Optional) Comma-separated list of variables to check (e.g., "tasmax,pr,hursmin")
                  If not provided, all matching files are compared
    max_files   : (Optional) Maximum number of files to randomly sample per variable
                  If not provided, all matching files are compared. Not recommended!
    --text-only : (Optional) Generate text report only, skip visualizations
    -o, --output: (Optional) Save text report to file instead of stdout
    --shapefile : (Optional) Custom shapefile path for boundary overlay in visualizations

Examples:
    # Generate visual HTML reports (default)
    python compare_datasets.py /data/old /data/new /output/qc

    # Text report to stdout
    python compare_datasets.py /data/old /data/new /output/qc --text-only

    # Text report saved to file
    python compare_datasets.py /data/old /data/new /output/qc --text-only -o qc_report.txt

    # Compare only specific variables with visualizations
    python compare_datasets.py /data/old /data/new /output/qc "tasmax,pr,hursmin"

    # Compare 20 randomly selected files per variable with visualizations
    python compare_datasets.py /data/old /data/new /output/qc "tasmax,pr" 20

    # Compare 10 files of all variables with text output
    python compare_datasets.py /data/old /data/new /output/qc "" 10 --text-only
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


def find_matching_files(dir1: Path, dir2: Path) -> dict:
    """
    Find files with identical names in two directories.

    Parameters
    ----------
    dir1 : Path
        First directory to search
    dir2 : Path
        Second directory to search

    Returns
    -------
    dict : Dictionary mapping filenames to tuples of (path_in_dir1, path_in_dir2)
    """
    # Get all NetCDF files recursively
    files1 = {f.name: f for f in dir1.rglob("*.nc")}
    files2 = {f.name: f for f in dir2.rglob("*.nc")}

    # Find common filenames
    common_names = set(files1.keys()) & set(files2.keys())

    if not common_names:
        print(f"ERROR: No matching files found between directories")
        print(f"  Directory 1: {dir1}")
        print(f"  Directory 2: {dir2}")
        sys.exit(1)

    return {name: (files1[name], files2[name]) for name in common_names}


def group_files_by_variable(file_pairs: dict) -> dict:
    """
    Group file pairs by variable name (extracted from filename).

    Parameters
    ----------
    file_pairs : dict
        Dictionary of filename: (path1, path2) pairs

    Returns
    -------
    dict : Dictionary mapping variable names to lists of file pairs
    """
    var_groups = defaultdict(list)

    for filename, (path1, path2) in file_pairs.items():
        # Extract variable name (assumes format: varname_*)
        var_name = filename.split("_")[0]
        var_groups[var_name].append((filename, path1, path2))

    return dict(var_groups)


def align_dimensions(ds1: xr.Dataset, ds2: xr.Dataset, var: str) -> tuple:
    """
    Attempt to align dimensions between two datasets if they have different ordering.

    Parameters
    ----------
    ds1 : xr.Dataset
        First dataset
    ds2 : xr.Dataset
        Second dataset
    var : str
        Variable name to align

    Returns
    -------
    tuple : (arr1, arr2, success, message)
    """
    arr1 = ds1[var]
    arr2 = ds2[var]

    # If shapes already match, return as-is
    if arr1.shape == arr2.shape:
        return arr1.values, arr2.values, True, ""

    # Try to align by dimension names
    dims1 = list(arr1.dims)
    dims2 = list(arr2.dims)

    # Check if they have the same dimensions, just in different order
    if set(dims1) == set(dims2):
        # Transpose arr2 to match arr1's dimension order
        try:
            arr2_aligned = arr2.transpose(*dims1)
            if arr1.shape == arr2_aligned.shape:
                return (
                    arr1.values,
                    arr2_aligned.values,
                    True,
                    f" (dimensions reordered from {dims2} to {dims1})",
                )
        except Exception as e:
            # If we can't align, return arrays as-is with failure flag
            return (
                arr1.values,
                arr2.values,
                False,
                f" (cannot align: {dims1} vs {dims2})",
            )
    # If we can't align, return arrays as-is with failure flag
    return (
        arr1.values,
        arr2.values,
        False,
        f" (different dimensions: {dims1} vs {dims2})",
    )


def compare_datasets(file1: Path, file2: Path) -> tuple:
    """
    Compare two NetCDF datasets and check if they are identical.

    Parameters
    ----------
    file1 : Path
        Path to first NetCDF file
    file2 : Path
        Path to second NetCDF file

    Returns
    -------
    tuple : (bool, str, dict) - (is_identical, message, stats)
    """
    try:
        ds1 = xr.open_dataset(file1)
        ds2 = xr.open_dataset(file2)
    except Exception as e:
        return False, f"Error opening files: {e}", {}

    # Get data variable names (exclude coordinates)
    data_vars1 = set(ds1.data_vars.keys())
    data_vars2 = set(ds2.data_vars.keys())

    if data_vars1 != data_vars2:
        ds1.close()
        ds2.close()
        return False, f"Data variables differ: {data_vars1} vs {data_vars2}", {}

    # Compare each data variable
    all_identical = True
    stats = {}
    messages = []

    for var in data_vars1:
        # Try to align dimensions if they differ
        arr1, arr2, aligned, align_msg = align_dimensions(ds1, ds2, var)

        # Check shapes match after alignment attempt
        if arr1.shape != arr2.shape:
            all_identical = False
            messages.append(
                f"  Variable '{var}': Shape mismatch ({arr1.shape} vs {arr2.shape}){align_msg}"
            )
            continue

        # Note if dimensions were reordered for comparison
        if aligned and align_msg:
            messages.append(f"  Variable '{var}': Note{align_msg}")

        # Handle NaN values
        nan_mask1 = np.isnan(arr1)
        nan_mask2 = np.isnan(arr2)

        if not np.array_equal(nan_mask1, nan_mask2):
            all_identical = False
            nan_diff = np.sum(nan_mask1 != nan_mask2)
            messages.append(
                f"  Variable '{var}': NaN patterns differ ({nan_diff} locations)"
            )

        # Compare non-NaN values
        valid_mask = ~nan_mask1 & ~nan_mask2
        if np.any(valid_mask):
            vals1 = arr1[valid_mask]
            vals2 = arr2[valid_mask]

            # Check for exact equality
            if not np.allclose(vals1, vals2, rtol=0, atol=0, equal_nan=True):
                all_identical = False
                max_diff = np.max(np.abs(vals1 - vals2))
                mean_diff = np.mean(np.abs(vals1 - vals2))
                n_diff = np.sum(vals1 != vals2)

                messages.append(
                    f"  Variable '{var}': Values differ\n"
                    f"    Max difference: {max_diff:.2e}\n"
                    f"    Mean difference: {mean_diff:.2e}\n"
                    f"    Number of different values: {n_diff}/{len(vals1)}"
                )

                stats[var] = {
                    "max_diff": float(max_diff),
                    "mean_diff": float(mean_diff),
                    "n_diff": int(n_diff),
                    "n_total": int(len(vals1)),
                }

    ds1.close()
    ds2.close()

    if all_identical:
        return True, "Datasets are identical", stats
    else:
        return False, "\n".join(messages), stats


def generate_visualizations(
    file_pairs_to_visualize: list, output_dir: Path, shapefile_path: Path
):
    """
    Generate visual comparison notebooks using papermill.

    Parameters
    ----------
    file_pairs_to_visualize : list
        List of (filename, ref_path, new_path, variables) tuples
    output_dir : Path
        Directory to save output notebooks and HTML
    shapefile_path : Path
        Path to shapefile for boundary overlay
    """
    if not file_pairs_to_visualize:
        print("No files to visualize.")
        return

    print("\n" + "=" * 80)
    print("GENERATING VISUAL COMPARISONS")
    print("=" * 80)
    print(f"Generating visualizations for {len(file_pairs_to_visualize)} file(s)...\n")

    # Path to notebook template and qc_utils.py
    template_nb = Path(__file__).parent / "visualize_comparison.ipynb"
    qc_utils_dir = str(
        Path(__file__).parent.resolve()
    )  # Absolute path to qc/ directory

    if not template_nb.exists():
        print(f"ERROR: Notebook template not found at {template_nb}")
        print("Skipping visualization generation.")
        return

    # Create output directory for notebooks
    nb_output_dir = output_dir / "notebooks"
    nb_output_dir.mkdir(parents=True, exist_ok=True)

    html_output_dir = output_dir / "html"
    html_output_dir.mkdir(parents=True, exist_ok=True)

    success_count = 0
    fail_count = 0

    for filename, ref_path, new_path, var_list in file_pairs_to_visualize:
        try:
            # Create output paths
            base_name = filename.replace(".nc", "")
            output_nb = nb_output_dir / f"{base_name}_comparison.ipynb"
            output_html = html_output_dir / f"{base_name}_comparison.html"

            print(f"  Processing: {filename}")

            # Prepare parameters for papermill
            parameters = {
                "ref_file": str(ref_path),
                "new_file": str(new_path),
                "variables": var_list,
                "shapefile_path": str(shapefile_path),
                "qc_utils_path": qc_utils_dir,
            }

            # Execute notebook with papermill
            result = subprocess.run(
                [
                    "papermill",
                    str(template_nb),
                    str(output_nb),
                    "-p",
                    "ref_file",
                    str(ref_path),
                    "-p",
                    "new_file",
                    str(new_path),
                    "-p",
                    "shapefile_path",
                    str(shapefile_path),
                    "-p",
                    "qc_utils_path",
                    qc_utils_dir,
                    "-y",
                    f"variables: {json.dumps(var_list)}",
                ],
                capture_output=True,
                text=True,
                timeout=600,  # 10 minute timeout
            )

            if result.returncode != 0:
                print(f"    ERROR executing notebook: {result.stderr}")
                fail_count += 1
                continue

            # Convert notebook to HTML
            html_filename = f"{base_name}_comparison.html"
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
                print(f"    WARNING: Failed to convert to HTML: {result_html.stderr}")
                print(f"    Notebook saved at: {output_nb}")
            else:
                print(f"    ✓ HTML report: {output_html}")

            success_count += 1

        except subprocess.TimeoutExpired:
            print(f"    ERROR: Visualization timed out for {filename}")
            fail_count += 1
        except Exception as e:
            print(f"    ERROR: {str(e)}")
            fail_count += 1

    print(f"\nVisualization summary: {success_count} succeeded, {fail_count} failed")
    print("=" * 80)


def run_qc(
    old_dir: str,
    new_dir: str,
    output_dir: str,
    variables: list = None,
    max_files: int = None,
    text_only: bool = False,
    text_output: str = None,
    shapefile_path: str = None,
):
    """
    Run quality control comparison between two dataset directories.

    Parameters
    ----------
    old_dir : str
        Path to reference/old dataset directory
    new_dir : str
        Path to new dataset directory to validate
    output_dir : str
        Directory where QC outputs (notebooks/html) will be saved
    variables : list, optional
        List of variable names to check. If None, all variables are checked
    max_files : int, optional
        Maximum number of files to randomly sample per variable. If None, all files are checked
    text_only : bool, optional
        If True, only generate text output (skip visualizations). Default: False
    text_output : str, optional
        Path to text output file. If None and text_only=True, prints to stdout
    shapefile_path : str, optional
        Path to shapefile for boundary overlay in visualizations
    """
    old_path = Path(old_dir)
    new_path = Path(new_dir)
    output_path = Path(output_dir)

    # Setup output stream
    output_file = None
    original_stdout = sys.stdout
    if text_output:
        output_file = open(text_output, "w")
        sys.stdout = output_file

    # Validate directories exist
    if not old_path.exists():
        print(f"ERROR: Old dataset directory does not exist: {old_dir}")
        sys.exit(1)
    if not new_path.exists():
        print(f"ERROR: New dataset directory does not exist: {new_dir}")
        sys.exit(1)

    # Create output directory if it doesn't exist
    output_path.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("BOREAL FIRE WEATHER PIPELINE - QUALITY CONTROL")
    print("=" * 80)
    print(f"Reference directory: {old_dir}")
    print(f"New directory:       {new_dir}")
    print(f"Output directory:    {output_dir}")
    if variables:
        print(f"Variables to check:  {', '.join(variables)}")
    else:
        print(f"Variables to check:  All")
    if max_files:
        print(f"Max files per var:   {max_files}")
    else:
        print(f"Max files per var:   All")
    print("=" * 80)
    print()

    # Find matching files
    print("Finding matching files...")
    file_pairs = find_matching_files(old_path, new_path)
    print(f"Found {len(file_pairs)} matching files\n")

    # Group by variable
    var_groups = group_files_by_variable(file_pairs)

    # Filter by requested variables
    if variables:
        var_groups = {
            var: files for var, files in var_groups.items() if var in variables
        }
        if not var_groups:
            print(f"ERROR: No files found for requested variables: {variables}")
            sys.exit(1)

    # Sample files if max_files specified
    if max_files:
        for var in var_groups:
            if len(var_groups[var]) > max_files:
                var_groups[var] = random.sample(var_groups[var], max_files)

    # Run comparisons
    total_files = sum(len(files) for files in var_groups.values())
    passed = 0
    failed = 0

    # Track files for visualization
    files_to_visualize = []

    print(
        f"Running QC on {total_files} files across {len(var_groups)} variable(s)...\n"
    )

    for var_name in sorted(var_groups.keys()):
        files = var_groups[var_name]
        print(f"Variable: {var_name} ({len(files)} files)")
        print("-" * 80)

        var_passed = 0
        var_failed = 0

        for filename, path1, path2 in sorted(files):
            is_identical, message, stats = compare_datasets(path1, path2)

            if is_identical:
                var_passed += 1
                passed += 1
                status = "✓ PASS"
            else:
                var_failed += 1
                failed += 1
                status = "✗ FAIL"

                # Collect failed files for visualization
                if not text_only:
                    # Determine which variables to visualize
                    if var_name == "cffdrs":
                        # For cffdrs files, extract all index variables
                        vars_to_viz = ["ffmc", "dmc", "dc", "isi", "bui", "fwi"]
                    else:
                        # For other files, use the variable name from the filename
                        vars_to_viz = [var_name]

                    files_to_visualize.append((filename, path1, path2, vars_to_viz))

            print(f"  {status}: {filename}")
            if not is_identical:
                print(message)
                print()

        print(f"  Summary: {var_passed} passed, {var_failed} failed")
        print()

    # Final summary
    print("=" * 80)
    print("QUALITY CONTROL SUMMARY")
    print("=" * 80)
    print(f"Total files compared: {total_files}")
    print(f"Passed:               {passed}")
    print(f"Failed:               {failed}")

    # Close text output file if needed
    if output_file:
        sys.stdout = original_stdout
        output_file.close()
        print(f"\nText report written to: {text_output}")

    # Generate visualizations for failed files (unless text_only mode)
    if not text_only and files_to_visualize:
        # Determine shapefile path
        if shapefile_path:
            shp_path = Path(shapefile_path)
        else:
            # Try default location relative to this script
            shp_path = Path(__file__).parent.parent / "shp" / "ecos.shp"

        if not shp_path.exists():
            print(f"\nWARNING: Shapefile not found at {shp_path}")
            print("Visualizations will be generated without shapefile overlay.")

        # Generate visualizations in the specified output directory
        generate_visualizations(files_to_visualize, output_path, shp_path)

    if failed == 0:
        if not output_file:  # Only print if not already in text file
            print("\n✓ All files identical - QC PASSED")
            print("=" * 80)
        return 0
    else:
        if not output_file:
            print(f"\n✗ {failed} file(s) differ - QC FAILED")
            print("=" * 80)
        return 1


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Compare NetCDF datasets between two processing pipeline runs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate visual HTML reports (default)
  python compare_datasets.py /data/old /data/new /output/qc

  # Text report to stdout
  python compare_datasets.py /data/old /data/new /output/qc --text-only

  # Text report to file
  python compare_datasets.py /data/old /data/new /output/qc --text-only -o qc_report.txt

  # Compare only specific variables with visualizations
  python compare_datasets.py /data/old /data/new /output/qc tasmax,pr,hursmin

  # Compare 20 randomly selected files per variable with text output
  python compare_datasets.py /data/old /data/new /output/qc tasmax,pr 20 --text-only

  # Compare 10 files of all variables with visualizations
  python compare_datasets.py /data/old /data/new /output/qc "" 10
        """,
    )

    parser.add_argument("old_dir", help="Path to reference/old dataset directory")
    parser.add_argument("new_dir", help="Path to new dataset directory to validate")
    parser.add_argument(
        "output_dir", help="Directory where QC outputs (notebooks/html) will be saved"
    )
    parser.add_argument(
        "variables",
        nargs="?",
        default="",
        help='Comma-separated list of variables to check (e.g., "tasmax,pr"). Leave empty for all variables.',
    )
    parser.add_argument(
        "max_files",
        nargs="?",
        type=int,
        default=None,
        help="Maximum number of files to randomly sample per variable. Leave empty to check all files.",
    )
    parser.add_argument(
        "--text-only",
        action="store_true",
        help="Generate text report only (skip visualization notebooks). Default: generate visualizations.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default=None,
        help="Path to text output file. If not specified with --text-only, prints to stdout.",
    )
    parser.add_argument(
        "--shapefile",
        type=str,
        default=None,
        help="Path to shapefile for boundary overlay in visualizations",
    )

    # If no arguments provided, print help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # Parse variables list
    variables = None
    if args.variables and args.variables.strip():
        variables = [v.strip() for v in args.variables.split(",") if v.strip()]

    # Run QC
    exit_code = run_qc(
        args.old_dir,
        args.new_dir,
        args.output_dir,
        variables,
        args.max_files,
        text_only=args.text_only,
        text_output=args.output,
        shapefile_path=args.shapefile,
    )
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
