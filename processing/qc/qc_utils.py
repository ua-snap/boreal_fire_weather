"""
Quality Control Utilities for Dataset Comparison

This module provides functions for comparing datasets and visualizing differences.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
import pandas as pd
from pathlib import Path
import xarray as xr
import geopandas as gpd


def format_time_value(t):
    """
    Convert time value to string, handling both standard datetime and cftime objects.

    Parameters
    ----------
    t : datetime-like
        Time value to format

    Returns
    -------
    str
        Formatted date string in YYYY-MM-DD format
    """
    # Check if it's a cftime object (has strftime method but not a pandas Timestamp)
    if hasattr(t, "strftime") and not isinstance(t, pd.Timestamp):
        return t.strftime("%Y-%m-%d")
    else:
        # Try pandas conversion for standard datetime objects
        try:
            return pd.to_datetime(t).strftime("%Y-%m-%d")
        except:
            # Fallback: just convert to string
            return str(t)


def compare_timestep(ds_ref, ds_new, var, time_idx):
    """
    Compare a single timestep between two datasets for a given variable.

    Parameters
    ----------
    ds_ref : xarray.Dataset
        Reference dataset
    ds_new : xarray.Dataset
        New dataset to compare
    var : str
        Variable name to compare
    time_idx : int
        Time index to compare

    Returns
    -------
    dict
        Dictionary containing comparison metrics:
        - n_nan_changed: Number of cells where NaN pattern differs
        - n_value_diff: Number of non-NaN cells where values differ
        - max_abs_diff: Maximum absolute difference in values
        - mean_ref: Mean of reference data
        - std_ref: Standard deviation of reference data
        - mean_new: Mean of new data
        - std_new: Standard deviation of new data
    """
    data_ref = ds_ref[var].isel(time=time_idx).values
    data_new = ds_new[var].isel(time=time_idx).values

    # NaN pattern analysis
    nan_ref = np.isnan(data_ref)
    nan_new = np.isnan(data_new)
    n_nan_changed = np.sum((nan_ref & ~nan_new) | (~nan_ref & nan_new))

    # Value comparison
    valid_mask = ~nan_ref & ~nan_new
    n_value_diff = (
        np.sum(data_ref[valid_mask] != data_new[valid_mask])
        if np.any(valid_mask)
        else 0
    )

    # Calculate difference statistics
    diff = data_new - data_ref
    max_abs_diff = np.nanmax(np.abs(diff)) if np.any(~np.isnan(diff)) else 0.0

    return {
        "n_nan_changed": n_nan_changed,
        "n_value_diff": n_value_diff,
        "max_abs_diff": max_abs_diff,
        "mean_ref": np.nanmean(data_ref),
        "std_ref": np.nanstd(data_ref),
        "mean_new": np.nanmean(data_new),
        "std_new": np.nanstd(data_new),
    }


def compare_all_timesteps(ds_ref, ds_new, var):
    """
    Compare all timesteps between two datasets for a given variable.

    Parameters
    ----------
    ds_ref : xarray.Dataset
        Reference dataset
    ds_new : xarray.Dataset
        New dataset to compare
    var : str
        Variable name to compare

    Returns
    -------
    pandas.DataFrame
        DataFrame with comparison metrics for each timestep
    """
    n_times = len(ds_ref.time)
    results = []

    for t_idx in range(n_times):
        time_val = ds_ref.time.values[t_idx]
        metrics = compare_timestep(ds_ref, ds_new, var, t_idx)
        metrics["time_idx"] = t_idx
        metrics["time"] = time_val
        results.append(metrics)

    df = pd.DataFrame(results)
    return df


def create_comparison_plot(
    ds_ref, ds_new, var, time_idx, gdf=None, bounds=None, title_suffix=""
):
    """
    Create a 4-panel comparison plot for a single timestep.

    Parameters
    ----------
    ds_ref : xarray.Dataset
        Reference dataset
    ds_new : xarray.Dataset
        New dataset to compare
    var : str
        Variable name to plot
    time_idx : int
        Time index to plot
    gdf : geopandas.GeoDataFrame, optional
        Shapefile boundary to overlay
    bounds : list, optional
        [minx, miny, maxx, maxy] for plot extent
    title_suffix : str, optional
        Additional text to append to the title

    Returns
    -------
    matplotlib.figure.Figure
        The created figure
    """
    # Extract data
    data_ref = ds_ref[var].isel(time=time_idx).values
    data_new = ds_new[var].isel(time=time_idx).values
    time_val = ds_ref.time.values[time_idx]

    # Calculate difference
    diff = data_new - data_ref

    # Create NaN pattern masks
    nan_ref = np.isnan(data_ref)
    nan_new = np.isnan(data_new)
    nan_pattern = np.zeros_like(data_ref)
    nan_pattern[nan_ref & ~nan_new] = 1  # NaN in ref only
    nan_pattern[~nan_ref & nan_new] = 2  # NaN in new only
    nan_pattern[nan_ref & nan_new] = 3  # NaN in both

    # Count differences
    n_nan_changed = np.sum((nan_ref & ~nan_new) | (~nan_ref & nan_new))
    valid_mask = ~nan_ref & ~nan_new
    n_value_diff = (
        np.sum(data_ref[valid_mask] != data_new[valid_mask])
        if np.any(valid_mask)
        else 0
    )

    # Get lat/lon for plotting
    lon = ds_ref.lon.values
    lat = ds_ref.lat.values

    # Set bounds if not provided
    if bounds is None:
        bounds = [lon.min(), lat.min(), lon.max(), lat.max()]

    # Create figure with 4 subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(
        f"{var.upper()} Comparison - {time_val}{title_suffix}\n"
        f"NaN changes: {n_nan_changed}, Value differences: {n_value_diff}",
        fontsize=14,
        fontweight="bold",
    )

    # Determine color limits using percentiles (excluding NaNs)
    valid_data = np.concatenate(
        [data_ref[~np.isnan(data_ref)], data_new[~np.isnan(data_new)]]
    )
    if len(valid_data) > 0:
        vmin = np.percentile(valid_data, 2)
        vmax = np.percentile(valid_data, 98)
    else:
        vmin, vmax = 0, 1

    # Plot 1: Original data
    ax = axes[0, 0]
    im1 = ax.pcolormesh(lon, lat, data_ref, cmap="YlOrRd", vmin=vmin, vmax=vmax)
    if gdf is not None:
        gdf.boundary.plot(ax=ax, edgecolor="black", linewidth=1)
    ax.set_xlim(bounds[0], bounds[2])
    ax.set_ylim(bounds[1], bounds[3])
    ax.set_title("Original Data")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    plt.colorbar(im1, ax=ax, label=var)

    # Plot 2: New data
    ax = axes[0, 1]
    im2 = ax.pcolormesh(lon, lat, data_new, cmap="YlOrRd", vmin=vmin, vmax=vmax)
    if gdf is not None:
        gdf.boundary.plot(ax=ax, edgecolor="black", linewidth=1)
    ax.set_xlim(bounds[0], bounds[2])
    ax.set_ylim(bounds[1], bounds[3])
    ax.set_title("New Data")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    plt.colorbar(im2, ax=ax, label=var)

    # Plot 3: Difference (new - original)
    ax = axes[1, 0]
    if np.any(~np.isnan(diff)):
        diff_abs_max = np.nanpercentile(np.abs(diff), 98)
        if diff_abs_max == 0:
            diff_abs_max = 1
        norm = TwoSlopeNorm(vmin=-diff_abs_max, vcenter=0, vmax=diff_abs_max)
        im3 = ax.pcolormesh(lon, lat, diff, cmap="RdBu_r", norm=norm)
    else:
        im3 = ax.pcolormesh(lon, lat, diff, cmap="RdBu_r")
    if gdf is not None:
        gdf.boundary.plot(ax=ax, edgecolor="black", linewidth=1)
    ax.set_xlim(bounds[0], bounds[2])
    ax.set_ylim(bounds[1], bounds[3])
    ax.set_title("Difference (New - Original)")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    plt.colorbar(im3, ax=ax, label=f"Δ{var}")

    # Plot 4: NaN pattern changes
    ax = axes[1, 1]
    colors = ["white", "red", "blue", "gray"]
    labels = ["Valid in both", "NaN in original only", "NaN in new only", "NaN in both"]
    im4 = ax.pcolormesh(
        lon, lat, nan_pattern, cmap=plt.cm.colors.ListedColormap(colors), vmin=0, vmax=3
    )
    if gdf is not None:
        gdf.boundary.plot(ax=ax, edgecolor="black", linewidth=1)
    ax.set_xlim(bounds[0], bounds[2])
    ax.set_ylim(bounds[1], bounds[3])
    ax.set_title("NaN Pattern Changes")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

    # Create custom legend
    patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(4)]
    ax.legend(handles=patches, loc="upper right", fontsize=8)

    plt.tight_layout()

    return fig


def create_summary_table(df_differences):
    """
    Create a formatted summary table from difference DataFrame.

    Parameters
    ----------
    df_differences : pandas.DataFrame
        DataFrame containing timesteps with differences

    Returns
    -------
    pandas.DataFrame
        Formatted table ready for display
    """
    if len(df_differences) == 0:
        return pd.DataFrame({"Message": ["No differences found across all timesteps"]})

    # Select and format columns
    table = df_differences[
        [
            "time",
            "mean_ref",
            "std_ref",
            "mean_new",
            "std_new",
            "max_abs_diff",
            "n_value_diff",
            "n_nan_changed",
        ]
    ].copy()

    # Format time column - handle both standard datetime and cftime objects
    table["time"] = table["time"].apply(format_time_value)

    # Round numerical columns
    for col in ["mean_ref", "std_ref", "mean_new", "std_new", "max_abs_diff"]:
        table[col] = table[col].round(3)

    # Rename columns for display
    table.columns = [
        "Date",
        "Mean (Ref)",
        "Std (Ref)",
        "Mean (New)",
        "Std (New)",
        "Max Abs Diff",
        "Value Diffs",
        "NaN Changes",
    ]

    return table


def generate_qc_summary(comparison_results, variables):
    """
    Generate a markdown-formatted QC summary.

    Parameters
    ----------
    comparison_results : dict
        Dictionary mapping variable names to comparison DataFrames
    variables : list
        List of variable names that were compared

    Returns
    -------
    str
        Markdown-formatted summary text
    """
    summary_lines = ["## QC Summary\n"]

    for var in variables:
        if var not in comparison_results:
            summary_lines.append(f"### {var.upper()}")
            summary_lines.append("- ⚠️ Variable not found in one or both datasets\n")
            continue

        df = comparison_results[var]
        df_diff = df[(df["n_value_diff"] > 0) | (df["n_nan_changed"] > 0)]

        summary_lines.append(f"### {var.upper()}")
        summary_lines.append(f"- **Total timesteps analyzed:** {len(df)}")
        summary_lines.append(f"- **Timesteps with differences:** {len(df_diff)}")

        if len(df_diff) == 0:
            summary_lines.append("- ✅ **All timesteps are identical**\n")
        else:
            summary_lines.append(f"- ⚠️ **Differences detected**")

            # Max value difference info
            max_val_idx = df_diff["max_abs_diff"].idxmax()
            max_val_row = df.loc[max_val_idx]
            max_val_time = format_time_value(max_val_row["time"])
            summary_lines.append(
                f"  - Maximum value difference: {max_val_row['max_abs_diff']:.6f} on {max_val_time}"
            )

            # Max NaN difference info
            max_nan_idx = df_diff["n_nan_changed"].idxmax()
            max_nan_row = df.loc[max_nan_idx]
            max_nan_time = format_time_value(max_nan_row["time"])
            summary_lines.append(
                f"  - Maximum NaN changes: {int(max_nan_row['n_nan_changed'])} cells on {max_nan_time}"
            )

            # Overall value diff count
            total_value_diffs = df_diff["n_value_diff"].sum()
            total_nan_changes = df_diff["n_nan_changed"].sum()
            summary_lines.append(
                f"  - Total value differences: {int(total_value_diffs)} across all timesteps"
            )
            summary_lines.append(
                f"  - Total NaN pattern changes: {int(total_nan_changes)} across all timesteps\n"
            )

    return "\n".join(summary_lines)


# High-level wrapper functions for notebook simplicity


def load_datasets_and_shapefile(ref_file, new_file, shapefile_path):
    """
    Load datasets and shapefile with all necessary setup.

    Parameters
    ----------
    ref_file : str
        Path to reference dataset
    new_file : str
        Path to new dataset
    shapefile_path : str
        Path to shapefile for boundary overlay

    Returns
    -------
    tuple
        (ds_ref, ds_new, gdf, bounds) where gdf may be None if shapefile not found
    """
    print("Loading datasets...")
    ds_ref = xr.open_dataset(ref_file)
    ds_new = xr.open_dataset(new_file)

    # Load shapefile for boreal region boundary
    if Path(shapefile_path).exists():
        gdf = gpd.read_file(shapefile_path)
        bounds = gdf.total_bounds  # [minx, miny, maxx, maxy]
        print(f"Shapefile bounds: {bounds}")
    else:
        print(f"Warning: Shapefile not found at {shapefile_path}")
        gdf = None
        # Use data extent as fallback
        bounds = [
            ds_ref.lon.min().values,
            ds_ref.lat.min().values,
            ds_ref.lon.max().values,
            ds_ref.lat.max().values,
        ]

    print(f"Dataset reference shape: {ds_ref.dims}")
    print(f"Dataset new shape: {ds_new.dims}")
    print(f"Number of timesteps: {len(ds_ref.time)}")

    return ds_ref, ds_new, gdf, bounds


def run_comparison_analysis(ds_ref, ds_new, variables):
    """
    Run comparison analysis for all variables across all timesteps.

    Parameters
    ----------
    ds_ref : xarray.Dataset
        Reference dataset
    ds_new : xarray.Dataset
        New dataset to compare
    variables : list
        List of variable names to compare

    Returns
    -------
    dict
        Dictionary mapping variable names to comparison DataFrames
    """
    print("Comparing all timesteps...")
    comparison_results = {}

    for var in variables:
        if var not in ds_ref.data_vars or var not in ds_new.data_vars:
            print(f"Skipping {var}: not found in one or both datasets")
            continue

        print(f"\nAnalyzing {var.upper()}...")
        df_comparison = compare_all_timesteps(ds_ref, ds_new, var)
        comparison_results[var] = df_comparison

        # Filter for timesteps with differences
        df_diff = df_comparison[
            (df_comparison["n_value_diff"] > 0) | (df_comparison["n_nan_changed"] > 0)
        ]

        print(f"  Total timesteps: {len(df_comparison)}")
        print(f"  Timesteps with differences: {len(df_diff)}")

        if len(df_diff) > 0:
            print(f"  Max absolute difference: {df_diff['max_abs_diff'].max():.6f}")
            print(f"  Max NaN changes: {int(df_diff['n_nan_changed'].max())}")

    print("\nComparison complete!")
    return comparison_results


def display_summary_tables(comparison_results, variables):
    """
    Display summary statistics tables for all variables.

    Parameters
    ----------
    comparison_results : dict
        Dictionary mapping variable names to comparison DataFrames
    variables : list
        List of variable names to display
    """
    for var in variables:
        if var not in comparison_results:
            continue

        df_comparison = comparison_results[var]
        df_diff = df_comparison[
            (df_comparison["n_value_diff"] > 0) | (df_comparison["n_nan_changed"] > 0)
        ]

        print(f"\n{'='*80}")
        print(f"{var.upper()} - Timesteps with Differences")
        print("=" * 80)

        if len(df_diff) == 0:
            print("✅ All timesteps are identical!")
        else:
            table = create_summary_table(df_diff)
            print(table.to_string(index=False))


def plot_max_value_differences(
    ds_ref, ds_new, comparison_results, variables, gdf, bounds
):
    """
    Plot timesteps with maximum value differences for all variables.

    Parameters
    ----------
    ds_ref : xarray.Dataset
        Reference dataset
    ds_new : xarray.Dataset
        New dataset to compare
    comparison_results : dict
        Dictionary mapping variable names to comparison DataFrames
    variables : list
        List of variable names to plot
    gdf : geopandas.GeoDataFrame or None
        Shapefile boundary to overlay
    bounds : list
        [minx, miny, maxx, maxy] for plot extent
    """
    for var in variables:
        if var not in comparison_results:
            continue

        df_comparison = comparison_results[var]
        df_diff = df_comparison[
            (df_comparison["n_value_diff"] > 0) | (df_comparison["n_nan_changed"] > 0)
        ]

        if len(df_diff) == 0:
            print(f"\n{var.upper()}: No differences to plot")
            continue

        # Find timestep with maximum value difference
        max_val_idx = df_diff["max_abs_diff"].idxmax()
        max_val_row = df_comparison.loc[max_val_idx]
        time_idx = int(max_val_row["time_idx"])

        print(f"\n{var.upper()}: Plotting timestep with maximum value difference")
        print(f"  Time: {max_val_row['time']}")
        print(f"  Max absolute difference: {max_val_row['max_abs_diff']:.6f}")

        fig = create_comparison_plot(
            ds_ref,
            ds_new,
            var,
            time_idx,
            gdf,
            bounds,
            title_suffix=" (Max Value Difference)",
        )
        fig.show()


def plot_max_nan_differences(
    ds_ref, ds_new, comparison_results, variables, gdf, bounds
):
    """
    Plot timesteps with maximum NaN pattern changes for all variables.

    Parameters
    ----------
    ds_ref : xarray.Dataset
        Reference dataset
    ds_new : xarray.Dataset
        New dataset to compare
    comparison_results : dict
        Dictionary mapping variable names to comparison DataFrames
    variables : list
        List of variable names to plot
    gdf : geopandas.GeoDataFrame or None
        Shapefile boundary to overlay
    bounds : list
        [minx, miny, maxx, maxy] for plot extent
    """
    for var in variables:
        if var not in comparison_results:
            continue

        df_comparison = comparison_results[var]
        df_diff = df_comparison[
            (df_comparison["n_value_diff"] > 0) | (df_comparison["n_nan_changed"] > 0)
        ]

        if len(df_diff) == 0:
            continue

        # Find timestep with maximum NaN changes
        max_nan_idx = df_diff["n_nan_changed"].idxmax()
        max_nan_row = df_comparison.loc[max_nan_idx]
        time_idx = int(max_nan_row["time_idx"])

        # Only plot if different from max value difference timestep
        max_val_idx = df_diff["max_abs_diff"].idxmax()
        if max_nan_idx == max_val_idx:
            print(
                f"\n{var.upper()}: Max NaN changes occur at same timestep as max value difference"
            )
            continue

        print(f"\n{var.upper()}: Plotting timestep with maximum NaN changes")
        print(f"  Time: {max_nan_row['time']}")
        print(f"  NaN changes: {int(max_nan_row['n_nan_changed'])}")

        fig = create_comparison_plot(
            ds_ref,
            ds_new,
            var,
            time_idx,
            gdf,
            bounds,
            title_suffix=" (Max NaN Changes)",
        )
        fig.show()
