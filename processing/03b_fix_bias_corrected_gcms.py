"""
Fix Bias-Corrected GCM Data - Apply CLIP_HURSMIN Correction

This script reads existing bias-corrected GCM data and applies the CLIP_HURSMIN
routine to clamp hursmin values to the valid range [0.0, 100.0]. This is useful
when working with bias-corrected data from the original dataset that was not
processed with the CLIP_HURSMIN correction.

The script reads from CMIP6_BIAS_CORRECTED (if set) or OUT_DIR/bias_corrected,
applies the correction, and writes the fixed files to OUT_DIR/bias_corrected.

Usage:
    Set environment variables:
        - CMIP6_BIAS_CORRECTED: Path to original bias-corrected data (optional)
        - OUT_DIR: Output directory for corrected files
        - CLIP_HURSMIN: Set to TRUE to apply clipping (default: TRUE)

    Run: python 03b_fix_bias_corrected_gcms.py
"""

from tqdm import tqdm
import xarray as xr
from pathlib import Path
from datetime import datetime
import time

from config import (
    OUT_DIR,
    CMIP6_BIAS_CORRECTED,
    CLIP_HURSMIN,
    gcm_list,
    hist_years,
    sim_periods,
)

import warnings

warnings.filterwarnings("ignore")


def clip_hursmin(data_array: xr.DataArray) -> xr.DataArray:
    """
    Clip hursmin values to valid range [0.0, 100.0].
    Values below 0 are set to 0, values above 100 are set to 100.

    Parameters
    ----------
    data_array: xarray.DataArray
        The data array to clip

    Returns
    -------
    xarray.DataArray with clipped values
    """
    if CLIP_HURSMIN:
        return data_array.clip(min=0.0, max=100.0)
    return data_array


def get_all_years():
    """Get all years that need processing (historical + all simulation periods)"""
    all_years = []
    # Historical years
    all_years.extend(range(hist_years[0], hist_years[1] + 1))
    # Simulation years
    for period in sim_periods:
        all_years.extend(range(period[0], period[1] + 1))
    return sorted(all_years)


if __name__ == "__main__":

    # Track script execution time
    start_time = time.time()
    start_datetime = datetime.now()
    print(f"\n{'='*60}")
    print(f"Fix Bias-Corrected GCMs Script Started")
    print(f"Start time: {start_datetime.strftime('%Y-%m-%d %H:%M:%S')}")

    if not CLIP_HURSMIN:
        print(f"⚠️  CLIP_HURSMIN DISABLED - No correction will be applied")
        print(f"   (Set CLIP_HURSMIN=TRUE to apply hursmin clipping)")

    # Determine input directory
    if CMIP6_BIAS_CORRECTED:
        input_base = Path(CMIP6_BIAS_CORRECTED)
        print(f"Input directory: {input_base}")
    else:
        input_base = Path(OUT_DIR).joinpath("bias_corrected")
        print(f"Input directory: {input_base} (using OUT_DIR/bias_corrected)")

    output_base = Path(OUT_DIR).joinpath("bias_corrected")
    print(f"Output directory: {output_base}")
    print(f"{'='*60}\n")

    # Get all years to process
    all_years = get_all_years()

    # Only process hursmin variable
    var = "hursmin"

    # Calculate total number of files to process
    N = len(gcm_list) * len(all_years)

    with tqdm(total=N, desc="Processing hursmin files") as pbar:
        for gcm in gcm_list:
            # Input and output directories for this GCM
            in_dir = input_base.joinpath(gcm)
            out_dir = output_base.joinpath(gcm)

            # Create output directory if it doesn't exist
            if not out_dir.exists():
                out_dir.mkdir(parents=True)

            for yr in all_years:
                pbar.set_postfix({"GCM": gcm, "year": yr})

                # Construct filename
                filename = f"{var}_{gcm}_{yr}.nc"
                in_file = in_dir.joinpath(filename)
                out_file = out_dir.joinpath(filename)

                # Check if input file exists
                if not in_file.exists():
                    print(f"\n⚠️  Warning: File not found: {in_file}")
                    pbar.update()
                    continue

                # Check if output file already exists (skip if same as input)
                if out_file.exists() and in_file == out_file:
                    # Skip if we're reading and writing to the same location
                    pbar.update()
                    continue

                try:
                    # Load the data
                    ds = xr.open_dataset(in_file, engine="h5netcdf")

                    # Get the hursmin data array
                    hursmin_data = ds[var]

                    # Apply clipping
                    hursmin_clipped = clip_hursmin(hursmin_data)

                    # Create new dataset with clipped data
                    ds_clipped = hursmin_clipped.to_dataset(name=var)

                    # Preserve attributes
                    ds_clipped.attrs = ds.attrs
                    ds_clipped[var].attrs = ds[var].attrs

                    # Convert to float32 for storage efficiency
                    ds_clipped = ds_clipped.astype("float32")

                    # Use float64 for time to avoid precision loss with timestamps
                    encoding = {
                        var: {"dtype": "float32"},
                        "time": {"dtype": "float64"},
                    }

                    # Write to output file
                    ds_clipped.to_netcdf(out_file, engine="h5netcdf", encoding=encoding)

                    # Close datasets
                    ds.close()
                    ds_clipped.close()

                except Exception as e:
                    print(f"\n❌ Error processing {filename}: {str(e)}")

                pbar.update()

    # Calculate and display elapsed time
    end_time = time.time()
    end_datetime = datetime.now()
    elapsed_seconds = end_time - start_time

    hours = int(elapsed_seconds // 3600)
    minutes = int((elapsed_seconds % 3600) // 60)
    seconds = int(elapsed_seconds % 60)

    print(f"\n{'='*60}")
    print(f"Fix Bias-Corrected GCMs Script Completed!")
    print(f"End time: {end_datetime.strftime('%Y-%m-%d %H:%M:%S')}")
    if hours > 0:
        print(f"Total elapsed time: {hours}h {minutes}m {seconds}s")
    elif minutes > 0:
        print(f"Total elapsed time: {minutes}m {seconds}s")
    else:
        print(f"Total elapsed time: {seconds}s")
    print(f"{'='*60}\n")
