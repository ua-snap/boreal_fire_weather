from tqdm import tqdm
import dask
from dask.diagnostics import ProgressBar
import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
import time
from xclim.core.units import str2pint
from xclim.sdba import QuantileDeltaMapping
from xclim.sdba.processing import jitter_under_thresh

from config import ERA5_PROCESSED, CMIP6_PROCESSED, OUT_DIR, SHP_MASK, gcm_list, metvars, hist_years, sim_periods, qdm
from utils import *

import warnings

# Filter out warning on all-nan slice operations, expected
warnings.filterwarnings(
    "ignore", message="All-NaN slice encountered in interp_on_quantiles"
)

# Suppress dask warnings on chunk size
dask.config.set({"array.slicing.split_large_chunks": False})


def same_vals(x) -> bool:
    """
    Check if all values in an interable object are equal to one another.
    """
    return x.count(x[0]) == len(x)


def check_var_names(*args):
    """
    Make sure all variables are the same and not mixing two different variables
    together.
    """
    var_names = [get_var_names(x)[0] for x in args]
    if not same_vals(var_names):
        raise Exception("Data variables need to be the same in all datasets.")
    return var_names[0]


def same_units(*args):
    """
    Check to make sure units are the same in each dataset. Raises exception
    if not the case.
    """
    units = [str2pint(x.units).units for x in args]
    if not same_vals(units):
        raise Exception("Units need to be the same in all datasets.")
    return None


def dimcheck_and_regrid(
    ref: xr.DataArray, hst: xr.DataArray, sim: xr.DataArray, regrid: str = "era2gcm"
) -> tuple:
    """
    Description
    -----------
    Make sure the grid sizes are the same for reference (ref), historical
    (hst), and simulation (sim) datasets.

    Parameters
    ----------
    ref: xarray.DataArray
        Reference dataset
    hst: xarray.DataArray
        Historical GCM output
    sim: xarray.DataArray
        Projected GCM Output
    regrid: str
        Which direction to do regridding/linear interpolation. Options are
        'era2gcm' which will regrid ERA5 data to scale of GCM. Other option is
        'gcm2era' which will regrid GCM to scale of ERA5.

    Returns
    -------
    tuple of xarray.DatasArrays for ref, hst, and sim
    """
    ref_shape = ref.shape
    hst_shape = hst.shape
    sim_shape = sim.shape

    shape_check_1 = same_vals([ref_shape, hst_shape, sim_shape])

    # If datasets are not the same size then regridding needs to be done
    if not shape_check_1:  # If datasets are not the same size ...
        if not any([regrid == "era2gcm", regrid == "gcm2era"]):
            raise Exception(
                "Array shapes not equal, need to do regridding and \
                'regrid' argument needs to be either 'era2gcm' or 'gcm2era'"
            )
        if regrid == "era2gcm":
            ref = regrid_geodata(ref, hst)
        elif regrid == "gcm2era":
            hst = regrid_geodata(hst, ref)
            sim = regrid_geodata(sim, ref)

        ref_shape = ref.shape
        hst_shape = ref.shape
        sim_shape = sim.shape

        # Check to make sure regridding worked and also double checks to make
        # sure time dims are the same size across datasets
        if not same_vals([ref_shape, hst_shape, sim_shape]):
            raise Exception(
                "ref,hst,sim array shapes are not equal! \
                Double check input. Time dimensions need to be same as well."
            )
    return (ref, hst, sim)


def mask_arrays(
    ref: xr.DataArray, hst: xr.DataArray, sim: xr.DataArray, mask=None
) -> tuple:
    """
    Description
    -----------
    Mask out spatial areas not of interest in ref, hst, and sim
    xarray.DataArrays.

    Parameters
    ----------
    ref: xarray.DataArray
        Reference dataset
    hst: xarray.DataArray
        Historical GCM output
    sim: xarray.DataArray
        Projected GCM Output
    mask: None or str or pathlib.Path
        If a spatial mask is to be applied, provide path to ESRI shapefile
        that defines the mask

    Returns
    -------
    tuple of xarray.DatasArrays for ref, hst, and sim
    """
    if mask is not None:
        if isinstance(mask, np.ndarray):
            mask_grd = mask.copy()
            del mask
        elif isinstance(mask, pathlib.PosixPath) or isinstance(mask, str):
            mask_grd = mask_from_shp(mask, ref)
        mask_3d = np.broadcast_to(mask_grd, ref.shape)

        ref = ref.where(mask_3d)
        hst = hst.where(mask_3d)
        sim = sim.where(mask_3d)
    return (ref, hst, sim)


def chunk_data_arrays(ref, hst, sim, frac=0.2):
    """
    Partition dataset into dask chunks for parallel processing if set to True.
    Currently set to 20% the size of each dimension, can modify this later to
    make it more dynamic. Time dimension cannot be broken up for doing quantile
    mapping.
    """
    axes = get_geoaxes(ref)
    n, m = (ref[axes["Y"]].size, ref[axes["X"]].size)

    ref = ref.chunk(chunks={"time": -1, "lat": round(n * frac), "lon": round(m * frac)})

    hst = hst.chunk(chunks={"time": -1, "lat": round(n * frac), "lon": round(m * frac)})

    sim = sim.chunk(chunks={"time": -1, "lat": round(n * frac), "lon": round(m * frac)})
    return (ref, hst, sim)


def quantile_delta_mapping(
    ref_src: list,
    hst_src: list,
    sim_src: list,
    min_thresh: float = None,
    return_hst=False,
    dask_load=False,
    dask_return=False,
    **kwargs
) -> tuple:
    """
    Description
    -----------
    Perform quantile delta mapping usin xclim package

    Parameters
    ----------
    ref_src: list
        List of paths to files of reference datasets
    hst_src: list
        List of paths to files of historical datasets
    sim_src: list
        List of paths to files of simulation/projected datasets
    min_thresh: float
        Minimum threshold below which all values are assumed equal to zero
    return_hst: bool
        Should the bias-corrected historical data also be returned
    dask_load: bool
        Should the dataset be read in using dask in parallel
    dask_return: bool
        Should a dask array be returned
    **kwargs: additional keyword arguments to be passed on to various functions

    Returns
    -------
    tuple of xarray.DataArrays containing resulsts from bias correcting
    """
    if dask_load is False:
        dask_return = False

    # Read in datasets for reference time period (ref), historical overlap of
    # gcm (hst), and future simulation/projection period (sim)
    ref = xr.open_mfdataset(
        ref_src,
        engine="h5netcdf",
        parallel=dask_load,
        chunks={"time": 365, "lon": -1, "lat": -1},
    )

    hst = xr.open_mfdataset(
        hst_src,
        engine="h5netcdf",
        parallel=dask_load,
        chunks={"time": 365, "lon": -1, "lat": -1},
    )

    sim = xr.open_mfdataset(
        sim_src,
        engine="h5netcdf",
        parallel=dask_load,
        chunks={"time": 365, "lon": -1, "lat": -1},
    )

    # Check to make sure the variable (e.g. precip) is the same for each dataset
    var = check_var_names(ref, hst, sim)

    # Convert from dataset to dataarray
    ref = ref[var]
    hst = hst[var]
    sim = sim[var]

    # Will raise Exception if not all the same.
    same_units(ref, hst, sim)

    # Align time coordinates - ensure all datasets use consistent time encoding
    # This handles cases where time units have different offsets (e.g., midnight vs noon)
    # Reindex hst to match ref's exact time coordinates (same dates)
    if not ref.time.equals(hst.time):
        hst = hst.reindex(time=ref.time, method="nearest", tolerance=pd.Timedelta("12h"))
    
    # For sim (which may be a different time period), normalize time-of-day to match ref
    # by flooring all times to the start of the day
    ref = ref.assign_coords(time=ref.time.dt.floor("D"))
    hst = hst.assign_coords(time=hst.time.dt.floor("D"))
    sim = sim.assign_coords(time=sim.time.dt.floor("D"))

    # Do regridding so all data arrays are aligned and have the same shape and
    # dimensions
    ref, hst, sim = dimcheck_and_regrid(
        ref, hst, sim, **get_kwargs(("regrid",), kwargs)
    )

    # Apply spatial mask, does nothing if 'mask=None'
    ref, hst, sim = mask_arrays(ref, hst, sim, **get_kwargs(("mask",), kwargs))

    # If working with dask arrays, set chunks so time dimension is not broken up
    # This is a requirement for using sdba.QuantileDeltaMapping
    if dask_load:
        ref, hst, sim = chunk_data_arrays(
            ref, hst, sim, **get_kwargs(("frac",), kwargs)
        )

    # Apply uniform random variable if value is less then preset amount. Mainly
    # used for precip (e.g., 0.5 mm/day)
    if min_thresh is not None:

        thresh_str = "%0.3f %s" % (min_thresh, ref.attrs["units"])

        ref = jitter_under_thresh(ref, thresh_str)
        hst = jitter_under_thresh(hst, thresh_str)
        sim = jitter_under_thresh(sim, thresh_str)

    # Training step in Quantile delta mapping
    QDM = QuantileDeltaMapping.train(
        ref, hst, **get_kwargs(("group", "nquantiles", "kind"), kwargs)
    )

    # Declare empty tuple to store output
    return_ds = ()

    # Only do historical bias corrections if return_hist=True
    if return_hst:

        # Bias correcting step for hst time period
        hst_ba = QDM.adjust(hst, **get_kwargs(("interp", "extrapolation"), kwargs))

        # Set all jittered values back to 0
        if min_thresh is not None:

            hst_ba = xr.where(hst_ba > min_thresh, hst_ba, 0.0)

        hst_ba = hst_ba.rename(var)
        hst_ba = hst_ba.assign_attrs(hst.attrs)

        return_ds = return_ds + tuple([hst_ba])

    # Bias correcting step for projected/simulated time period
    sim_ba = QDM.adjust(sim, **get_kwargs(("interp", "extrapolation"), kwargs))

    # Set all jittered values back to 0
    if min_thresh is not None:

        sim_ba = xr.where(sim_ba > min_thresh, sim_ba, 0.0)

    # The xclim function renames the variable, this is just setting it back to
    # the original name.
    sim_ba = sim_ba.rename(var)

    # Makes sure the proper attributes remain witht the data array after
    # using the xclim function
    sim_ba = sim_ba.assign_attrs(sim.attrs)

    # Add bias adjusted data array to the export tuple
    return_ds = return_ds + tuple([sim_ba])

    # Export results
    if (not dask_return) & (dask.is_dask_collection(return_ds[0])):
        return tuple([x.compute() for x in return_ds])
    else:
        return return_ds


if __name__ == "__main__":

    # Track script execution time
    start_time = time.time()
    start_datetime = datetime.now()
    print(f"\n{'='*60}")
    print(f"Bias Correction Script Started")
    print(f"Start time: {start_datetime.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*60}\n")

    # check if SHP_MASK is defined
    if SHP_MASK is None:
        shpfile = None
    else:
        shpfile = Path(SHP_MASK)

    era5_dir = Path(ERA5_PROCESSED)

    oper = qdm["oper"]
    min_thresh = qdm["min_thresh"]
    quantile_vals = np.array(qdm["quantile_vals"], dtype=np.float32)

    N = len(gcm_list) * len(metvars) * len(sim_periods)
    with tqdm(total=N) as pbar:
        for gcm in gcm_list:
            in_dir = Path(CMIP6_PROCESSED).joinpath(gcm)
            out_dir = Path(OUT_DIR).joinpath("bias_corrected", gcm)
            if out_dir.exists() is False:
                out_dir.mkdir(parents=True)
            for var in metvars:
                for sim_year in sim_periods:
                    ref_src = [
                        list(era5_dir.glob("%s*%d*" % (var, i)))[0]
                        for i in range(hist_years[0], hist_years[1] + 1)
                    ]
                    hst_src = [
                        list(in_dir.glob("%s*%d*" % (var, i)))[0]
                        for i in range(hist_years[0], hist_years[1] + 1)
                    ]
                    sim_src = [
                        list(in_dir.glob("%s*%d*" % (var, i)))[0]
                        for i in range(sim_year[0], sim_year[1] + 1)
                    ]

                    if sim_year == sim_periods[0]:
                        hst_return_bool = True
                    else:
                        hst_return_bool = False

                    qdm_arrays = quantile_delta_mapping(
                        ref_src,
                        hst_src,
                        sim_src,
                        dask_load=True,
                        regrid="gcm2era",
                        mask=shpfile,
                        kind=oper[var],
                        nquantiles=quantile_vals,
                        group="time.month",
                        min_thresh=min_thresh[var],
                        return_hst=hst_return_bool,
                        interp="linear",
                        dask_return=True,
                    )

                    if hst_return_bool:

                        hst_ba = qdm_arrays[0]
                        with ProgressBar():
                            hst_ba = hst_ba.compute()

                        for yr in range(hist_years[0], hist_years[1] + 1):
                            fn = out_dir.joinpath("%s_%s_%d.nc" % (var, gcm, yr))
                            yr_slice = slice(str(yr), str(yr))
                            export_ds = hst_ba.sel(time=yr_slice)
                            export_ds = export_ds.astype("float32")
                            # Use float64 for time to avoid precision loss with timestamps
                            encoding = {
                                export_ds.name: {"dtype": "float32"},
                                "time": {"dtype": "float64"}
                            }
                            export_ds.to_netcdf(fn, engine="h5netcdf", encoding=encoding)
                        del hst_ba

                    sim_ba = qdm_arrays[-1]
                    with ProgressBar():
                        sim_ba = sim_ba.compute()

                    for yr in range(sim_year[0], sim_year[1] + 1):

                        fn = out_dir.joinpath("%s_%s_%d.nc" % (var, gcm, yr))
                        yr_slice = slice(str(yr), str(yr))
                        export_ds = sim_ba.sel(time=yr_slice)
                        export_ds = export_ds.astype("float32")
                        # Use float64 for time to avoid precision loss with timestamps
                        encoding = {
                            export_ds.name: {"dtype": "float32"},
                            "time": {"dtype": "float64"}
                        }
                        export_ds.to_netcdf(fn, engine="h5netcdf", encoding=encoding)

                    pbar.update()

    # Calculate and display elapsed time
    end_time = time.time()
    end_datetime = datetime.now()
    elapsed_seconds = end_time - start_time
    
    hours = int(elapsed_seconds // 3600)
    minutes = int((elapsed_seconds % 3600) // 60)
    seconds = int(elapsed_seconds % 60)
    
    print(f"\n{'='*60}")
    print(f"Bias Correction Completed!")
    print(f"End time: {end_datetime.strftime('%Y-%m-%d %H:%M:%S')}")
    if hours > 0:
        print(f"Total elapsed time: {hours}h {minutes}m {seconds}s")
    elif minutes > 0:
        print(f"Total elapsed time: {minutes}m {seconds}s")
    else:
        print(f"Total elapsed time: {seconds}s")
    print(f"{'='*60}\n")