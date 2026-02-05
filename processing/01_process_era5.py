import tqdm
import dask
import xarray as xr
import xclim as xc
import numpy as np
from pathlib import Path
from config import DATA_DIR, OUT_DIR, era5_years, geog_bbox
from utils import *

# Suppress dask warnings on chunk size
dask.config.set({"array.slicing.split_large_chunks": False})


# Function to calculate daily maximum temperature from hourly data
def _calc_tasmax(da: xr.DataArray, dask_compute=False) -> xr.DataArray:
    da = xc.units.convert_units_to(da, "degC")  # [degC]
    tasmax = da.resample(time="D").max()
    tasmax = tasmax.rename("tasmax")
    tasmax.attrs["long_name"] = "Daily maximum 2 metre temperature"
    if dask.is_dask_collection(tasmax) & dask_compute:
        tasmax = tasmax.compute()
    return tasmax


# Function to calculate daily total precip from hourly data
def _calc_pr(da: xr.DataArray, dask_compute=False) -> xr.DataArray:
    da = xc.units.convert_units_to(da, "mm")  # [mm]
    pr = da.resample(time="D").sum()
    pr = pr.rename("pr")
    pr = pr.assign_attrs({"units": "mm/day"})
    if (dask.is_dask_collection(pr)) & dask_compute:
        pr = pr.compute()
    return pr


# Function to calculate daily mean horizontal near-surface wind speed from
# hourly data. Thhis uses u10 and v10 vectors to calculate mean wind speed
def _calc_sfcWind(ds: xr.Dataset, dask_compute=False) -> xr.DataArray:
    u10 = ds["u10"]
    v10 = ds["v10"]

    u10_daily = u10.resample(time="D").mean()
    v10_daily = v10.resample(time="D").mean()

    sfcWind = np.sqrt(u10_daily**2 + v10_daily**2)  # [m s**-1]
    sfcWind.name = "sfcWind"
    sfcWind.attrs = {
        "units": u10.units,
        "long_name": "Average 10 metre horizontal wind speed",
    }

    sfcWind = xc.units.convert_units_to(sfcWind, "km/hour")  # [km hr**-1]

    if (dask.is_dask_collection(sfcWind)) & dask_compute:
        sfcWind = sfcWind.compute()
    return sfcWind


# Function to calculate vapor pressure from temperature. Used to calculate relative humidity.
def _vapor_pressure(T: xr.DataArray) -> xr.DataArray:
    if xc.units.str2pint(T.units).units == "kelvin":
        T = xc.units.convert_units_to(T, "degC")

    vp = 611.2 * np.exp((17.62 * T) / (243.12 + T))  # [Pa]
    vp = vp * 0.001  # [kPa]

    vp.name = "esat"
    vp.attrs = {"units": "kPa", "long_name": "Vapor pressure"}
    return vp


# Calculate near-surface daily minimum relative humidity from hourly air and
# dewpoint temperature. Clip values to 0-100%.
def _calc_hursmin(ds: xr.Dataset, dask_compute=False) -> xr.DataArray:
    tas = ds["t2m"]
    tdps = ds["d2m"]

    hurs = 100 * _vapor_pressure(tdps) / _vapor_pressure(tas)  # [%]
    hurs.attrs["long_name"] = "2 metre relative humidity"
    hurs.name = "hurs"
    hurs.attrs = {"units": "%", "long_name": "2 metre minimum relative humidity"}

    hursmin = hurs.resample(time="D").min()
    hursmin = hursmin.clip(min=0.0, max=100.0)
    hursmin.name = "hursmin"

    if (dask.is_dask_collection(hursmin)) & dask_compute:
        hursmin = hursmin.compute()
    return hursmin


# Organize ERA5 dataset
def _organize_era5(ds: xr.Dataset) -> xr.Dataset:
    global_attrs = {
        "source_id": "ERA5",
        "url": "https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5",
        "history": "Original data downloaded Jan 7-13,2022 in netcdf format. \
        Converted to netcdf4 using netCDF library version 4.8.1.",
    }

    ds = ds.rename(
        {
            "longitude": "lon",
            "latitude": "lat",
        }
    )

    # Clip geographic limits
    ds = trim_geolims(ds, geog_bbox)

    # Convert calendar to 365 day years or 'noleap'
    ds = ds.convert_calendar("noleap")

    # Assign axes labels
    ds["lon"].attrs.update({"standard_name": "longitude", "axis": "X"})
    ds["lat"].attrs.update({"standard_name": "latitude", "axis": "Y"})
    ds["time"].attrs.update({"standard_name": "time", "axis": "T"})

    ds.attrs = global_attrs

    return ds


def process_era5(src: str, dask_load=False) -> xr.Dataset:

    ds = xr.open_mfdataset(
        src,
        parallel=dask_load,
        engine="h5netcdf",
        mask_and_scale=True,
    )

    # Get variable names for given dataset
    varnames = get_var_names(ds)

    # Identify which variable is in current dataset
    tasmax_bool = False
    pr_bool = False
    sfcWind_bool = False
    hursmin_bool = False

    if any([x == "t2m" for x in varnames]):
        tasmax_bool = True

    if any([x == "tp" for x in varnames]):
        pr_bool = True

    if any([x == "u10" for x in varnames]) and any([x == "v10" for x in varnames]):
        sfcWind_bool = True

    if any([x == "t2m" for x in varnames]) and any([x == "d2m" for x in varnames]):
        hursmin_bool = True

    export_dict = {}

    # Create xarray dataset object, easier to do these daily summaries using groupby.
    if tasmax_bool:
        tasmax = _calc_tasmax(ds["t2m"], dask_compute=dask_load)
        export_dict.update({"tasmax": tasmax})

    if pr_bool:
        pr = _calc_pr(ds["tp"], dask_compute=dask_load)
        export_dict.update({"pr": pr})

    if sfcWind_bool:
        sfcWind = _calc_sfcWind(ds[["u10", "v10"]], dask_compute=dask_load)
        export_dict.update({"sfcWind": sfcWind})

    if hursmin_bool:
        hursmin = _calc_hursmin(ds[["t2m", "d2m"]], dask_compute=dask_load)
        export_dict.update({"hursmin": hursmin})

    daily_ds = xr.merge([export_dict])  # Put dict in brackets to make iterable
    daily_ds = _organize_era5(daily_ds)  # Final re-arrangement of coords

    return daily_ds


if __name__ == "__main__":

    era5_dir = Path(DATA_DIR) / "era5"
    out_dir = Path(OUT_DIR) / "processed" / "era5"

    start_year, end_year = era5_years
    year_range = range(start_year, end_year + 1)

    with tqdm(total=len()) as pbar:
        for year in year_range:
            src = list(era5_dir.glob("%d*.nc" % year))
            ds = process_era5(src, dask_load=True)
            for var in list(ds.keys()):
                fn = out_dir / ("%s_era5_%d.nc" % (var, year))
                ds[var].to_netcdf(fn, engine="h5netcdf")
            ds.close()
            pbar.update()
