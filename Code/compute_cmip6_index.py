# JS 2025-04-07

# load python modules
import numpy as np
import xarray as xr
import pandas as pd
import csv
import matplotlib.dates as mdates
# from datetime import datetime
import netCDF4 as nc
import scipy 
import dask
from dask.diagnostics import ProgressBar
import os
import re
import glob
import cftime
from timeit import default_timer as timer
import sys
# load required modules
from datetime import timedelta
import argparse

def parser_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--m",
        type=int,
        help="Model index from the modlist to process")
    return p.parse_args()

args = parser_args()
m = args.m

# print run time
def print_run_time(time):
	print('Elapsed (wall) time: ' + str(timedelta(seconds=time)))

def timestamp_to_cftime(ds):
    if isinstance(ds, pd.Timestamp):
        # Create a cftime.DatetimeProlepticGregorian object from a Timestamp
        return cftime.DatetimeProlepticGregorian(ds.year, ds.month, ds.day, ds.hour, ds.minute, ds.second)
    return ds  # Return as is if not a Timestamp
    
# allow large chunks and silence the warning
# https://docs.dask.org/en/latest/array-slicing.html
dask.config.set({"array.slicing.split_large_chunks": False})

# set paths and filenames
inpath = '/cmip6_rawdata/'   # location of model source list
outpath = '/cmip6_sa_summer/'  # output directory

# define which scenarios to use
scen_h = 'historical'
scen_f = 'ssp585'        # ssp126 or ssp585

# read model lists
# model list .txt files should have OBS path as first entry
modfile_h = open(inpath + 'cmip6_daily_' + scen_h + '.txt', "r")
modfile_f = open(inpath + 'cmip6_daily_' + scen_f + '.txt', "r")
modlist_h = modfile_h.read().splitlines()
modlist_f = modfile_f.read().splitlines()

# ds_land = xr.open_dataset('/g/data/ng72/js5018/sst/cmip6/aus_land.nc')
# land_mask = (ds_land.AUS == 0)

## ============== CHANGE ARGUMENT ABOVE!!
# split path names to retrieve model parameters
modcode_h = modlist_h[m].split("/")  # [0] as we just use the first one path currently
modcode_f = modlist_f[m].split("/")
mod_name = modcode_f[8]
variant = modcode_f[10]

modcode = mod_name + '.' + scen_h + '+' + scen_f + '.' + variant   ## CHANGE THIS!!! 
print('Processing... ' + modcode)

# read netcdf filenames from historical and scenario paths
filenames_h = sorted(glob.glob(modlist_h[m] + '*.nc'))
filenames_f = sorted(glob.glob(modlist_f[m] + '*.nc'))
infiles = filenames_h + filenames_f
infiles = list(dict.fromkeys(infiles))  # remove duplicates, for OBS?

# load data
print('Loading data... ')
ds = xr.open_mfdataset(infiles, combine='nested', concat_dim='time')

_, index = np.unique(ds['time'], return_index=True)   # set of unique times
ds = ds.isel(time=index)   # keep only unique times

# rename tos variable in CMIP6 models to sst
ds = ds.rename_vars({'tos':'sst'})   ### CHANGE THIS VARIABLE!!!!

idx_reg = {'W': 100, 'E': 165, 'S': -50, 'N': -20}

# rename coords if not lat/lon
if hasattr(ds,'longitude'):
    if hasattr(ds,'lon'):
        print('Warning: this model has both lon and longitude arrays..')
        print('Deleting longitude and latitude...')
        ds = ds.drop_vars({'latitude','longitude'})
    else:
        ds = ds.rename_vars({'longitude':'lon','latitude':'lat'})
elif hasattr(ds,'nav_lon'):
    ds = ds.rename_vars({'nav_lon':'lon','nav_lat':'lat'})

# store min and max lon, may be required for rotated lon coords
min_lon = np.min(ds.lon.values)
max_lon = np.max(ds.lon.values)

print(min_lon)
print(max_lon)

# ===================== trim to required region =========================
print('select the region... ')
# first case: lon is [0:360]
if (max_lon > idx_reg['E']):
    mask = ((ds.lon > idx_reg['W']) & (ds.lon < idx_reg['E'])& (ds.lat > idx_reg['S']) & (ds.lat < idx_reg['N'])).compute()
    ds_region = ds.where(mask, drop=True)

# second case: lon is [-180:180]
elif (max_lon > idx_reg['W']) & (max_lon < 200):
    mask = (((ds.lon > idx_reg['W']) & (ds.lon < idx_reg['E']-360)) & (ds.lat > idx_reg['S']) & (ds.lat < idx_reg['N'])).compute()
    ds_region = ds.where(mask, drop=True)
    ds_region = ds_region.assign_coords(lon=(ds_region.lon % 360))  # wrap to [0:360]

# third case: lon is [-300:60]
elif (max_lon < idx_reg['W']):  
    mask = ((ds.lon > idx_reg['W']-360) & (ds.lon < idx_reg['E']-360) & (ds.lat > idx_reg['S']) & (ds.lat < idx_reg['N'])).compute()
    ds_region = ds.where(mask, drop=True)
    ds_region = ds_region.assign_coords(lon=(ds_region.lon % 360))  # wrap to [0:360]

# clear all variables except 'sst'
varlist = list(ds_region.data_vars)   # first get the variable list
varlist.remove('sst')  # remove 'sst' from the list, which will be kept
ds_region = ds_region.drop_vars(varlist)   


# ==================== trim time, and store start and end of time array as labels ======================
# first case: time array is numpy.datetime64
if np.issubdtype(ds_region.time.dtype, np.datetime64):
    print('Time array: numpy.datetime64')
    ds_region_period = ds_region.where(ds_region.time >= np.datetime64('1989-12-01T00:00:00'), drop=True)
    ds_region_period = ds_region_period.where(ds_region_period.time <= np.datetime64('2100-12-31T23:59:59'), drop=True)

    tlab1 = np.datetime_as_string(ds_region_period.time.values[0], unit='Y')
    tlab2 = np.datetime_as_string(ds_region_period.time.values[-1], unit='Y')
    date_start = np.datetime_as_string(ds_region_period.time.values[0])
    date_end = np.datetime_as_string(ds_region_period.time.values[-1])

else:
    # second case: time array is 365Day calendar
    if type(ds_region.time.values[0]) is cftime.DatetimeNoLeap:
        print('Time array: 365Day calendar')
        ds_region_period = ds_region.where(ds_region.time >= cftime.DatetimeNoLeap(1989,12,1,0,0,0), drop=True)
        ds_region_period = ds_region_period.where(ds_region_period.time <= cftime.DatetimeNoLeap(2100,12,31,23,59,59), drop=True)

    # third case: time array is 360Day calendar
    elif type(ds_region.time.values[0]) is cftime.Datetime360Day:
        print('Time array: 360Day calendar')
        ds_region_period = ds_region.where(ds_region.time >= cftime.Datetime360Day(1989,12,1,0,0,0), drop=True)
        ds_region_period = ds_region_period.where(ds_region_period.time <= cftime.Datetime360Day(2100,12,30,23,59,59), drop=True)

    elif type(ds_region.time.values[0]) is pd._libs.tslibs.timestamps.Timestamp:
        print('Time array: pandas timestamps')
        # Convert the time coordinate to cftime if it contains Timestamp objects
        ds_region['time'] = xr.apply_ufunc(timestamp_to_cftime,
                                           ds_region['time'], 
                                           vectorize=True) # Apply the function element-wise
        ds_region_period = ds_region.where(ds_region.time >= cftime._cftime.DatetimeProlepticGregorian(1989,12,1,0,0,0), drop=True)
        ds_region_period = ds_region_period.where(ds_region_period.time <= cftime._cftime.DatetimeProlepticGregorian(2100,12,30,23,59,59), drop=True)

    else:
        print('Time array: Not determined (ERROR!)')

    tlab1 = ds_region_period.time.values[0].strftime('%Y')
    tlab2 = ds_region_period.time.values[-1].strftime('%Y')
    date_start = ds_region_period.time.values[0].strftime()
    date_end = ds_region_period.time.values[-1].strftime()

# year span label
tlab = tlab1 + '-' + tlab2
print(tlab)
print('Start: ' + date_start)
print('End: ' + date_end)
## ==================== mask the land ================================
print('select the period and mask the land... ')
ds_region_period['sst'] = ds_region_period.sst.where(ds_region_period.sst > 0)

# =======================Select summer months ================================
print('Computing the summer average... ')
ds_region_period_summer = ds_region_period.sel(time=ds_region_period.time.dt.month.isin([12,1,2,3,4]))

## ======================= write to file ======================================
# ds_final = xr.Dataset({"sst": ds_region_average})
print('Computing index and writing to file... ')
outfile = outpath + 'summerdailysst.' + modcode + '.nc'  # CHANGE THIS!
ds_region_period_summer.to_netcdf(outfile)


del ds, ds_region, ds_region_period, mask, ds_region_period_summer
import gc; gc.collect()