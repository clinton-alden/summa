# %% [markdown]
# ## Create forcing up to current date from SNOTEL obs

# %%
# %%
import sys
sys.path.append('/home/cdalden/summa_setup/model/')
from utils import lw_clr
from utils import forcing_filler as ff
from utils import summa_check as sc


snotel = input('Enter the desired SNOTEL site code (ie. 1107:WA): ') + ':SNTL'
water_year = int(input())
out_name = input('Enter the output file name (ie. buck_WY16): ')
# out_path = input('Enter the output path (ie. ../model/forcings/): ')
out_path = '/home/cdalden/summa_setup/model/forcings/'

print('********** generating the forcing file, please be patient **********')
print('********** should take ~3 minutes to run **********')
# %% [markdown]
# ## Use metloom API to pull snotel data

# %%
import warnings
# pysumma has many depreciated packages, this ignores their warnings
warnings.filterwarnings("ignore", category=UserWarning, module='scipy')

from datetime import datetime, timedelta
from metloom.pointdata import SnotelPointData
import pandas as pd
import geopandas as gpd
import xarray as xr
from metsim import MetSim
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from metpy.units import units
import metpy.calc as mpcalc
import math
import scipy
import os
import shutil
from pytz import UTC

#  Create needed directories to store metsim run and snotel CSVs
if not os.path.exists('/home/cdalden/summa_setup/twitter_api/input/'):
    os.makedirs('/home/cdalden/summa_setup/twitter_api/input/')

if not os.path.exists('/home/cdalden/summa_setup/twitter_api/output/'):
    os.makedirs('/home/cdalden/summa_setup/twitter_api/output/')

# metsim and metloom require different formats of time and ranges, reformat here

current_day = datetime.now().day
current_day_str = str(current_day).zfill(2)
current_month = datetime.now().month
current_month_str = str(current_month).zfill(2)
current_year = datetime.now().year

water_year_str = str(current_year)
start_year = water_year - 1
start_year_str = str(start_year)

start_date = datetime(start_year, 7, 3)
end_date = datetime(current_year, current_month, current_day)

spinstart = pd.to_datetime(start_year_str + '-07-03').tz_localize('UTC')
spinend = pd.to_datetime(start_year_str + '-09-30').tz_localize('UTC')

start_loc = datetime(start_year, 10, 1).replace(tzinfo=UTC)
mask_date = datetime(start_year, 10, 2).replace(tzinfo=UTC)

dates = pd.date_range('10/01/' + start_year_str, current_month_str+'/'+current_day_str+'/' + water_year_str)

spin_range = pd.date_range('07/03/' + start_year_str, '09/30/' + start_year_str)

# %%
# Pull desired variables from snotel to dataframe
snotel_point = SnotelPointData(snotel, "MyStation")
df = snotel_point.get_hourly_data(
    start_date, end_date,
    [snotel_point.ALLOWED_VARIABLES.PRECIPITATIONACCUM, snotel_point.ALLOWED_VARIABLES.TEMP, 
     snotel_point.ALLOWED_VARIABLES.SWE, snotel_point.ALLOWED_VARIABLES.SNOWDEPTH]
)

# Specify latitude, longitude, and elevation from station metadata
lat = snotel_point.metadata.y
lon = snotel_point.metadata.x
elev = snotel_point.metadata.z

# Clean up the dataframe
df.reset_index(inplace=True)

# Rename columns
replace = {'ACCUMULATED PRECIPITATION':'accppt','AIR TEMP':'airtemp', 'datetime':'time'}
df.rename(columns=replace, inplace=True)
df.set_index('time', inplace=True)

# df.to_csv('./snotel_csvs/'+out_name+'.csv')

# add 'SNOWDEPTH' and 'SNOWDEPTH_units' to the droplist if it decides to work again
try:
    df.drop(columns=['site', 'ACCUMULATED PRECIPITATION_units', 'geometry', 'AIR TEMP_units', 'datasource', 
                 'SWE', 'SWE_units', 'SNOWDEPTH', 'SNOWDEPTH_units'], inplace=True)
except:
    df.drop(columns=['site', 'ACCUMULATED PRECIPITATION_units', 'geometry', 'AIR TEMP_units', 'datasource', 
                 'SWE', 'SWE_units'], inplace=True)
    print('SNOTEL csv has no snowdepth for this run')

# %% [markdown]
# ## Fill missing timesteps from snotel data

# %%
# Convert the index of the dataframe to a DatetimeIndex
df.index = pd.to_datetime(df.index)

# Create a date range from the first to the last timestep
date_range = pd.date_range(start=df.index.min(), end=df.index.max(), freq='h')

# Find the missing dates
missing_dates = date_range[~date_range.isin(df.index)]

# Print the missing dates
# print(missing_dates)

# Reindex the data DataFrame with the missing dates
# Concatenate the original DataFrame with a DataFrame containing the missing dates
new_df = pd.concat([df, pd.DataFrame(index=missing_dates)], axis=0)

# Sort the new DataFrame by the index
new_df = new_df.sort_index()
df = new_df

# Fill NaNs for every other column
df = df.fillna(np.nan)

# Rename index
df.index.name = 'time'

# %% [markdown]
# ## Unit Conversions

# %%
# Covert air temperature to celsius
df['airtemp'] = (df['airtemp'] - 32) * 5.0/9.0

# Convert precipitation to mm
df['accppt'] = df['accppt'] * 25.4

# Convert from geodataframe to dataframe
df = pd.DataFrame(df)


# %% [markdown]
# ## Split up data into spinup state and desired date range for MetSim

# %%
# Interpolate the missing values
df.interpolate(inplace=True)

# Seperate the data into two dataframes, before and after October 1
# spinstart = pd.to_datetime('2014-07-03').tz_localize('UTC')
# spinend = pd.to_datetime('2014-09-30').tz_localize('UTC')
spinup = df.loc[spinstart:spinend].copy()
data = df.loc[start_loc:]

# Copy the dataframe a2 to a2_copy
data_copy = data.copy()

# Create a mask to identify rows where the index is less than or equal to October 2, 2023
mask = data_copy.index <= mask_date

# Set the 'precip_accum' column to 0 for rows that satisfy the mask condition
data_copy.loc[mask, 'accppt'] = 0

# Update the value of a2 to the modified copy
data = data_copy

# Calculate the difference between the maximum value of 'precip_accum' and the previous value
spinup['pptrate'] = spinup['accppt'].cummax().diff()
data['pptrate'] = data['accppt'].cummax().diff()

# Drop accppt column
spinup.drop(columns=['accppt'], inplace=True)
data.drop(columns=['accppt'], inplace=True)

# %% [markdown]
# ## Generate SW from MetSim

# %%
# Create empty dataset
shape = (len(dates), 1, 1, )
dims = ('time', 'lat', 'lon', )

# We are running only one site, at these coordinates
lats = [lat]
lons = [lon]
elev = elev # meters
coords = {'time': dates, 'lat': lats, 'lon': lons}

# Create the initial met data input data structure
met_data = xr.Dataset(coords=coords)

# %%
for varname in ['prec', 't_min', 't_max']:
    met_data[varname] = xr.DataArray(data=np.full(shape, np.nan),
                                     coords=coords, dims=dims,
                                     name=varname)

# %%
# Resample the data to daily frequency and calculate the maximum and minimum temperatures
tmax_vals = data['airtemp'].resample('D').max()
tmin_vals = data['airtemp'].resample('D').min()

# Calculate the daily precipitation values
prec_vals = data['pptrate'].resample('D').sum()

# Interpolate the temperature values to fill in any missing days
# tmax_vals = tmax_vals.interpolate(method='linear')
# tmin_vals = tmin_vals.interpolate(method='linear')

met_data['prec'].values[:, 0, 0] = prec_vals

# Assign the daily maximum and minimum temperatures to the met_data xarray, converting to Celsius
met_data['t_min'].values[:, 0, 0] = tmin_vals
met_data['t_max'].values[:, 0, 0] = tmax_vals

met_data.to_netcdf('/home/cdalden/summa_setup/twitter_api/input/rc_forcing.nc')

# %%
# We form the domain in a similar fashion
# First, by creating the data structure
coords = {'lat': lats, 'lon': lons}
domain = xr.Dataset(coords=coords)
domain['elev'] = xr.DataArray(data=np.full((1,1,), np.nan),
                          coords=coords,
                          dims=('lat', 'lon', ))
domain['mask'] = xr.DataArray(data=np.full((1,1,), np.nan),
                          coords=coords,
                          dims=('lat', 'lon', ))

# Add the data
domain['elev'][0, 0] = elev
domain['mask'][0, 0] = 1
domain.to_netcdf('/home/cdalden/summa_setup/twitter_api/input/rc_domain.nc')

# %%
# Finally, we create the state file - the dates are 90 days prior to 
# the MetSim run dates - as usual, create an empty data structure to
# read the data into
shape = (len(spin_range), 1, 1, )
dims = ('time', 'lat', 'lon', )
coords = {'time': spin_range, 'lat': lats, 'lon': lons}
state = xr.Dataset(coords=coords)
for varname in ['prec', 't_min', 't_max']:
    state[varname] = xr.DataArray(data=np.full(shape, np.nan),
                               coords=coords, dims=dims,
                               name=varname)
    
# Resample precip to daily
prec_vals = spinup['pptrate'].resample('D').sum()

# Resample the data to daily frequency and calculate the maximum and minimum temperatures
tmax_vals = spinup['airtemp'].resample('D').max()
tmin_vals = spinup['airtemp'].resample('D').min()

# Do precip data
state['prec'].values[:, 0, 0] = prec_vals

# And now temp data and convert to C
state['t_min'].values[:, 0, 0] = tmin_vals
state['t_max'].values[:, 0, 0] = tmax_vals
state.to_netcdf('/home/cdalden/summa_setup/twitter_api/input/rc_state.nc')

# %%
# dates = pd.date_range('10/01/2014', '09/30/2015')
params = {
    'time_step'    : "60",       
    'start'        : dates[0],
    'stop'         : dates[-1],
    'forcing'      : '/home/cdalden/summa_setup/twitter_api/input/rc_forcing.nc',     
    'domain'       : '/home/cdalden/summa_setup/twitter_api/input/rc_domain.nc',
    'state'        : '/home/cdalden/summa_setup/twitter_api/input/rc_state.nc',
    'forcing_fmt'  : 'netcdf',
    'out_dir'      : '/home/cdalden/summa_setup/twitter_api/output/',
    'out_prefix': out_name,
    'scheduler'    : 'threading',
    'chunks'       : 
        {'lat': 1, 'lon': 1},
    'forcing_vars' : 
        {'prec' : 'prec', 't_max': 't_max', 't_min': 't_min'},
    'state_vars'   : 
        {'prec' : 'prec', 't_max': 't_max', 't_min': 't_min'},
    'domain_vars'  : 
        {'elev': 'elev', 'lat': 'lat', 'lon': 'lon', 'mask': 'mask'}
    }               

# Run MetSim
ms = MetSim(params)
ms.run()
output = ms.open_output().load()

# Delete MetSim input and output directories to declutter, they are unnecessary
if os.path.exists('/home/cdalden/summa_setup/twitter_api/input/'):
    shutil.rmtree('/home/cdalden/summa_setup/twitter_api/input/')

if os.path.exists('/home/cdalden/summa_setup/twitter_api/output/'):
    shutil.rmtree('/home/cdalden/summa_setup/twitter_api/output/')

# %% [markdown]
# ## Create SUMMA forcing netCDF

# %%
out_df = output.to_dataframe()
out_df.reset_index(inplace=True)
out_df.set_index('time', inplace=True)

# %%
# Remove timezone from index
data.index = data.index.tz_convert(None)

# Convert precipitation rate from m hr^-1 to kg m^-2 s^-1
data['pptrate'] = data['pptrate']/3600

# Generate relative humidity assuming T_d is overnight low temperature
# Used to calculate specific humidity and longwave radiation
ff.fill_rel_hum(data)

# Convert airtemp to Kelvin
data['airtemp'] = (1.03*(data['airtemp']-0.9)) + 273.15 # Currier snotel temp correction

# Generate pressure from hypsometric equation and site elevation (1981m)
ff.fill_pressure(data, elev)
data["airpres"] = 80000 # Set to constant value, Pa

# Generate specific humidity
ff.fill_spec_hum(data)
data['spechum'] = data['spechum'].clip(lower=0.001)


# Set shortwave radiation to MetSim output
data['SWRadAtm'] = out_df['shortwave']

# Generate longwave radiation
data['LWRadAtm'] = lw_clr.dilleyobrien1998(data['airtemp'], data['rh'])

# Can alternatively use the MetSim LW radiation
# data['LWRadAtm'] = out_df['longwave']

# Set wind to 2 m/s
data['windspd'] = 2

# Fill in missing values
data['pptrate'] = data['pptrate'].fillna(0)

# Drop unnecessary columns
data = data.drop(columns=['rh'])

# Interpolate the missing values
data.interpolate(inplace=True)

# %% [markdown]
# ## HRRR/GFS Forecast

# %%
import openmeteo_requests

import requests_cache
import pandas as pd
from retry_requests import retry

# Setup the Open-Meteo API client with cache and retry on error
cache_session = requests_cache.CachedSession('.cache', expire_after = 3600)
retry_session = retry(cache_session, retries = 5, backoff_factor = 0.2)
openmeteo = openmeteo_requests.Client(session = retry_session)

# Make sure all required weather variables are listed here
# The order of variables in hourly or daily is important to assign them correctly below
url = "https://api.open-meteo.com/v1/forecast"
params = {
	"latitude": lat,
	"longitude": lon,
	"hourly": ["temperature_2m", "relative_humidity_2m", "precipitation", "surface_pressure", "wind_speed_10m", "shortwave_radiation"],
    "wind_speed_unit": "ms",
	"timezone": "America/Los_Angeles",
	"forecast_days": 2,
	"models": "gfs_hrrr"
}
responses = openmeteo.weather_api(url, params=params)

# Process first location. Add a for-loop for multiple locations or weather models
response = responses[0]
print(f"Coordinates {response.Latitude()}°N {response.Longitude()}°E")
print(f"Elevation {response.Elevation()} m asl")
print(f"Timezone {response.Timezone()} {response.TimezoneAbbreviation()}")
print(f"Timezone difference to GMT+0 {response.UtcOffsetSeconds()} s")

# %%
# Process hourly data. The order of variables needs to be the same as requested.
hourly = response.Hourly()
hourly_temperature_2m = hourly.Variables(0).ValuesAsNumpy()
hourly_relative_humidity_2m = hourly.Variables(1).ValuesAsNumpy()
hourly_precipitation = hourly.Variables(2).ValuesAsNumpy()
hourly_surface_pressure = hourly.Variables(3).ValuesAsNumpy()
hourly_wind_speed_10m = hourly.Variables(4).ValuesAsNumpy()
hourly_shortwave_radiation = hourly.Variables(5).ValuesAsNumpy()

hourly_data = {"time": pd.date_range(
	start = pd.to_datetime(hourly.Time(), unit = "s", utc = True),
	end = pd.to_datetime(hourly.TimeEnd(), unit = "s", utc = True),
	freq = pd.Timedelta(seconds = hourly.Interval()),
	inclusive = "left"
)}
hourly_data["airtemp"] = hourly_temperature_2m
hourly_data["rh"] = hourly_relative_humidity_2m
hourly_data["precip"] = hourly_precipitation
# hourly_data["airpres"] = hourly_surface_pressure*100 # hPa to Pa
hourly_data["airpres"] = 80000 #constant air pressure, Pa
hourly_data["windspd"] = hourly_wind_speed_10m
hourly_data["SWRadAtm"] = hourly_shortwave_radiation

hourly_dataframe = pd.DataFrame(data = hourly_data)

ff.CtoK(hourly_dataframe)
ff.fill_spec_hum(hourly_dataframe)
hourly_dataframe['pptrate'] = hourly_dataframe['precip']/3600

# Generate longwave radiation
hourly_dataframe['LWRadAtm'] = lw_clr.dilleyobrien1998(hourly_dataframe['airtemp'], hourly_dataframe['rh'])

# Drop unnecessary columns
hourly_dataframe = hourly_dataframe.drop(columns=['rh', 'precip'])
# Interpolate the missing values
hourly_dataframe.interpolate(inplace=True)

# Convert the timezone-aware datetime index to naive datetime index
hourly_dataframe.set_index('time', inplace=True)
hourly_dataframe.index = hourly_dataframe.index.tz_localize(None)

# print(hourly_dataframe)

# %%
# Assuming data and hourly_dataframe are your DataFrames
# Combine them, using values from 'data' where the time index overlaps
combined_dataframe = hourly_dataframe.combine_first(data)

# If you want to ensure the index is sorted after combining
combined_dataframe = combined_dataframe.sort_index()

# %%
# Load template forcing file to preserve attributes
template = xr.open_dataset('/home/cdalden/summa_setup/model/summa_forcing_template.nc')

# Convert dataframe to xarray
dsx = combined_dataframe.to_xarray()

# Loop through variables and add attributes from template forcing file
for data_var in dsx:
    dsx[data_var].attrs = template[data_var].attrs
    
# Add hru dimension
dsx = dsx.expand_dims(dim={'hru':1})

# Add gap-filled and datastep variables
dsx['gap_filled'] = xr.DataArray(np.ones((1,dsx.time.shape[0])),dims = ['hru','time'])
dsx['data_step'] = 3600 # 3600 seconds for 1hr timesteps

# Convert all to float64
for var in dsx.data_vars:
    dsx[var] = dsx[var].astype(np.float64)

# Set hruID based on template
dsx['hruId'] = (xr.DataArray(np.ones((1))*template['hruId'].values,dims = ['hru'])).astype(np.int32)

# Transpose all variables to match SUMMA dimensions
count = 0
for var in dsx.data_vars:
    # print(var,count)
    count += 1
    if count <= 7:
        attribs = dsx[var].attrs
        arr_t = dsx[var].values.T
        dsx[var] = xr.DataArray(dims = ['time','hru'],data = arr_t)
        dsx[var].attrs = attribs

# Set hruID based on template
dsx['hruId'] = (xr.DataArray(np.ones((1))*template['hruId'].values,dims = ['hru'])).astype(np.float64).fillna(0).astype(np.int32)

# Set time to match SUMMA format and save
dsx.to_netcdf(out_path+out_name+'.nc',
                        encoding = {"time":
                                        {'dtype' : 'float64',
                                         'units' : 'hours since 1990-01-01 00:00:00',
                                         'calendar' : 'standard'}})