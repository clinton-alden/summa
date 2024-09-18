import numpy as np
import pandas as pd
import xarray as xr
import datetime

run_name ='harts_+4K_WY23'
harts = xr.open_dataset('../model/output/harts_pass/template_output_'+run_name+'_timestep.nc')

depth = harts.isel(hru=0)['iLayerHeight']
var = harts.isel(hru=0)['mLayerVolFracWat']
temp = harts.isel(hru=0)['mLayerTemp']
vmask = var != -9999
dmask = depth != -9999
tmask = temp != -9999
depth.values = justify(depth.where(dmask).values)
var.values = justify(var.where(vmask).values)
temp.values = justify(temp.where(tmask).values)

# Calculate the average at all layers
average = temp.mean(dim='midToto')

# Filter var where the average is less than 273.15
filtered_var = var.where(average < 273)

# Calculate the vertical derivative
derivative = filtered_var.diff(dim='midToto')

# Initialize an empty list to store the counts
counts = []

# Loop over the 'time' dimension
for t in var.time.values:
    # Select the derivative for the current timestep
    derivative_t = derivative.sel(time=t)

    # Filter values that are greater than or equal to 0.2 or less than or equal to -0.2
    # filtered = derivative_t.where((derivative_t >= 0.1) | (derivative_t <= -0.1))
    filtered = derivative_t.where(derivative_t >= 0.1)

    # Count the number of layers with at least one such value
    count = np.isfinite(filtered).sum().values

    # Append the count to the list
    counts.append(count)

# Convert the list to a numpy array
counts = np.array(counts)

crust_days = counts.sum()/24
mean_crusts = counts.mean()

# Append netcdf
ds = xr.open_dataset('/Users/clintonalden/Documents/Research/summa_work/analysis/crust_stats.nc')

# Split the string at the underscores
parts = run_name.split("_")

# Extract the parts
site = parts[0]
model_run = parts[1]

# Extract the year and convert it to a datetime
year_str = parts[2][2:]  # Remove the 'WY' prefix
year = int(year_str) + 2000  # Convert to an integer and add 2000
date = datetime.datetime(year, 1, 1)  # Create a datetime object for the first day of the year

# Assign a value to the 'crust_days' variable at the specified coordinates
ds['crust_days'].loc[dict(time=date, model_run=model_run, site=site)] = crust_days

# Assign a value to the 'mean_crusts' variable at the specified coordinates
ds['mean_crusts'].loc[dict(time=date, model_run=model_run, site=site)] = mean_crusts

ds.to_netcdf('/Users/clintonalden/Documents/Research/summa_work/analysis/crust_stats.nc')