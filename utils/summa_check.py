''' 
Script for checking SUMMA forcing files before model initialization 
to ensure no errors in the meteorological data.

To use - run summa_check.final('name_of_site', 'state_abbreviation') to run all checks
The name_of_site variable should be the name of the directory you are running from
ie: Running from `/home/buckinghorse/ ` in Washington state
--> `summa_check.final('buckinghorse', 'WA')


Includes:
- Check for NaN values
- Check to ensure meteorological variables are within realistic ranges
- Check local attibutes file to ensure latitude, longitude, and elevation are 
within acceptable ranges for the state the site is in


Clinton Alden
University of Washington
Last edit - January 2024
'''


import numpy as np
import pandas as pd
import xarray as xr
import io
import sys


def find_forcing(site):
    # Input site name as string to return the forcing file xarray
    
    directory = '/home/jovyan/data/' + site
    forcing = open(directory + '/forcings/forcing_file_list.txt').readline().strip("'")
    data = xr.open_dataset(directory + '/forcings/' + forcing)
    
    return data


def nan_check(site):
    # Check forcing file for NaN values and return number of NaNs per variable

    data = find_forcing(site)
    
    nan_variables = {}
    
    for variable in data.variables:
        nan_count = data[variable].isnull().sum().item()
        if nan_count > 0:
            nan_variables[variable] = nan_count
            # flag = 1
    
    if nan_variables:
        print("***Variables with NaN values:")
        for variable, count in nan_variables.items():
            print(f"{variable}: {count} NaN values")
        print('To fill NaN values with mean values, run `summa_check.fill_nans(site_name)`')
        print(' ')
    else:
        print("That's some lonely tikka masala.")
        print("No Na(a)N values found in any variable.")
        print(" ")
        

        
        
def fill_nans(site):
    # Fill NaN values with the mean value for each variable
    # *** Want to update later to fill with realistic values for time of year***
    
    data = find_forcing(site)
    
    for variable in data.variables:
        mean_value = data[variable].mean().item()
        data[variable] = data[variable].fillna(mean_value)
        
    return data
        

def check_range(site, variable_name, min_value, max_value):
    
    data = find_forcing(site)
    
    variable = data[variable_name].values

    # Check if values are within the desired range
    within_range = (variable >= min_value) & (variable <= max_value)
    outside_range_count = len(within_range) - within_range.sum()

    if within_range.all():
        print(f"All values of '{variable_name}' are physically realistic.")
    elif within_range.any():
        print(f"***Some values of '{variable_name}' are NOT physically realistic.")
        print(f"Number of values outside the range: {outside_range_count}")
        flag = 1
    else:
        print(f"***No values of '{variable_name}' are within the desired range.")
        flag = 1
        
        
def met_range_check(site):
    
    
    data = find_forcing(site)
    
    # Check airtemp, bounded between 227K (-50F) and 322K (120F)
    # Coldest/warmest temp in WA on record = -48F/120F
    # May need to be updated for other climates (ie. Alaska)
    check_range(site, 'airtemp', 227, 322)
    
    # Check specific humidity, bounded where 0.1g/kg (~-40F dewpoint) and 20g/kg (~75F dewpoint)
    check_range(site, 'spechum', 0.0001, 0.020)
    
    # Check air pressure, bounded between 1080hPa (~max recorded on earth) and 500hPa (~18,000' ASL)
    check_range(site, 'airpres', 50000, 108000)
    
    # Check wind speed, bounded between 0 m/s and 50 m/s
    check_range(site, 'windspd', 0, 50)
    
    # Check precipitation rate, bounded between 0 and 0.0212 kg m^-2 s^-1 (3 in/hr)
    # 3 in/hr chosen as roughly the 1/1000 year 1 hr event for comprable maritime mountain climates in W. US
    # from NOAA Precipitation Frequency Data Server: https://hdsc.nws.noaa.gov/pfds/pfds_map_cont.html?bkmrk=ca
    check_range(site, 'pptrate', 0, 0.0212)
    
    # Check incoming shortwave radiation, bounded between 0 and 1361 W m^2 (solar constant, probabaly too high)
    check_range(site, 'SWRadAtm', 0, 1361)
    
    # Check incoming longwave radiation, bounded between 0 and 1000 W m^2 (not totally sure what this should be)
    check_range(site, 'LWRadAtm', 0, 1000)
    print(' ')
    
    
def attrs_check(site, state):
    # Check local attributes file to ensure latitude, longitude, and elevation are set correctly
    # Elevation check is not great as it is just a range for each state from min to max but will check for NaN, -, etc.
    
    directory = '/home/jovyan/data/' + site
    attrs = xr.open_dataset(directory + '/params/local_attributes.nc')
    
    lat = attrs['latitude'].values[0]
    lon = attrs['longitude'].values[0]
    elev = attrs['elevation'].values[0]
    
    state_info = pd.read_csv('/home/jovyan/data/lib/attr_check_states.csv')
    
    max_lon_value = state_info.loc[state_info['state'] == state, 'max_lon'].values[0]
    min_lon_value = state_info.loc[state_info['state'] == state, 'min_lon'].values[0]
    max_lat_value = state_info.loc[state_info['state'] == state, 'max_lat'].values[0]
    min_lat_value = state_info.loc[state_info['state'] == state, 'min_lat'].values[0]
    max_z_value = state_info.loc[state_info['state'] == state, 'max_z'].values[0]
    min_z_value = state_info.loc[state_info['state'] == state, 'min_z'].values[0]
    

    # Check if values are within the desired range
    within_lat = (lat >= min_lat_value) & (lat <= max_lat_value)
    within_lon = (lon >= min_lon_value) & (lon <= max_lon_value)
    
    if (within_lon.all() & within_lat.all()):
        print('Lat/Lon coordinates match specified state')
    else:
        print('***Lat/Lon coordinates do NOT match specified state')
        
    within_z = (elev >= min_z_value) & (elev <= max_z_value)
    if within_z.all():
        print('Elevation matches specified state')
    else:
        print('***Elevation does NOT match specified state')
    
    
def final(site, state):
    # Run all checks before initializing model
    # Specify site name (current directory) and state abbreviation both as strings
    
    # Redirect stdout to capture the printed output
    captured_output = io.StringIO()
    sys.stdout = captured_output
    
    nan_check(site)
    
    met_range_check(site)
    
    attrs_check(site, state)
    
    # Reset stdout to the original value
    sys.stdout = sys.__stdout__
    
    # Get the captured output as a string
    output_string = captured_output.getvalue()
    print(output_string)

    # Check if the character '*' is found in the captured output
    if '*' in output_string:
        print('\033[1m' + "Hold up, check the attributes and/or forcing files to fix aforementioned issues")
        print("If there are issues with missing data, run `summa_check.fill_nans('site')`")
    else:
        print("No issues with the forcing file found, SUMMA to your heart's content")
    

    




    
