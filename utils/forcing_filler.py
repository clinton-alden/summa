''' 
Define functions to convert units and fill missing data for SUMMA pre-processing

List of functions:
    - Relative humidity calculation from temperature (Tdew from Tmin)
        rel_hum
    - Fahrenheit to Celsius conversion
        FtoC
    - Celsius to Kelvin conversion
        CtoK
    - Fahrenheit to Kelvin converstion
        FtoK
    - Specific humidity calculation from relative humidity, air temperature, and air pressure
        spec_hum
    - Air pressure calculation from air temperature, relative humidity, and site elevation
        pressure
    - Precipitation rate and unit conversion from accumulated precipitation
        precip_rate
    - Miles per hour to meters per second conversion
        MPHtoMPS
        

Clinton Alden - last updated March 2024
University of Washington
'''
import numpy as np
import pandas as pd
import math

def fill_rel_hum(df):
    '''
    Create a function to filling missing relative humidity values from estimated dewpoint using method 
    described in Lawrence (2005). Dewpoint is defined as minimum temperature in over the preceeding 12 
    hours or following 12 hours (whichever value is lower) similar to findings from Running et al. (1987).
    
    Lawrence, M. G. (2005). The relationship between relative humidity 
    and the dewpoint temperature in moist air: A simple conversion and applications. 
    Bulletin of the American Meteorological Society, 86(2), 225-234.
    
    Running, S. W., Nemani, R. R., & Hungerford, R. D. (1987). 
    Extrapolation of synoptic meteorological data in mountainous terrain and its use for simulating 
    forest evapotranspiration and photosynthesis. Canadian Journal of Forest Research, 17(6), 472-483.
    
    NOTE - This inputs air temperature in degrees Celsius

    Clinton Alden - October 2023
    '''
    
    preceding_rows = 12 # Observations are hourly, so previous 12 obs is previous 12 hours
    following_rows = 12

    # Calculate the minimum temperature in the preceding and following rows
    df['min_temp_preceding'] = df['airtemp'].shift(1).rolling(window=preceding_rows, min_periods=1).min()
    df['min_temp_following'] = df['airtemp'].shift(-1).rolling(window=following_rows, min_periods=1).min()

    # Initialize the t_d and rh columns
    df['t_d'] = 0
    df['rh'] = 0

    # Determine the dewpoint temperature from the minimum temperature between preceding and following 12 hours.
    df['t_d'] = df[['min_temp_preceding', 'min_temp_following']].min(axis=1)

    # Fill missing relative humidity obs based on air and dewpoint temperature using Lawrence (2005) method.
    df['rh'] = 100 - 5 * (df['airtemp'] - df['t_d'])

    df.drop(columns=['min_temp_preceding', 'min_temp_following', 't_d'], inplace=True)
    
    return(df)
            

def FtoC(df):
    # Create a function to convert temperature from Fahrenheit to Celsius
    df['airtemp'] = (df['airtemp']-32)*(5/9)
    return(df)


def CtoK(df):
    # Create a function to convert temperature from Celsius to Kelvin
    df['airtemp'] = df['airtemp']+273.15
    return(df)
    
    
def FtoK(df):
    # Create a function to convert temperature from Fahrenheit to Celsius
    df['airtemp'] = (df['airtemp']-32)*(5/9)+273.15
    return(df)

    
def fill_spec_hum(df):
    # Calculate specific humidity from relative humidity, air temperature in K, and air pressure.
     
    rh = df['rh']
    T = df['airtemp']
    p = df['airpres']

    T0 = 273.15

    spechum = rh * np.exp((17.67*(T-T0))/(T-29.65)) / (0.263*p)
    df['spechum'] = spechum

    return(df)


def fill_pressure(df, elevation):
    '''This function first derives a typical pressure value for the elevation of the 
    grid using the hypsometric equation and atmospheric scale height for midlatitudes.
    Then, using this typical pressure value it computes mixing ratio and virtual
    temperature. Using these calculated values, it then recomputes air pressure using
    the hypsometric equation.
    
    Clinton Alden - last updated March 2024
    '''
    
    # Define constants
    p_0 = 101325 # Standard SLP in pascals
    g = 9.81 # m s^-2
    R_d = 287.053 # J K^-1 kg^-1
    z = elevation
    T_0 = 273 # K
    R_v = 461 # J K^-1 kg^-1
    L_v = 2.5 * (10**6) # J kg^-1
    e_0 = 611
    H = 8000 # Scale height of atmosphere for mid-latitudes in meters
    
    # Define variables from dataframe
    T = df.airtemp
    rh = df.rh
    
    # Calculate typical pressure value from hypsometric equation using scale height
    p_typ = p_0 * math.exp(-z/H)
    #print(p_typ)
    
    # Initialize new columns
    df['e_s'] = np.zeros(len(df), dtype=float)
    df['e'] = np.zeros(len(df), dtype=float)
    df['w'] = np.zeros(len(df), dtype=float)
    df['T_v'] = np.zeros(len(df), dtype=float)
    df['airpres'] = np.zeros(len(df), dtype=float)
    
    for i in df.index:
        # Calculate vapor pressure
        df.at[i, 'e_s'] = e_0 * math.exp((L_v/R_v)*((1/T_0)-(1/df.at[i, 'airtemp'])))
        df.at[i, 'e'] = (df.at[i, 'rh'] * df.at[i, 'e_s'])/100
    
        # Calculate mixing ratio
        df.at[i, 'w'] = (0.622*df.at[i, 'e'])/(p_typ-df.at[i, 'e'])
    
        # Calculate virtual temperature
        df.at[i, 'T_v'] = df.at[i, 'airtemp'] * (1 + 0.608*df.at[i, 'w'])
    
        # Compute air pressure
        df.at[i, 'airpres'] = p_0 * math.exp((-z*g)/(R_v*df.at[i, 'T_v']))
        
    df.drop(columns=['e_s', 'e', 'w', 'T_v'], inplace=True)
    
    return(df)
    

def precip_rate(df):
    # Calculate precip rate from accumulated precip and convert to valid SUMMA units (kg m^-2 s^-1)
    
    # Set erroneous negative precip values to zero
    df.acc_precip[df.acc_precip < 0] = 0

    # Calculate hourly precip from accumulated
    df['pptrate'] = df.acc_precip.diff()

    # Set negative accumulations to zero
    df.pptrate[df.pptrate < 0] = 0
    
    # Convert precip from in/hr to m/hr
    df['pptrate'] = df['pptrate']/39.37

    # Convert precipitation rate from m hr^-1 to kg m^-2 s^-1
    df['pptrate'] = df['pptrate']/3.6
    
    df.drop(columns='acc_precip', inplace=True)
    
    return(df)


# def sw_rad(df):
    ''' 
    Calculate incoming shortwave radiation using the method described in Hargreaves and Samani (1982). Shortwave
    radiation is computed as a function of latitude, day of year, and daily air temperature range.
    
    Hargreaves, G. H., & Samani, Z. A. (1982). Estimating potential evapotranspiration. Journal of the irrigation and
    Drainage Division, 108(3), 225-230.
    '''
    
    
def MPHtoMPS(df):
    ''' Convert wind speed in miles per hour to meters per second '''
    
    df['windspd'] = df['windspd']/2.237
    
    return(df)
    


