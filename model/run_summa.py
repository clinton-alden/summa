# %%
import xarray as xr
import os
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import pysumma as ps
import pysumma.plotting as psp
import warnings

# pysumma has many depreciated packages, this ignores their warnings
warnings.filterwarnings("ignore", category=UserWarning)

# %%
executable = 'summa.exe'
filemanager = '/home/cdalden/summa_setup/model/settings/file_manager_summa.txt'

# Create a pySUMMA simulation object
s = ps.Simulation(executable, filemanager)


# Set the simulation start and end times from forcing file
with open('./forcings/forcing_file_list.txt', 'r') as file:
    # Read the contents
    forcing_file = file.read().replace("'", "")
forcing = xr.open_dataset('./forcings/'+forcing_file.strip(), engine='netcdf4')
time = forcing['time']

dt64 = np.datetime64(time.isel(time=0).values)
dt = pd.to_datetime(dt64)
start = dt.strftime('%Y-%m-%d %H:%M')
dt64 = np.datetime64(time.isel(time=-1).values)
dt = pd.to_datetime(dt64)
end = dt.strftime('%Y-%m-%d %H:%M')

s.manager['simStartTime'] = start
s.manager['simEndTime'] = end

# Add in some additional variables so we can demonstrate plotting capabilities
output_settings = {'period': 1, 'instant': 1, 'sum': 0, 
              'mean': 0, 'variance': 0, 'min': 0, 'max': 0}
layer_vars = ['mLayerTemp', 'mLayerDepth', 'mLayerHeight',
              'mLayerLiqFluxSoil', 'mLayerVolFracIce', 'mLayerVolFracLiq', 
              'mLayerVolFracWat','mLayerMatricHead', 'iLayerHeight', 'scalarSnowDepth', 'nSnow']

# Create the new variables
for var in layer_vars:
    s.output_control[var] = output_settings

# Ensure all variables have the same statistics
all_vars = set(layer_vars + [o.name for o in s.output_control.options])
for var in all_vars:
    s.output_control[var] = output_settings

# Set params
s.decisions['snowLayers'] = 'jrdn1991'
s.decisions['thCondSnow'] = 'tyen1965'
s.decisions['snowDenNew'] = 'hedAndPom'
s.decisions['compaction'] = 'consettl'
s.decisions['astability'] = 'mahrtexp'

s.global_hru_params['tempCritRain'] = 273.15
s.global_hru_params['newSnowDenMin'] = 50
s.global_hru_params['densScalGrowth'] = 0.10
s.global_hru_params['densScalOvrbdn'] = 0.025
s.global_hru_params['fixedThermalCond_snow'] = 0.35
s.global_hru_params['Fcapil'] = 0.04 # intially 0.06 for salmon
# s.global_hru_params['albedoDecayRate'] = 1.0d+5

out_name = os.path.splitext(forcing_file)[0]
# Run the model, specify the output suffix
print('********** MODEL INITIALIZED **********')
print('********** estimated runtime ~30 seconds **********')
s.run('local', run_suffix=out_name)

# %%
print('Model status:', s.status)


# Create and save density and temp profile plots
out = s.output

# define variables to be plotted from output
depth = out.isel(hru=0)['iLayerHeight']
temp = out.isel(hru=0)['mLayerTemp']
frac_wat = out.isel(hru=0)['mLayerVolFracWat']

# # plot snow depth and layer temperature
# plot = psp.layers(temp, depth, colormap='viridis', plot_soil=False, plot_snow=True, add_colorbar=True)
# out['scalarSnowDepth'].plot(color='red', linewidth=2)
# plt.title(out_name+' Temperature and Snow Depth')
# plt.xlabel('Date')
# plt.ylabel('Snow Depth [m]')
# plt.grid(ls='--', alpha=0.5)

# # plt.savefig('./output/plots/'+out_name+'_layer_temp.png', dpi=200)
# plt.close()

# # plot snow depth and layer density
plot = psp.layers(frac_wat, depth, colormap='viridis', plot_soil=False, plot_snow=True, add_colorbar=True)
out['scalarSnowDepth'].plot(color='red', linewidth=2)
plt.title(out_name+' Density and Snow Depth')
plt.xlabel('Date')
plt.ylabel('Snow Depth [m]')
plt.grid(ls='--', alpha=0.5)

plt.savefig('./output/plots/'+out_name+'_layer_density.png', dpi=200)
plt.close()
