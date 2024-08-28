#!/bin/bash

# Loop through all model output years, 
# comparing control to warmed runs and saving as pngs

# Declare an empty array
declare -a years

# Loop from 0 to 23 (water years 2000 to 2023)
for i in {14..23}; do
    # Add leading zero if necessary and append to the array
    years+=($(printf "%02d" $i))
done

# Loop over each year in the array
for years in "${years[@]}"; do
    # Pass the hour as input to your Python script
    echo $years | python density_plots.py
done