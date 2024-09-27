#!/bin/bash

# Directory containing the .nc files
dir="../model/output/harts_pass/"

# Loop through all .nc files in the directory
for file in "$dir"/*.nc; do
    # Extract the filename without the path and extension
    filename=$(basename "$file" .nc)

    # Split the filename at the underscores
    IFS='_' read -ra parts <<< "$filename"

    # Extract the required parts and format the input string
    input="${parts[2]}_${parts[3]}_${parts[4]}"
    echo $input

    # Execute the Python script with the input string
    echo "$input" | python crust_stats.py 
done