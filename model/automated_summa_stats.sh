#!/bin/bash

# Read the CSV file line by line, skipping the first line
tail -n +2 /home/cdalden/summa_setup/analysis/sntl_list.csv | while IFS=, read -r _ _ _ site_name _ _ _ _ _ _ _ site_code
do
    # Define your sets of inputs
    site_codes=($(for i in $(seq 2000 2024); do echo "$site_code:WA"; done))
    echo $site_codes
    water_years=($(seq 2000 2024))
    site="$site_name"
    echo $site

# Define your sets of inputs
# site_codes=($(for i in $(seq 2000 2024); do echo "791:WA"; done))
# water_years=($(seq 2000 2024))
# site="stevenspass"

    # loop through warming scenarios
    declare -a temp_runs=("current" "+1K" "+2K" "+3K" "+4K")

    for run in "${temp_runs[@]}"; do
        output_names=($(for i in $(seq 2000 2024); do echo "${site}_${run}_WY${i:2:2}"; done))

        # Get the length of the arrays
        length=${#site_codes[@]}

        # Loop over the arrays
        for ((i=0; i<$length; i++)); do
        # Get the inputs for this iteration
        input1=${site_codes[$i]}
        input2=${water_years[$i]}
        input3=${output_names[$i]}

        # Provide input to the first Python script
        echo -e "$input1\n$input2\n$input3" | python3 snotel_to_pysumma_${run}.py

        echo "forcing file created, now running summa"

        # Replace the content of the text file with input3
        echo "'$input3.nc'" > forcings/forcing_file_list.txt

        # Run the second Python script
        python3 run_summa.py

        echo "summa run for $input3 complete, check output folder for output file"

        done
    done

    # Directory containing the output files
    dir="../model/output/"

    # Loop through all .nc files in the directory
    for file in "$dir"/*.nc; do
        # Extract the filename without the path and extension
        filename=$(basename "$file" .nc)

        # Split the filename at the underscores
        IFS='_' read -ra parts <<< "$filename"

        # Extract the required parts and format the input string
        input="${parts[1]}_${parts[2]}_${parts[3]}"
        echo $input

        # Execute the Python script with the input string
        echo "$input" | python /home/cdalden/summa_setup/analysis/crust_stats.py 
    done

    # remove all output files for this site after stats are calculated
    rm ../model/output/output_${site}*
    
done