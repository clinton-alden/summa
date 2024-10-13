#!/bin/bash

# Define your sets of inputs
site_codes=($(for i in $(seq 2000 2024); do echo "728:WA"; done))
water_years=($(seq 2000 2023))
output_names=($(for i in $(seq 2000 2024); do echo "salmon_+3K_WY${i:2:2}"; done))

# Get the length of the arrays
length=${#site_codes[@]}

# Loop over the arrays
for ((i=0; i<$length; i++)); do
  # Get the inputs for this iteration
  input1=${site_codes[$i]}
  input2=${water_years[$i]}
  input3=${output_names[$i]}

  # Provide input to the first Python script
  echo -e "$input1\n$input2\n$input3" | python3 snotel_to_pysumma_+3K.py

  echo "forcing file created, now running summa"

  # Replace the content of the text file with input3
  echo "'$input3.nc'" > forcings/forcing_file_list.txt

  # Run the second Python script
  python3 run_summa.py

  echo "summa run for $input3 complete, check output folder for output file"

done
