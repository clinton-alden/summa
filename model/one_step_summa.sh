#!/bin/bash

# Prompt for inputs
echo "Enter the desired SNOTEL site code (ie. 1107:WA): "
read input1
echo "Enter the water year: "
read input2
echo "Enter the output file name (ie. buck_WY16): "
read input3

# Provide input to the first Python script
echo -e "$input1\n$input2\n$input3" | python3 snotel_to_pysumma_shortWY.py

echo "forcing file created, now running summa"

# Replace the content of the text file with input3
echo "'$input3.nc'" > forcings/forcing_file_list.txt

# Run the second Python script
# mkdir ./output
# mkdir ./output/plots
python3 run_summa.py

echo "summa run complete, check output folder for output file"
