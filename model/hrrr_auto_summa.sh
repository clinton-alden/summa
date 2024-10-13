#!/bin/bash

# Remove all plots in the specified directory first so twitter api is clean
rm -f /home/cdalden/summa_setup/twitter_api/plots/*.png

# Source the conda.sh script to enable conda commands
source /home/cdalden/miniforge3/etc/profile.d/conda.sh

# Activate the first conda environment and run the first Python script
conda activate pysumma

# Prompt for inputs
# echo "Enter the desired SNOTEL site code (ie. 1107:WA): "
# read input1
input1=679:WA # paradise snotel

input2=2025

# Automatically generate the output file name based on the current date
current_date=$(date +'%d%b%Y')
input3="paradise_$current_date"

# Provide input to the first Python script
echo -e "$input1\n$input2\n$input3" | python3 hrrr_forcings.py

echo "forcing file created, now running summa"

# Replace the content of the text file with input3
echo "'$input3.nc'" > forcings/forcing_file_list.txt

# Run the second Python script
python3 hrrr_run_summa.py

echo "summa run complete, check output folder for output file"

conda deactivate

# activate new env and post tweets
conda activate new_tweets

# Loop through the inputs "density" and "temp"
for var in density temp
do
    # Execute the Python script with the current input
    echo $var | python /home/cdalden/summa_setup/twitter_api/tweet_plots.py
done

# move the plots to archive
mv /home/cdalden/summa_setup/twitter_api/plots/*.png /home/cdalden/summa_setup/twitter_api/plots/archive/

conda deactivate


