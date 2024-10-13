#!/bin/bash

# Define the remote server and directory
REMOTE_USER="cdalden"
REMOTE_HOST="j-lundquist-3.ce.washington.edu"
REMOTE_DIR="/home/cdalden/summa_setup/twitter_api/plots/"

# Define the local directory
LOCAL_DIR="/Users/clintonalden/Documents/Research/summa_work/hrrr_initialized_realtime/plots/"

# Use scp to copy files from the remote server to the local directory
scp ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_DIR}*.png ${LOCAL_DIR}

# Loop through the inputs "density" and "temp"
for var in density temp
do
    # Execute the Python script with the current input
    echo $var | python tweet_plots.py
done

# move the plots to archive
mv ${LOCAL_DIR}*.png ${LOCAL_DIR}archive/