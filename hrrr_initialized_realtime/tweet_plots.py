# %%
var = input("density or temp")
if var == "temp":
    var_long = "temperature"
elif var == "density":
    var_long = "density"

import tweepy
# the below import is for calculating date. Not needed for you but I needed it.
from datetime import date, datetime, timedelta
import shutil, pathlib, os

# take these from developer.twitter.com
ACCESS_KEY = "4544611537-MH3ykKPtEflWNJ3bVr4JQxug7pdppoD9lIfFCFR" # access token
ACCESS_SECRET = "lP4Oe0TksQw3DDsYU9GUBPxMSn6Yicd8gpoP2o1gHmQ4i" # access token secret
CONSUMER_KEY = "B5n3mp5OP3cgn7tNgtQK1XjVQ" # api key
CONSUMER_SECRET = "B1NfF3u2AV0an6CPGgkyEXNmeD6IFkUbVO0vhPjJeJDazBhcf2" # api secret
BEARER_TOKEN = "AAAAAAAAAAAAAAAAAAAAAHxywQEAAAAAvuJUZ9w3QM1GE2BBG%2FBaSAdO3zQ%3DTyAm5E12w6s8Wh0ucUiirjXxdorMdGcQu827btvPiA9qaGIMzH"

# Authenticate to Twitter
auth = tweepy.OAuthHandler(CONSUMER_KEY, CONSUMER_SECRET)
auth.set_access_token(
    ACCESS_KEY,
    ACCESS_SECRET,
)
# this is the syntax for twitter API 2.0. It uses the client credentials that we created
newapi = tweepy.Client(
    bearer_token=BEARER_TOKEN,
    access_token=ACCESS_KEY,
    access_token_secret=ACCESS_SECRET,
    consumer_key=CONSUMER_KEY,
    consumer_secret=CONSUMER_SECRET,
)

# Create API object using the old twitter APIv1.1
api = tweepy.API(auth)

# adding the tweet content in a multiline string. The {mydays} part is updated dynamically as the number of days from 6th Nov, 2023
# Get the current date
today = datetime.today()
# Format the date as 'DDMMM'
date_str = today.strftime('%d%b%Y')
# Calculate the date 2 days later
two_days_later = today + timedelta(days=2)
# Format the date as 'DDMMM'
date_str_two_days_later = two_days_later.strftime('%d%b%Y')

sampletweet = f'SUMMA-HRRR Snow Depth and {var_long.capitalize()} 48hr Forecast initialized {date_str} for Paradise, WA.\n\nThis tool is automated and purely experimental and should not be relied on for decision making. Please refer to nwac.us for the avalanche forecast.'

# Define the directory to search
search_dir = os.path.join(os.getcwd(), 'plots')

# Find the first PNG image in the ./plots directory
img = None
for images in os.listdir(search_dir):
    # check if the image ends with "density.png" and take the first image that you find
    if images.endswith(f"{var}.png"):
        img = images
        break

if img is None:
    raise FileNotFoundError(f"No image ending with '{var}.png' found in the ./plots directory")

# Upload the media using the old API
media = api.media_upload(os.path.join(search_dir, img))

# Create the tweet using the new API, mentioning the image uploaded via the old API
post_result = newapi.create_tweet(text=sampletweet, media_ids=[media.media_id])

# Print the response from the API
print(post_result)


