import ee
import numpy as np
service_account = "jims-948@jimssatellites.iam.gserviceaccount.com"
key_path = 'jimssatellites-25b0240d5cf2.json'
credentials = ee.ServiceAccountCredentials(service_account, key_path)
ee.Initialize(credentials=credentials)
import geemap
import folium
# Set some points of interest
# Mackay, Qld
lat = -21.1434
lon = 149.1868
# Get URL for our image
Data_Set_Clearing = ee.ImageCollection("WRI/GFW/FORMA/raw_output_firms")
#Click on the down arrows to see the dataset.
Data_Set_Clearing


# Create an interactive map
Map = geemap.Map(center=[lat, lon], zoom=7)
# Define the dataset
dataset = ee.ImageCollection('WRI/GFW/FORMA/raw_output_firms')\
.filter(ee.Filter.date('2018-08-01', '2018-08-15'))
# Select the 'nday' band
percentageOfClearing = dataset.select('nday')
# Visualization parameters
visParams = {
'min': 0.0,
'max': 0.01
}
# Add layer to map
Map.addLayer(percentageOfClearing, visParams, 'Percentage of clearing')
# Show the map
Map