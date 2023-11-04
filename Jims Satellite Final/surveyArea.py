import numpy as np
import cv2
import ee
import geemap
import os
from PIL import Image
import matplotlib.pyplot as plt

def surveyArea(town_name, mosaic, historic_mosaic, minlat, maxLat, minLong, maxLong):
    
    #Save all files within the current folder
    current_directory = os.getcwd()
    out_dir = os.path.join(current_directory)
    #Check Comleroy
    bounds = [minLong ,minlat, maxLong, maxLat]
    #Save mapped mosaic
    filename = os.path.join(out_dir, 'eeResults/' + town_name + '.tif')
    roi = ee.Geometry.Rectangle(bounds, None, False)
    image = mosaic.clip(roi).unmask()
    geemap.ee_export_image(
        image, filename=filename, scale=50, region=roi, file_per_band=False
    )

    #Save mapped mosaic with historical data
    filename = os.path.join(out_dir, 'eeResults/' + town_name + '_historic.tif')
    roi = ee.Geometry.Rectangle(bounds, None, False)
    image = historic_mosaic.clip(roi).unmask()
    geemap.ee_export_image(
        image, filename=filename, scale=50, region=roi, file_per_band=False
    )


    # Open the TIFF file
    tiff_file = 'eeResults/' + town_name + '.tif'
    image = Image.open(tiff_file)
    image = np.array(image)

    #Convert to jpg
    cv2.imwrite("eeResults/" + town_name + "comleroy.jpg", cv2.cvtColor(image, cv2.COLOR_RGB2BGR))
    # Apply filter
    lower1 = np.array([50, 0, 0])
    upper1 = np.array([90, 2, 2])
    mask = cv2.inRange(image, lower1, upper1)
    white = 0
    total = 0
    for row in mask:
        for pixel in row:
            if (pixel == 255):
                white+=1
            total+=1
    ratio =white/total * 100

    # Create a figure with two subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    # Display the first image on the left subplot
    ax1.imshow(image, extent = [bounds[0], bounds[2], bounds[1], bounds[3]])
    ax1.set_title(town_name)

    # Display the second image on the right subplot
    ax2.imshow(mask, extent = [bounds[0], bounds[2], bounds[1], bounds[3]])
    ax2.set_title(town_name + " Masked")

    
    plt.show()
    print(town_name + " danger percentage:",ratio)
