import rasterio

# Open the GeoTIFF file
tif_file_path = 'GEEBAM_v3p1_dNBR_Classes.tif'
with rasterio.open(tif_file_path) as src:
    # Read the desired band (e.g., band 1)
    band_number = 1  # Adjust as needed
    band = src.read(band_number)
    
    # Get the dimensions of the raster
    rows, cols = band.shape

    # Loop through the rows and columns and print each pixel's value
    for row in range(rows):
        for col in range(cols):
            pixel_value = band[row, col]
            if pixel_value != 253:
                print(f'Row: {row}, Col: {col}, Value: {pixel_value}')