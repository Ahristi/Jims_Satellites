import cv2
import numpy as np





from PIL import Image

# Open the TIFF file
tiff_file = "eeResults/sydney.tif"
image = Image.open(tiff_file)


image = np.array(image)
cv2.imwrite("eeResults/sydney.jpg", cv2.cvtColor(image, cv2.COLOR_RGB2BGR))

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
print(ratio)

# Display the image with detected blobs
cv2.imshow('Image with Blobs', mask)
cv2.waitKey(0)
cv2.destroyAllWindows()


