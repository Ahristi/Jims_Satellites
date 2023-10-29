import cv2
import numpy as np





from PIL import Image

# Open the TIFF file
tiff_file = "landsat.tif"
image = Image.open(tiff_file)

# You can now manipulate the image or extract information from it.
# For example, you can display the image:
image.show()

# To access pixel data, you can convert the image to a NumPy array
import numpy as np
image_data = np.array(image)

# To access image information, such as dimensions and mode (e.g., "RGB" or "L" for grayscale), you can use the following attributes:
width, height = image.size
image_mode = image.mode

print(f"Image dimensions: {width} x {height}")
print(f"Image mode: {image_mode}")








# Load the image
image = cv2.imread('Results/onelayer1.jpg')

# Apply filter
lower1 = np.array([0, 0, 66])
upper1 = np.array([5, 5, 97])
mask = cv2.inRange(image, lower1, upper1)

 
# Display the image with detected blobs
cv2.imshow('Image with Blobs', mask)
cv2.waitKey(0)
cv2.destroyAllWindows()


