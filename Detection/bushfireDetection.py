import cv2
import numpy as np

# Load the image
image = cv2.imread('Fires.jpg')

# Convert the image to the HSV color space
hsv_image = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)

# Filter out values outside the range of the "inferno" colormap
lower_color = np.array([40, 11, 84])
upper_color = np.array([255, 255, 255])

# Create a binary mask using inRange
mask = cv2.inRange(hsv_image, lower_color, upper_color)

# Apply the mask to the original image
result = cv2.bitwise_and(image, image, mask=mask)

# Convert the result to grayscale
gray_result = cv2.cvtColor(result, cv2.COLOR_BGR2GRAY)

# Apply a binary threshold
_, binary_result = cv2.threshold(gray_result, 1, 255, cv2.THRESH_BINARY)

# Find contours in the binary result
contours, _ = cv2.findContours(binary_result, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

# Draw the detected blobs on the original image
image_with_blobs = image.copy()
cv2.drawContours(image_with_blobs, contours, -1, (0, 0, 255), 2)  # Draw contours in red

# Display the image with detected blobs
cv2.imshow('Image with Blobs', image_with_blobs)
cv2.waitKey(0)
cv2.destroyAllWindows()
