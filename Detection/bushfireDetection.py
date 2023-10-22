import cv2
import numpy as np

# Load the image
image = cv2.imread('Fires.jpg')

# Convert the image to the HSV color space
hsv_image = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)

#Filter out values outside the range of the "inferno" colormap
lower_color = np.array([40, 11, 84])  
upper_color = np.array([255, 255, 255])  

# Create a binary mask using inRange
mask = cv2.inRange(hsv_image, lower_color, upper_color)

# Apply the mask to the original image
result = cv2.bitwise_and(image, image, mask=mask)


# Create a detector with the parameters you want (you can adjust these)
params = cv2.SimpleBlobDetector_Params()
 
# Filter by Area (you can adjust these values)
params.filterByArea = True
params.minArea = 4  # Minimum blob area
params.maxArea = 200000  # Maximum blob area

# Filter by Circularity (you can adjust these values)
params.filterByCircularity = False
params.minCircularity = 0.2  # Minimum circularity

# Filter by Convexity (you can adjust these values)
params.filterByConvexity = False
params.minConvexity = 0.87  # Minimum convexity

# Filter by Inertia (you can adjust these values)
params.filterByInertia = False
params.minInertiaRatio = 0.01  # Minimum inertia ratio

# Create the detector
detector = cv2.SimpleBlobDetector_create(params)

# Detect blobs
keypoints = detector.detect(result)

# Draw detected blobs on the image
result_with_keypoints = cv2.drawKeypoints(result, keypoints, np.array([]), (0, 0, 255),
                                         cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)

# Display the result with detected blobs
cv2.imshow('Result with Blobs', result_with_keypoints)
cv2.waitKey(0)
cv2.destroyAllWindows()