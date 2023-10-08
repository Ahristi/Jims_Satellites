#5104327867 AERO3760 Assignment 3 Position Determination
import ephem
import math
import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

# Create a GNSS observer
observer = ephem.Observer()

# Set the observer's location (latitude, longitude, and elevation)
observer.lat = -33  # Replace with the observer's latitude in radians
observer.lon = 151  # Replace with the observer's longitude in radians
observer.elev = 10  # Replace with the observer's elevation in meters

# Initialize the GNSS satellites (e.g., GPS, GLONASS, Galileo, etc.)
# You can specify the constellation you want to use.
# For example, for GPS:
gnss_constellation = ephem.readtle(
    'GPS',
    '1 28361U 04049A   14295.79166667  .00000086  00000-0  11606-4 0  9983',
    '2 28361  55.1506  86.2190 0004256  43.2610 316.8001  2.00560251  7529'
)

# Set the current time (UTC)
observer.date = ephem.now()

# Compute the satellite's position
gnss_constellation.compute(observer)

# Extract satellite information
satellite_azimuth = math.degrees(gnss_constellation.az)
satellite_elevation = math.degrees(gnss_constellation.alt)
satellite_range = gnss_constellation.range

print(f"Satellite Azimuth: {satellite_azimuth} degrees")
print(f"Satellite Elevation: {satellite_elevation} degrees")
print(f"Satellite Range: {satellite_range} meters")

# Simulated data (simplified)
# True initial state [x, y, vx, vy]
true_state = np.array([1000, 2000, 1, -2])

# Simulated star tracker measurements (simplified)
num_measurements = 10
star_positions = np.random.normal(true_state[:2], 10, (num_measurements, 2))  # Add noise to true positions

# Function to compute residuals (difference between predicted and measured star positions)
def residuals(state):
    predicted_positions = state[:2]  # Simplified motion model (constant velocity)
    return (predicted_positions - star_positions).flatten()

# Initial guess for the state
initial_guess = np.array([0, 0, 0, 0])

# Weight for each measurement (simplified)
measurement_weights = np.ones(num_measurements)

# Optimization using weighted non-linear least squares
result = least_squares(residuals, initial_guess, method='lm', ftol=1e-6, args=(), kwargs={})
estimated_state = result.x

# Plot the true and estimated positions
plt.figure(figsize=(10, 6))
plt.scatter(true_state[0], true_state[1], label="True Position", color="blue", marker="o")
plt.scatter(estimated_state[0], estimated_state[1], label="Estimated Position", color="red", marker="x")
plt.scatter(star_positions[:, 0], star_positions[:, 1], label="Star Tracker Measurements", color="green", marker="*")
plt.xlabel("X Position")
plt.ylabel("Y Position")
plt.legend()
plt.grid(True)
plt.title("True vs. Estimated Satellite Position")
plt.show()

print("True State:", true_state)
print("Estimated State:", estimated_state)
