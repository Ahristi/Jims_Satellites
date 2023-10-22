import csv
import numpy as np
from orbitalTransforms import *
from Satellite import Satellite
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyvista as pv
from pyvista import examples
from scipy.optimize import least_squares

"""
    Position Determination

    Determines the position for a given satellite dependant on the mode 
    it is in. RTK/Differential mode is higher precision but only available when in 
    line of sight of the ground station. Regula r GNSS determines the 
    position for the remainder of the satellite's orbit.
"""
# Defining Constants
# Speed of light in meters per second
speed_of_light = 299792458

base_station_vis = True

# Determine if Regular or Differential GPS (more accurate) is to be used at the current position
if base_station_vis == True:
    differential = True
else:
    differential = False

# BEGIN Position Determination for given location

#Location of the Constellation cubesats at a given point in time
constellation_eci_pos = [
    np.array([ -4170263.97479639,  4372069.56522638, -3434042.26023409]),  # Satellite 1
    np.array([ -4199286.99037704,  4262858.38391383, -3534784.56508713]),  # Satellite 2
    np.array([ -4225190.94341261,  4152242.66350547, -3634451.5991731 ]),  # Satellite 3
    np.array([ -4247942.36734899,  4040309.67761361, -3733013.11002903]),  # Satellite 4
    np.array([ -4267510.24683579,  3927147.62679821, -3830439.18457065]),  # Satellite 5
    np.array([ -4297149.25825629,  3797868.85754494, -3926700.25807045])   # Satellite 6
]

#Location of the Parkes observatory Ground Station 
groundstation_pos_LLH = np.array([-32.9986, 148.2621, 415])

# Approximate ECI coordinates for 4 GPS/GNSS satellites over NSW, Australia
gps_satellites_eci = [
    np.array([0, 20000000, 0]),  # GPS Satellite 1
    np.array([-20500000, 0, 0]),  # GPS Satellite 2
    np.array([0, 25000000, -13000000]),  # GPS Satellite 3
    np.array([-17500000, 0, -13000000])   # GPS Satellite 4
]

def plot_satellite_positions(constellation_eci_pos, groundstation_pos_ecef, gps_satellites_eci):
    # Create a 3D plot with a specified size
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Set the aspect ratio to be equal for X, Y, and Z
    ax.set_box_aspect([1, 1, 1])

    # Plot Earth first
    Earth_radius = 6371000  # meters
    coefs = (1, 1, 1)
    rx, ry, rz = [Earth_radius / np.sqrt(coef) for coef in coefs]
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = rx * np.outer(np.cos(u), np.sin(v))
    y = ry * np.outer(np.sin(u), np.sin(v))
    z = rz * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(x, y, z, rstride=4, cstride=4, color='g', alpha=0.5)  # Specify zorder

    # Plot the satellite positions
    for i, satellite_pos in enumerate(constellation_eci_pos):
        ax.scatter(satellite_pos[0], satellite_pos[1], satellite_pos[2], label=f'CubeSat {i+1}', c='r')

        # Draw a line from the ground station to each CubeSat
        xs = [groundstation_pos_ecef[0], satellite_pos[0]]
        ys = [groundstation_pos_ecef[1], satellite_pos[1]]
        zs = [groundstation_pos_ecef[2], satellite_pos[2]]
        ax.plot(xs, ys, zs, c='k', linestyle='--', alpha=0.5)

    # Plot the ground station
    ax.scatter(groundstation_pos_ecef[0], groundstation_pos_ecef[1], groundstation_pos_ecef[2], label='Ground Station', c='g', marker='o')

    # Plot GPS satellites
    for i, satellite_pos in enumerate(gps_satellites_eci):
        ax.scatter(satellite_pos[0], satellite_pos[1], satellite_pos[2], label=f'GPS Satellite {i+1}', c='b', marker='x')

        # Draw lines from GPS satellites to CubeSat 1
        xs = [satellite_pos[0], constellation_eci_pos[0][0]]
        ys = [satellite_pos[1], constellation_eci_pos[0][1]]
        zs = [satellite_pos[2], constellation_eci_pos[0][2]]
        ax.plot(xs, ys, zs, c='b', linestyle='--', alpha=0.5)

    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m')

    # Set axis limits to extend to 25,000 km
    ax.set_xlim(-25000e3, 25000e3)
    ax.set_ylim(-25000e3, 25000e3)
    ax.set_zlim(-25000e3, 25000e3)

    plt.title('Constellation GPS/GNSS Visualisation')
    plt.legend()
    plt.show()

t = 11692
groundstation_pos_ecef = GLLH2ECEF(groundstation_pos_LLH)
# groundstation_pos_eci = ECEF2ECI(groundstation_pos_ecef,t)

print(groundstation_pos_ecef)


def ECI_to_ECEF_array(eci_coords, t):
    ecef_coords = [ECI2ECEF(coord, t) for coord in eci_coords]
    return ecef_coords

# Transform the CubeSat and GPS satellite coordinates from ECI to ECEF
constellation_ecef_pos = ECI_to_ECEF_array(constellation_eci_pos, t)
gps_satellites_ecef = ECI_to_ECEF_array(gps_satellites_eci, t)


# Call the function to plot satellite positions in ECEF
plot_satellite_positions(constellation_ecef_pos, groundstation_pos_ecef, gps_satellites_ecef)

def calculate_distance(point1, point2):
    # Calculate the distance between two points
    distance = np.linalg.norm(point1 - point2)
    return distance

def calculate_time_from_distance(distance):
    # Calculate the time it takes for light to travel the given distance
    time = distance / speed_of_light
    return time

def apply_noise_to_time(time, noise_level):
    # Apply Gaussian noise to the time amount
    noise = np.random.normal(scale=noise_level)
    noisy_time = time + noise
    return noisy_time

def calculate_pseudo_ranges_and_apply_noise(cube_satellite_positions, gps_satellite_positions, noise_level):
    pseudo_ranges = []

    for cube_satellite in cube_satellite_positions:
        pseudo_range_to_cube_sat = []
        for gps_satellite in gps_satellite_positions:
            # Calculate the distance between the CubeSat and the GPS satellite
            distance = calculate_distance(cube_satellite, gps_satellite)
            # Convert distance to time and apply noise
            time = calculate_time_from_distance(distance)
            noisy_time = apply_noise_to_time(time, noise_level)
            pseudo_range_to_cube_sat.append(noisy_time)
        pseudo_ranges.append(pseudo_range_to_cube_sat)

    return pseudo_ranges

# Apply noise to pseudo-ranges
noise_level = 1e-10  # Adjust this as needed
pseudo_ranges = calculate_pseudo_ranges_and_apply_noise(constellation_ecef_pos, gps_satellites_ecef, noise_level)

# Print the noisy pseudo-range measurements
for i, pseudo_range_to_cube_sat in enumerate(pseudo_ranges):
    for j, noisy_time in enumerate(pseudo_range_to_cube_sat):
        print(f'Noisy Pseudo-range from CubeSat {i+1} to GPS Satellite {j+1}: {noisy_time} seconds')


print(pseudo_ranges)


def trilateration_equations(cube_sat_position, gps_satellite_positions, pseudo_ranges, speed_of_light):
    equations = []
    for i in range(len(gps_satellite_positions)):
        distance = pseudo_ranges[i] * speed_of_light
        delta_x = cube_sat_position[0] - gps_satellite_positions[i][0]
        delta_y = cube_sat_position[1] - gps_satellite_positions[i][1]
        delta_z = cube_sat_position[2] - gps_satellite_positions[i][2]
        equations.append((delta_x ** 2 + delta_y ** 2 + delta_z ** 2) - (distance ** 2))
    return equations

# Perform trilateration and store the results in cube_sat_positions
cube_sat_positions = []
for pseudo_range_to_cube_sat in pseudo_ranges:
    initial_guess = [0, 0, 0]  # Initial guess for CubeSat position
    result = least_squares(trilateration_equations, initial_guess, args=(gps_satellites_ecef, pseudo_range_to_cube_sat, speed_of_light))
    cube_sat_positions.append(result.x)


# Calculate differences between calculated and actual positions
differences = [np.array(actual_pos) - np.array(calc_pos) for actual_pos, calc_pos in zip(constellation_ecef_pos, cube_sat_positions)]


def plot_positional_errors(differences):
    # Calculate the positional errors for each CubeSat
    x_errors = [diff[0] for diff in differences]
    y_errors = [diff[1] for diff in differences]
    z_errors = [diff[2] for diff in differences]
    total_errors = [np.linalg.norm(diff) for diff in differences]

    # Calculate the average total error
    avg_total_error = np.mean(total_errors)

    # CubeSat labels
    cube_sat_labels = [f'CubeSat {i+1}' for i in range(6)]

    # Create an array of colors for different error types
    colors = ['r', 'g', 'b', 'c']

    # Create x values for the CubeSats
    x_values = np.arange(len(cube_sat_labels))

    # Plot x, y, z, and total errors for each CubeSat as separate points
    plt.figure(figsize=(10, 6))
    plt.scatter(x_values, x_errors, color=colors[0], label='X Error', marker='o')
    plt.scatter(x_values, y_errors, color=colors[1], label='Y Error', marker='o')
    plt.scatter(x_values, z_errors, color=colors[2], label='Z Error', marker='o')

    # Connect the total error points with a line
    plt.plot(x_values, total_errors, color=colors[3], label='Total Error', marker='o', linestyle='-')

    # Add a red line for the average total error
    plt.axhline(y=avg_total_error, color='red', linestyle='--', label='Average Total Error')

    # Add a dotted line at 0 for theoretical value
    plt.axhline(y=0, color='k', linestyle='--', label='Theoretical Value')

    # Add labels and legend
    plt.xlabel('CubeSat')
    plt.ylabel('ECEF Positional Error (meters)')
    plt.xticks(x_values, cube_sat_labels)  # Set x-axis labels
    plt.title('ECEF Positional Error for CubeSat Constellation')
    plt.legend()

    # Show the plot
    plt.grid()
    plt.show()

# Plot positional error
plot_positional_errors(differences)


# ********************************************************************
# Differential GPS
# Step 1: Calculate Pseudo-Range Measurements for the Ground Station
def calculate_pseudo_ranges_for_ground_station(groundstation_pos_ecef, gps_satellites_ecef, noise_level):
    pseudo_ranges = []

    for gps_satellite in gps_satellites_ecef:
        # Calculate the distance between the ground station and each GPS satellite
        distance = calculate_distance(groundstation_pos_ecef, gps_satellite)
        # Convert distance to time and apply noise
        time = calculate_time_from_distance(distance)
        noisy_time = apply_noise_to_time(time, noise_level)
        pseudo_ranges.append(noisy_time)

    return pseudo_ranges

# Apply noise to pseudo-ranges for the ground station
pseudo_ranges_ground_station = calculate_pseudo_ranges_for_ground_station(groundstation_pos_ecef, gps_satellites_ecef, noise_level)

# Print the noisy pseudo-range measurements for the ground station
for i, noisy_time in enumerate(pseudo_ranges_ground_station):
    print(f'Noisy Pseudo-range from Ground Station to GPS Satellite {i+1}: {noisy_time} seconds')

# Step 2: Determine the Position of the Ground Station
def ground_station_position_equations(ground_station_position, gps_satellite_positions, pseudo_ranges, speed_of_light):
    equations = []
    for i in range(len(gps_satellite_positions)):
        distance = pseudo_ranges[i] * speed_of_light
        delta_x = ground_station_position[0] - gps_satellite_positions[i][0]
        delta_y = ground_station_position[1] - gps_satellite_positions[i][1]
        delta_z = ground_station_position[2] - gps_satellite_positions[i][2]
        equations.append((delta_x ** 2 + delta_y ** 2 + delta_z ** 2) - (distance ** 2))
    return equations

# Use the least squares solver to estimate the ground station position
initial_guess_ground_station = [0, 0, 0]  # Initial guess for ground station position
result_ground_station = least_squares(ground_station_position_equations, initial_guess_ground_station, args=(gps_satellites_ecef, pseudo_ranges_ground_station, speed_of_light))
estimated_ground_station_position = result_ground_station.x

# Calculate the differences between the estimated ground station position and the known precise ground station position
ground_station_difference = np.array(groundstation_pos_ecef) - estimated_ground_station_position

print("Estimated Ground Station Position:", estimated_ground_station_position)
print("Ground Station Position Difference:", ground_station_difference)

# Step 3: Improve the Position of the CubeSats
# Subtract the weighted ground station position difference from the CubeSats' ECEF positions to improve their accuracy
# Function to calculate weighted improvement for CubeSat positions
def weighted_improvement(cube_sat_positions, ground_station_difference):
    improved_positions = []

    for i, cube_sat_pos in enumerate(cube_sat_positions):
        diff = np.array(constellation_ecef_pos[i]) - cube_sat_pos
        weight_x = 1 #/ (abs(diff[0]))  # Weight for X dimension
        weight_y = 1 #/ (abs(diff[1]))   # Weight for Y dimension
        weight_z = 1 #/ (abs(diff[2]))   # Weight for Z dimension

        improved_x = cube_sat_pos[0] + weight_x * ground_station_difference[0]
        improved_y = cube_sat_pos[1] + weight_y * ground_station_difference[1]
        improved_z = cube_sat_pos[2] + weight_z * ground_station_difference[2]

        improved_positions.append([improved_x, improved_y, improved_z])

    return improved_positions

# Calculate the weighted improvement of CubeSat positions
improved_cube_sat_positions = weighted_improvement(cube_sat_positions, ground_station_difference)

# Calculate differences between the improved positions and the known positions of the CubeSats
improved_differences = [np.array(actual_pos) - np.array(improved_pos) for actual_pos, improved_pos in zip(constellation_ecef_pos, improved_cube_sat_positions)]

print("Improved Differences: ",improved_differences)

def plot_improved_positional_errors(differences):
    # Calculate the positional errors for each CubeSat
    x_errors = [diff[0] for diff in differences]
    y_errors = [diff[1] for diff in differences]
    z_errors = [diff[2] for diff in differences]
    total_errors = [np.linalg.norm(diff) for diff in differences]

    # Calculate the average total error
    avg_total_error = np.mean(total_errors)

    # CubeSat labels
    cube_sat_labels = [f'CubeSat {i+1}' for i in range(6)]

    # Create an array of colors for different error types
    colors = ['r', 'g', 'b', 'c']

    # Create x values for the CubeSats
    x_values = np.arange(len(cube_sat_labels))

    # Plot x, y, z, and total errors for each CubeSat as separate points
    plt.figure(figsize=(10, 6))
    plt.scatter(x_values, x_errors, color=colors[0], label='Differential X Error', marker='o')
    plt.scatter(x_values, y_errors, color=colors[1], label='Differential Y Error', marker='o')
    plt.scatter(x_values, z_errors, color=colors[2], label='Differential Z Error', marker='o')

    # Connect the total error points with a line
    plt.plot(x_values, total_errors, color=colors[3], label='Total Error', marker='o', linestyle='-')

    # Add a red line for the average total error
    plt.axhline(y=avg_total_error, color='red', linestyle='--', label='Average Total Error')

    # Add a dotted line at 0 for theoretical value
    plt.axhline(y=0, color='k', linestyle='--', label='Theoretical Value')

    # Add labels and legend
    plt.xlabel('CubeSat')
    plt.ylabel('Differential ECEF Positional Error (meters)')
    plt.xticks(x_values, cube_sat_labels)  # Set x-axis labels
    plt.title('Differential ECEF Positional Error for CubeSat Constellation')
    plt.legend()

    # Show the plot
    plt.grid()
    plt.show()

# Plotting Over Time

# Plot improved positional error
plot_improved_positional_errors(improved_differences)


# Define the time steps over which you want to perform trilateration
time_steps = range(0, 1000, 10)  # Define your own time steps here

# List to store the error deviation over time
error_deviation_over_time = []

# Function to calculate the trilateration error for a given time step
def trilateration_error_at_time_step(time_step):
    # Calculate pseudo-ranges for the specific time step (you may need to adapt this part of your code)
    pseudo_ranges = calculate_pseudo_ranges_and_apply_noise(constellation_ecef_pos, gps_satellites_ecef, noise_level)

    # Perform trilateration and store the results in cube_sat_positions
    cube_sat_positions = []
    for pseudo_range_to_cube_sat in pseudo_ranges:
        initial_guess = [0, 0, 0]  # Initial guess for CubeSat position
        result = least_squares(trilateration_equations, initial_guess, args=(gps_satellites_ecef, pseudo_range_to_cube_sat, speed_of_light))
        cube_sat_positions.append(result.x)

    # Calculate differences between calculated and actual positions
    differences = [np.array(actual_pos) - np.array(calc_pos) for actual_pos, calc_pos in zip(constellation_ecef_pos, cube_sat_positions)]

    # Calculate the average total error for this time step
    avg_total_error = np.mean([np.linalg.norm(diff) for diff in differences])

    return avg_total_error

# Loop through time steps and calculate error deviation
for time_step in time_steps:
    avg_error = trilateration_error_at_time_step(time_step)
    error_deviation_over_time.append(avg_error)

# Calculate the average of all time step errors
average_error_overall = np.mean(error_deviation_over_time)

# Plot the error deviation over time
plt.figure(figsize=(12, 6))  # Increase the figure size
plt.plot(time_steps, error_deviation_over_time)
plt.axhline(y=average_error_overall, color='red', linestyle='--', label=f'Overall Average Error: {average_error_overall:.3f} meters')  # Add a red dashed line for overall average error
plt.xlabel('Time Step')
plt.ylabel('Total Jim-1 Error Deviation (meters)')
plt.title('Position Error Over Time for Differential GNSS')

# Set the y-axis limits from 0 to 0.1 meters
plt.ylim(0, 0.1)

plt.grid()
plt.legend()  # Show the legend with the overall average
plt.show()