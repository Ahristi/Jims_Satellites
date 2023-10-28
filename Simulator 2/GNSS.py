

"""
    GNSS

    Implementation of a GNSS receiver in the satellite.
    Utilises a Hybrid Differential GNSS strategy.
    If the Ground Station is in view and mapping is enabled,
    Differential GNSS is used, otherwise regular GNSS.


"""

import numpy as np
from orbitalTransforms import *
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Satellite States
SAFE    =   0       #Used when not in LOS and regular GNSS is requires.
IMAGING =   1       #Used when mapping and high precision GNSS required.
READOUT =   2



class GNSS:
    def __init__(self, sat):
        self.satellite = sat
        self.positionEstimates = []  #Array of the estimated positions as determined using RTK method
        self.position          = []  #Calculated position
        self.groundStations    = []  #Array of groundstation objects
        self.noise_level = 1e-10  # Noise associated with GNSS Satellite

    def estimatePosition(self):
        gnssConstellation = self.satellite.gnssConstellation        #List of GNSS satellite objects
        actualPosition = self.satellite.states[-1][0:3]             #Actual position of the satellite in ECI frame
        t =  self.satellite.times[-1] + self.satellite.tSinceVernal #Time since vernal equinox for ECEF conversion
        state = self.satellite.state                                #The state of the satellite 
        speed_of_light = 299792458                                  #Speed of light used for trilateration

        #**********************************************************
        #Frame conversions to ECEF

        #convert ECI Actual to ECEF Actual
        actualECEF = ECI2ECEF(actualPosition, t)

        gnssConstellationECEF = []

        #convert GNSS ECI to ECEF
        for gnssECI in gnssConstellation:
            gnssECEF = ECI2ECEF(gnssECI, t)
            gnssConstellationECEF.append(gnssECEF)

        #convert Groundstation GLLH to ECEF
        gsECEF = GLLH2ECEF(self.groundStations)

        #**********************************************************

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

        def calculate_pseudo_ranges_and_apply_noise(cube_satellite_position, gps_satellite_positions, noise_level):

            pseudo_range_to_cube_sat = []
            for gps_satellite in gps_satellite_positions:
                # Calculate the distance between the CubeSat and the GPS satellite
                distance = calculate_distance(cube_satellite_position, gps_satellite)
                # Convert distance to time and apply noise
                time = calculate_time_from_distance(distance)
                noisy_time = apply_noise_to_time(time, noise_level)
                pseudo_range_to_cube_sat.append(noisy_time)

            return pseudo_range_to_cube_sat
        

        def trilateration_equations(cube_sat_position, gps_satellite_positions, pseudo_ranges, speed_of_light):
            equations = []
            for i in range(len(gps_satellite_positions)):
                distance = pseudo_ranges[i] * speed_of_light
                delta_x = cube_sat_position[0] - gps_satellite_positions[i][0]
                delta_y = cube_sat_position[1] - gps_satellite_positions[i][1]
                delta_z = cube_sat_position[2] - gps_satellite_positions[i][2]
                equations.append((delta_x ** 2 + delta_y ** 2 + delta_z ** 2) - (distance ** 2))
            return equations




        # # Determine GNSS Differential Correction based on the state
        # # Check the state of the satellite
        # if state == "IMAGING":
        #     # Code to execute when the state is IMAGING (high precision differential GNSS)
        #     # For example, you can calculate and store the estimated position
        #     self.position = actualPosition
        #     self.positionEstimates.append(self.position)
        # elif state == "SAFE":
        #     # Code to execute when the state is SAFE (lower precision regular GNSS)
        #     # You can do something else here
        #     pass  # Placeholder for SAFE state actions






        #***********************************************************************
        #Postioning Main Code

        # Generating Pseudo Ranges
        pseudo_ranges = calculate_pseudo_ranges_and_apply_noise(actualECEF, gnssConstellationECEF, self.noise_level)
        print(pseudo_ranges)

        # Perform trilateration and store the results in cube_sat_positions
        initial_guess = [0, 0, 0]  # Initial guess for CubeSat position
        result = least_squares(trilateration_equations, initial_guess, args=(gnssConstellationECEF, pseudo_ranges, speed_of_light))
        cube_sat_position = result.x
        
        print("Calculated Position Found: ", cube_sat_position)
 

        #Just returning the actual position for now 
        self.position = cube_sat_position                  #Return the calculated GNSS Position
        self.positionEstimates.append(self.position)    #Append the calculated GNSS Position
        #***********************************************************************

if __name__ == "__main__":
    print("Running GNSS test")

    class Satellite:
        def __init__(self):
            self.gnssConstellation = []  # List of GNSS satellite objects
            self.states = []  # List of satellite states in ECI frame
            self.state = "SAFE"  # The state of the satellite, initialized to "SAFE"
            self.times = []  # List of times
            self.tSinceVernal = 0  # Time since vernal equinox

    # Create a Satellite object and set its attributes with fake data
    satellite = Satellite()

    # Set the GNSS coordinates for the satellite (example data)
    satellite.gnssConstellation = [
    np.array([0, 20000000, 0]),  # GNSS Satellite 1
    np.array([-20500000, 0, 0]),  # GNSS Satellite 2
    np.array([0, 25000000, -13000000]),  # GNSS Satellite 3
    np.array([-17500000, 0, -13000000])  # GNSS Satellite 4
    ]

    # Set the Satellite Times
    satellite.times = [11692]

    # Set the state of the satellite (assuming a single state for all GNSS satellites)
    satellite.state = "IMAGING"

    satellite.states = [np.array([ -4170263.97479639,  4372069.56522638, -3434042.26023409])]

    # Create an instance of the GNSS class and call the estimatePosition method
    gnss_receiver = GNSS(satellite)
    gnss_receiver.groundStations = (-32.9986, 148.2621, 415)
    gnss_receiver.estimatePosition()

    # # Print or use the results
    # print("Actual Position:", GNSS.position)
    # print("Position Estimates:", GNSS.positionEstimates)

    print("Success")
