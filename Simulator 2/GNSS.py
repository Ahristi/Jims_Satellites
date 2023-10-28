import numpy as np
from orbitalTransforms import *

"""
    GNSS

    Implementation of a GNSS receiver in the satellite.
    Utilises a Hybrid Differential GNSS strategy.
    If the Ground Station is in view and mapping is enabled,
    Differential GNSS is used, otherwise regular GNSS.


"""
#Satellite States
SAFE    =   0       #Used when not in LOS and regular GNSS is requires.
IMAGING =   1       #Used when mapping and high precision GNSS required.
READOUT =   2



class GNSS:
    def __init__(self, sat):
        self.satellite = sat
        self.positionEstimates = []  #Array of the estimated positions as determined using RTK method
        self.position          = None  #Calculated position
        self.groundStations    = []  #Array of groundstation objects

    def estimatePosition(self):
        gnssConstellation = self.satellite.gnssConstellation        #List of GNSS satellite objects
        actualPosition = self.satellite.states[-1][0:3]             #Actual position of the satellite in ECI frame
        t =  self.satellite.times[-1] + self.satellite.tSinceVernal #Time since vernal equinox for ECEF conversion
        state = self.satellite.state                                #The state of the satellite 
        speed_of_light = 299792458                                  #Speed of light used for trilateration

        #convert ECI Actual to ECEF Actual
        actualECEF = ECI2ECEF(actualPosition)

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
        





        # Determine GNSS Differential Correction based on the state
        # Check the state of the satellite
        if state == "IMAGING":
            # Code to execute when the state is IMAGING (high precision differential GNSS)
            # For example, you can calculate and store the estimated position
            self.position = actualPosition
            self.positionEstimates.append(self.position)
        elif state == "SAFE":
            # Code to execute when the state is SAFE (lower precision regular GNSS)
            # You can do something else here
            pass  # Placeholder for SAFE state actions







        #Postioning Main Code
        


        #Just returning the actual position for now 
        self.position = actualPosition                  #Return the calculated GNSS Position
        self.positionEstimates.append(self.position)    #Append the calculated GNSS Position


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

    # Set the state of the satellite (assuming a single state for all GNSS satellites)
    satellite.state = "IMAGING"


    # Call the estimatePosition method to test it
    gnss.estimatePosition()

    # Print or use the results
    print("Estimated Position:", gnss.position)
    print("Position Estimates:", gnss.positionEstimates)
