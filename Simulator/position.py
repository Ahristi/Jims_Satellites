import csv
import numpy as np
from orbitalTransforms import *
from Satellite import Satellite

"""
    Position Determination

    Determines the position for a given satellite dependant on the mode 
    it is in. RTK mode is higher precision but only available when in 
    line of sight of the ground station. Regular GNSS determines the 
    position for the remainder of the satellite's orbit.
"""

def getPositionRTK(satellite, gnss_satellites, noise_std_dev=0.0):
    """
    Calculate the position of the satellite relative to GNSS satellites and apply noise.

    Parameters:
    - satellite: The satellite for which you want to determine the position.
    - gnss_satellites: A list of GNSS satellite objects, each with known positions.
    - noise_std_dev: Standard deviation of measurement noise (default is 0.0 for no noise).

    Returns:
    - Estimated position of the satellite relative to GNSS satellites with noise.
    """

    # Calculate pseudorange measurements to GNSS satellites
    pseudorange_measurements = []
    for gnss_satellite in gnss_satellites:
        # Calculate pseudorange for each GNSS satellite
        pseudorange = calculatePseudorange(satellite, gnss_satellite)
        # Apply measurement noise
        if noise_std_dev > 0.0:
            pseudorange += np.random.normal(0, noise_std_dev)
        pseudorange_measurements.append(pseudorange)

    # Perform trilateration to estimate the relative position
    estimated_position = trilateration(gnss_satellites, pseudorange_measurements)

    return estimated_position

def calculatePseudorange(satellite, gnss_satellite):
    """
    Calculate the pseudorange between the satellite and a GNSS satellite.

    Parameters:
    - satellite: The satellite for which you want to calculate the pseudorange.
    - gnss_satellite: A GNSS satellite with a known position.

    Returns:
    - Pseudorange measurement from the satellite to the GNSS satellite.
    """
    # Implement the pseudorange calculation here
    # You'll need to account for the satellite's position, GNSS satellite's position, and signal propagation delay.
    # For example, you can use the formula:
    # pseudorange = c * (t_receiver - t_transmitter)
    # where c is the speed of light, t_receiver is the time when the signal is received,
    # and t_transmitter is the time when the signal was transmitted.

def trilateration(gnss_satellites, pseudorange_measurements):
    """
    Perform trilateration to estimate the position of the satellite relative to GNSS satellites.

    Parameters:
    - gnss_satellites: A list of GNSS satellite objects, each with known positions.
    - pseudorange_measurements: A list of pseudorange measurements to GNSS satellites.

    Returns:
    - Estimated relative position of the satellite.
    """
    # Implement trilateration to estimate the relative position
    # You'll need to use the known positions of GNSS satellites and pseudorange measurements.

# Usage example
if __name__ == "__main__":
    # Create a satellite and GNSS satellites with known positions
    satellite = Satellite("SatelliteTLE.txt", "SatelliteName")
    gnss_satellites = [GNSSSatellite("GNSS1", position=[x1, y1, z1]),
                      GNSSSatellite("GNSS2", position=[x2, y2, z2]),
                      # Add more GNSS satellites with known positions
                     ]

    # Call getPositionRTK to estimate the position of the satellite relative to GNSS satellites
    estimated_position = getPositionRTK(satellite, gnss_satellites, noise_std_dev=1.0)

    print("Estimated Position:", estimated_position)