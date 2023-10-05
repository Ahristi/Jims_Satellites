"""
    ADCS (Attitude Determination and Control System)
    
    Object oriented implementation of cubeADCS platform in software simulation.
     
    Note: Control is not implemented, the name is just used for convenience. 
"""

from satellite import Satellite
from orbitalTransforms import *
from magnetometer import Magnetometer
from starTracker import starTracker


class ADCS:
    def __init__(self, satellite, magnetometer = None, startracker = None):
        self.magnetometer = magnetometer
        self.starTracker = startracker
        self.satellite = satellite

    def determineAttitude(self):
        """
            determineAttitude:

            Determines the satellite attitude based on the sensor readings
            by applying weighted non linear least squares.        
        """
        return
    

if __name__ == "__main__":
    adcs = ADCS()
    print(adcs.starTracker)

        
    

