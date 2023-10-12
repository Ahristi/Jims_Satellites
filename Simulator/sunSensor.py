"""
    Very basic skeleton of a sun sensor. 
    Doesn't change the reading at the moment, but still implementable for ADCS.
"""

import numpy as np
from orbitalTransforms import *


class sunSensor:
    def __init__(self, accuracy):
        self.accuracy = np.deg2rad(accuracy)
    
    def getActualReading(self):
        return np.array([1.496e+11,0,0])
    
    def getReading(self, attitude):
        C = directionalCosine(attitude[0],attitude[1], attitude[2])

        inBody = C @ self.getActualReading()
        inBodyPolar = CART2POLAR(inBody)
        el = inBodyPolar[0]
        az = inBodyPolar[1]
        R = inBodyPolar[2]

        #Add Gausian noise
        el = np.random.normal(el, self.accuracy)
        az = np.random.normal(az, self.accuracy)
        
        #Convert back to cartesian in the body frame of the satellite
        inBodyPolar =  np.array([el,az,R])
        sunInBody  =  POLAR2CART(inBodyPolar)

        return sunInBody
    

if __name__ == "__main__":
    sunSens = sunSensor(10)
    roll = np.deg2rad(45)
    pitch = np.deg2rad(45)
    yaw = np.deg2rad(45)


    attitude = np.array([roll,pitch,yaw])
    C = directionalCosine(attitude[0],attitude[1], attitude[2])
    print(C @ sunSens.getActualReading())
    print(sunSens.getReading(attitude))

    