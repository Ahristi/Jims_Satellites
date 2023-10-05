"""
    magnetometer.py


    Implementation of an ADCS magnetometer as an object.
    Uses the IGRF model to determine the vector of earth's 
    magnetic field at the satellite's position.

"""
import ppigrf
import numpy as np
from datetime import datetime
from orbitalTransforms import ECI2ECEF
from orbitalTransforms import ECEF2GLLH
from orbitalTransforms import GLLH2ECEF 
from orbitalTransforms import ECEF2ECI
from orbitalTransforms import CART2POLAR
from orbitalTransforms import POLAR2CART
from satellite import Satellite
R = 6378137


np.random.seed(4)
class Magnetometer:
    def __init__(self, epochTime, accuracy):
        self.epochTime = epochTime
        self.accuracy  = accuracy
        return
    
    def getActualField(self, date, satPos):
        """
            getActualField

            Returns the magnetic field at the satellite in ECI frame.

            date        -   datetime of the observation
            satPos      -   position of the satellite ECI
        """

        dT = self.epochTime - date

        # Get the number of seconds as an integer
        seconds_difference = int(dT.total_seconds())

        satPosECEF = ECI2ECEF(satPos, seconds_difference)
        satPosLLH  = ECEF2GLLH(satPosECEF)

        theta   =   90-satPosLLH[0]
        phi     =   satPos[1]
        r       =   (R + satPosLLH[2])/1000

        Br, Btheta, Bphi = ppigrf.igrf_gc(r, theta, phi, date) # returns radial, south, east (R, lat long)
        
        bFieldLLH  = np.array([Btheta[0], Bphi[0], Br[0] - R])
        bFieldECEF = GLLH2ECEF(bFieldLLH)
        
        bFieldECI  = ECEF2ECI(bFieldECEF, seconds_difference)
        
        return bFieldECI
    
    def getReading(self, date, satPos, satAttitude):
        """
            getReading
            Returns a reading of the magnetic field at the satellite in ECI frame accounting for noise.

            date        -   datetime of the observation
            satPos      -   position of the satellite ECI
            satAttitude    -   true attitude of the satellite  

        """
        #Calculate directional cosine matrix
        psi = satAttitude[0]
        theta = satAttitude[1]
        phi = satAttitude[2]

        C = np.array([[np.cos(theta)*np.cos(psi), np.cos(theta)*np.sin(psi), -np.sin(theta)],
                     [np.sin(phi)*np.sin(theta)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.sin(phi)*np.sin(theta)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.sin(phi)*np.cos(theta)],
                     [np.cos(phi)*np.sin(theta)*np.cos(psi) + np.sin(phi)*np.sin(psi), np.cos(phi)*np.sin(theta)*np.sin(psi)- np.sin(phi)*np.cos(psi), np.cos(phi)*np.cos(theta)]])
        

        bFieldECI  = self.getActualField(date, satPos)
        
        #Get the magnetic field in the body frame
        bFieldBody = C @ bFieldECI
        
        bFieldBodyPolar = CART2POLAR(bFieldBody)
        
        #Add Gausian noise
        el = np.random.normal(bFieldBodyPolar[0], self.accuracy)
        az = np.random.normal(bFieldBodyPolar[1], self.accuracy)
        R  = bFieldBodyPolar[2]

        #Convert back to cartesian in the body frame of the satellite
        bFieldBodyPolar =  np.array([el,az,R])
        
        bFieldBody  =  POLAR2CART(bFieldBodyPolar)
        
        #Convert to ECI
        bFieldReading     =  C @ bFieldBody

        return bFieldReading

if __name__ == "__main__":

    satellite = Satellite("ISS.txt")
    J2000 = datetime(2000, 1, 1, 12)
    mg = Magnetometer(J2000, 0.5)

    pitch = np.arctan2(satellite.X[1], satellite.X[2])
    yaw   = np.arctan2(satellite.X[0], satellite.X[1])
    roll = np.arctan2(np.sin(pitch)*np.cos(yaw), np.cos(pitch)*np.cos(yaw))
    satAttitude = np.array([yaw,pitch,roll])
    print(mg.getReading(J2000, satellite.X, satAttitude))
