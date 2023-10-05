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
R = 6378137

class magnetometer:
    def __init__(self, epochTime):
        self.epochTime = epochTime
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

        Br, Btheta, Bphi = ppigrf.igrf_gc(r, theta, phi, date) # returns radial, south, east
        bFieldLLH  = np.array([Btheta, Bphi, Br - R])
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

        C = np.array([np.cos(theta)*np.cos(psi), np.cos(theta)*np.sin(psi), -np.sin(theta)],
                     [np.sin(phi)*np.sin(theta)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.sin(phi)*np.sin(theta)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.sin(phi)*np.cos(theta)],
                     [np.cos(phi)*np.sin(theta)*np.cos(psi) + np.sin(phi)*np.sin(psi), np.cos(phi)*np.sin(theta)*np.sin(psi)- np.sin(phi)*np.cos(psi), np.cos(phi)*np.cos(theta)])
        

        bFieldECI  = self.getActualField(date, satPos)
        #Get the magnetic field in the body frame
        bFieldBody = C @ bFieldECI
        bFieldBodyPolar = CART2POLAR(bFieldBody)

        el = bFieldBodyPolar[0]
        az = bFieldBodyPolar[1]
        R  = bFieldBodyPolar[2]

        #Convert back to cartesian in the body frame of the satellite
        bFieldBodyPolar =  np.array([el,az,R])
        bFieldBody  =  POLAR2CART(bFieldBodyPolar)

        #Convert to ECI
        bFieldReading     =  C @ bFieldBody

        return bFieldReading