"""
    Satellite.py
    
    Object implentation of a satellite. Contains information about position and attitude
    
"""
import numpy as np
from datetime import datetime
from orbitalTransforms import PFF2ECEF
from Orbit import Orbit

u = 3.986004418*10**14  
R = 6378137

class Satellite:
    def __init__(self, TLEfile, epoch = None, ADCS = None):
        self.orbit = orbitfromTLE(TLEfile)
        self.params = self.orbit.eccentricity, np.radians(self.orbit.inclination), np.radians(self.orbit.argumentPerigee), np.radians(self.orbit.rightAscension), np.radians(self.orbit.meanAnomaly), self.orbit.meanMotion
        self.X, self.V = kepler_orbit(self.params, 0)
        self.ADCS = ADCS
        self.time = epoch
        self.attitude = np.array([0,0,0])
        self.states = [np.concatenate((self.X, self.V), axis=0)]
        self.times = []

    def setAttitude(self):
        """
            Assumes that the control is perfect and sets the satellite
            attitude to always be nadir pointing
        """
        x = self.states[-1][0]
        y = self.states[-1][1]
        z = self.states[-1][2]

        yaw = np.arctan2(y,x)
        pitch = np.arctan2(z, np.sqrt(x**2 +  y**2))
        roll = np.pi/2

        self.attitude = np.array([roll,pitch,yaw])










def calculateState(e,h,trueAnom):
    """
        Calculates position and velocity vector from TLE data

        Inputs:
        e: Eccentric anomaly
        h: Angular Momentum
        trueAnom: True Anomaly

        Returns:

        r: Position vector of the satellite
        v: Velocity vector of the satellite

    """
    r = (h**2 / u)*(1/(1+e*np.cos(trueAnom))) * np.array([np.cos(trueAnom), np.sin(trueAnom), 0])
    v = (u/h)* np.array([-np.sin(trueAnom), e + np.cos(trueAnom), 0])
    return r,v




def calculateEccentricAnomaly(M_e,e):
    """
        Calcualtes the eccentric anomaly of a satellite

        Inputs:
        M_e: Mean Anomaly of a satellite
        e:   eccentricity of satellite's orbit

        Returns:
        E: eccentric anomaly

    """
    #Assume that eccentric anomaly is equal to mean anomaly
    E = M_e
    while (np.abs(M_e + e*np.sin(E) - E) > 0.000000000001):
        E = M_e + e*np.sin(E)
    return E


# Define the equation of motion for a Keplerian orbit
def kepler_orbit(params, t):
    e, i, w, RAA, M0, meanMotion = params

    #Get number of revolutions per second
    meanMotionSeconds   =   meanMotion/(24*60*60)
    #Calculate period in seconds
    T   =   1/meanMotionSeconds
    n = 2*np.pi/T
    M_e = M0 + n*t

    #Calculate Eccentric Anomaly
    E   =   calculateEccentricAnomaly(M_e,e)
    #Calculate True Anomaly

    trueAnom    =   2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2))

    #Calculate perigee
    a   =   (((T*np.sqrt(u))/(2*np.pi))**(2/3))

    #Calculate the distance from the focus
    r   =   (a*(1-e**2)/(1+e*np.cos(trueAnom)))
    
    #Calculate the Angular Momentum
    h   =    np.sqrt(a*u*(1-e**2))
    
    #Calculate radius of perigee
    r_p =  (h**2/u)* (1/(1+e))
    
    #Calculate semi-major axis
    r_a =  2*a - r_p

    #Get initial state of the satellite in perifocal frame
    state = calculateState(e,h, trueAnom)
    #Convert perifocal to ECI
    X, V = PFF2ECEF(state[0],state[1],RAA,w,i)
    
    return X,V  # Return as a NumPy array




def orbitfromTLE(TLEfile):
    tleFile = open(TLEfile, "r")
    tleArray = []
    for i in range(2):
        tleArray.append(tleFile.readline().split())
    satelliteName              =   tleArray[0][1] + tleArray[1][1]

    #Line1
    EpochYearandJulianDay      =   tleArray[0][3]
    year = float(EpochYearandJulianDay[0:2])
    day = float(EpochYearandJulianDay[2:])
    daysJ2000 = year*362.25
    epoch = (daysJ2000 + day)*24*60*60 #Time passed since J2000 in seconds
    #Find time elapsed at vernal equinox
    vernal = datetime(2023,3,21,21,25) - datetime(2000,1,1,12,0)
    vernalSeconds = vernal.total_seconds()
    tSinceVernal = epoch - vernalSeconds
    #Line2
    inclination                =   float(tleArray[1][2])
    rightAscension             =   float(tleArray[1][3])
    eccentricity               =   float(tleArray[1][4])*10**(-7)
    argumentPerigee            =   float(tleArray[1][5])
    meanAnomaly                =   float(tleArray[1][6])
    meanMotion                 =   float(tleArray[1][7])
    
    return Orbit(inclination, rightAscension, eccentricity, argumentPerigee,meanAnomaly, meanMotion, tSinceVernal)
