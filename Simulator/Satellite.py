"""
    Satellite

    Python implementation of a satellite object containing the 
    integration of different subsystems.

    NOTE:   The positions and attitudes used are the ACTUAL ones which we use for simulation and to add noise.

            The estimated positions, attitudes etc are the ones that are in the subsystems,
            Ie if you wanted the estimated attitude, then you would use the attitude from the 
            ADCS subsystem.

"""
from datetime import datetime

from orbitalTransforms import *


class Satellite:

    def __init__(self, TLEfile: str, name: str):
        """
            Satellite Object

            Inputs: 
            TLEfile - Name of the .txt file containing the TLE of the satellite
            name    - name of the satellite        
        """
        #Orbital Parameters
        params = orbitfromTLE(TLEfile)
        self.inclination, self.rightAscension, self.eccentricity, self.argumentPerigee, self.meanAnomaly, self.meanMotion, self.tSinceVernal  = params
        
        #Coordinates
        X,V = keplerOrbit(params, 0)
        state0  = np.concatenate((X, V))

        self.states     = [state0]      #Each state (pos,vel) of the satellite in ECI frame
        self.ECEF       =   []          #Each position of the satellite in ECEF frame
        self.GLLH       =   []          #Each position of the satellite in GLLH frame
        self.times      =   []          #Each time during the simulation

        #Attitude
        self.attitudes  =   []          #Each satellite attitude in ECI frame.
        self.setDesiredAttitude()       #Set the initial attitude
        
        
        #Subsystems
        self.ADCS       =   None        #Satellite ADCS subsystem
        self.payload    =   None        #Satellite camera payload
        self.GNSS       =   None        #Satellite GNSS
        self.EPS        =   None        #Satellite EPS
        
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


    def tick(self):
        """
            Peforms the satellite OBC functions.
            Called every iteration of the simulation.       
        
        """

        #ADCS routines
        self.setDesiredAttitude()

        #Position routine

        #Payload routines
     

        #Power routines



    def setDesiredAttitude(self):
        """
            Sets the desired attitude of the satellite. 
            For the time being attitude change is assumed to be instant.
        """
        x = self.states[-1][0]
        y = self.states[-1][1]
        z = self.states[-1][2]

        yaw = np.arctan2(y,x)
        pitch = np.arctan2(z, np.sqrt(x**2 +  y**2))
        roll = 0
        self.attitudes.append(np.array([roll,pitch,yaw]))


def orbitfromTLE(TLEfile):
    """  
        Output orbital parameters from a TLE text file.    
    """
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
    inclination                =  np.deg2rad(float(tleArray[1][2]))
    rightAscension             =  np.deg2rad(float(tleArray[1][3]))
    eccentricity               =             float(tleArray[1][4])*10**(-7)
    argumentPerigee            =  np.deg2rad(float(tleArray[1][5]))
    meanAnomaly                =  np.deg2rad(float(tleArray[1][6]))
    meanMotion                 =  float(tleArray[1][7])
    
    return inclination, rightAscension, eccentricity, argumentPerigee,meanAnomaly, meanMotion, tSinceVernal


if __name__ == "__main__":
    sat = Satellite("ISS.txt", "ISS")