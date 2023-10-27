"""
    GNSS

    Python implementation of a  GNSS satellite object. Removes some computationally heavy routines of the satellite object

    NOTE:   The positions and attitudes used are the ACTUAL ones which we use for simulation and to add noise.

            The estimated positions, attitudes etc are the ones that are in the subsystems,
            Ie if you wanted the estimated attitude, then you would use the attitude from the 
            ADCS subsystem.

"""
from datetime import datetime
from orbitalTransforms import *
from EPS import *
from ADCS import *
from Payload import *



BATTERY_CAPACITY = 22#Whr
NSW_BOUNDING = [[-37.10333642191903, 140.89503281528468], [-28.07557676995773, 154.61690531304407]] #Bounding box for NSW 

#Satellite States
SAFE    =   0
IMAGING =   1
READOUT =   2



class navSatellite:
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
        self.sunPos     =   None        #The actual position of the sun in ECI frame.

    def tick(self):
        #Get the amount of time passed since last tick 
        h = self.times[-1] - self.times[-2]


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

    tSinceVernal = (daysJ2000 + day)*24*60*60 #Time passed since J2000 in seconds

    #Hard coding this for now to match Mathu's code
    date1 = datetime(year=2023, month=3, day=20, hour=21, minute=24, second=0) # Date and time of most recent vernal equinox
    date2 = datetime(year=2023, month=10, day=4, hour=12, minute=14, second=0) # Date and time of TLE epoch
    tSinceVernal = (date2 - date1).total_seconds()

    #Line2
    inclination                =  np.deg2rad(float(tleArray[1][2]))
    rightAscension             =  np.deg2rad(float(tleArray[1][3]))
    eccentricity               =             float(tleArray[1][4])*10**(-7)
    argumentPerigee            =  np.deg2rad(float(tleArray[1][5]))
    meanAnomaly                =  np.deg2rad(float(tleArray[1][6]))
    meanMotion                 =  float(tleArray[1][7])
    
    return inclination, rightAscension, eccentricity, argumentPerigee,meanAnomaly, meanMotion, tSinceVernal


if __name__ == "__main__":
    sat = navSatellite("ISS.txt", "ISS")