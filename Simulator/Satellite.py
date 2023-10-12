"""
    Satellite

    Python implementation of a satellite object containing the 
    integration of different subsystems.

"""
from datetime import datetime



class Satellite:

    def __init__(self, TLEfile: str, name: str):
        """
            Satellite Object

            Inputs: 
            TLEfile - Name of the .txt file containing the TLE of the satellite
            name    - name of the satellite        
        """
        #Orbital Parameters
        self.inclination, self.rightAscension, self.eccentricity, self.argumentPerigee, self.meanAnomaly, self.meanMotion, self.tSinceVernal  = orbitfromTLE(TLEfile)

        #Coordinates
        self.states     =   []          #Each state (pos,vel) of the satellite in ECI frame
        self.ECEF       =   []          #Each position of the satellite in ECEF frame
        self.GLLH       =   []          #Each position of the satellite in GLLH frame
        self.times      =   []          #Each time during the simulation

        #Subsystems
        self.ADCS       =   None        #Satellite ADCS subsystem
        self.payload    =   None        #Satellite camera payload
        self.GNSS       =   None        #Satellite GNSS
        self.EPS        =   None        #Satellite EPS


        

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
    inclination                =   float(tleArray[1][2])
    rightAscension             =   float(tleArray[1][3])
    eccentricity               =   float(tleArray[1][4])*10**(-7)
    argumentPerigee            =   float(tleArray[1][5])
    meanAnomaly                =   float(tleArray[1][6])
    meanMotion                 =   float(tleArray[1][7])
    
    return inclination, rightAscension, eccentricity, argumentPerigee,meanAnomaly, meanMotion, tSinceVernal


if __name__ == "__main__":
    sat = Satellite("ISS.txt", "ISS")