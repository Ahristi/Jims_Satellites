"""
    Satellite

    Python implementation of a satellite object containing the 
    integration of different subsystems and the motion of a satellite from a simulation.

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
from GNSS import *



BATTERY_CAPACITY = 22#Whr
NSW_BOUNDING = [[-37.10333642191903, 140.89503281528468], [-28.07557676995773, 154.61690531304407]] #Bounding box for NSW 

#Satellite States
SAFE    =   0
IMAGING =   1
READOUT =   2



class Satellite:

    def __init__(self, TLEfile: str, name: str):
        """
            Satellite Object

            Inputs: 
            TLEfile - Name of the .txt file containing the TLE of the satellite
            name    - name of the satellite        
        """
        self.name = name

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

        #Misc
        self.eclipses   =   [False]                  #Boolean containing whether or not the satellite was under an eclipse
        self.state      =   SAFE                     #Current mode of operation
        self.imaging    =   [SAFE]                       #Array containing boolean values if the satellite was imaging or not
        self.time       =   datetime(2000, 1, 1, 12) #Used for magnetometer. Just hard coding this for now but will change later.
        self.mappingErrors = [] 
        self.gnssConstellation  = []                 #Array of GNSS Satellites in the simulation


        #Attitude
        self.attitudes  = []            #Each satellite attitude in ECI frame.
        self.attitude   = []            #Current attitude in ECI frame
        self.setDesiredAttitude()       #Set the initial attitude
        
        
        #Subsystems
        self.ADCS       =   ADCS("star_config.csv", self)     #Satellite ADCS subsystem
        self.payload    =   Payload(self)                     #Satellite camera payload
        self.GNSS       =   GNSS(self)                        #Satellite GNSS (Not yet implemented)
        self.EPS        =   EPS(BATTERY_CAPACITY,self)        #Satellite EPS


    def tick(self):
        """
            Peforms the satellite OBC functions.
            Called every iteration of the simulation.
        """

        #Get the amount of time passed since last tick 
        h = self.times[-1] - self.times[-2]

        # #ADCS routines
        # self.setDesiredAttitude()
        # self.ADCS.determineAttitude()
        # self.ADCS.starTracker1.getReading(self.states[-1][0:3], self.attitude)

        #Position routine
        self.GNSS.estimatePosition()

        # #Payload routines
        self.updateState()
        # self.payload.obtainPointing()

        # #Power routines
        # self.eclipses.append(self.checkEclipse()) #Check if eclipse is occuring
        # self.imaging.append(self.state==IMAGING)
        # self.EPS.calculateCharge(h)  #Calculate amount of charge in the battery


    def setDesiredAttitude(self):
        """
            Sets the desired attitude of the satellite. 
            For the time being attitude change is assumed to be instant.

            #TODO: Add control from magnetometer and reaction wheels (not urgent)
        """
        x = self.states[-1][0]
        y = self.states[-1][1]
        z = self.states[-1][2]

        yaw = np.arctan2(y,x)
        pitch = np.arctan2(z, np.sqrt(x**2 +  y**2))
        roll = 0
        currentAttitude = np.array([roll,pitch,yaw])
        self.attitude = currentAttitude
        self.attitudes.append(currentAttitude)
        


    def checkEclipse(self):  
        """
            Checks if an eclipse is occurring based on the satellite position and sun position.
            Assumes caller is updating the sun position.

            #TODO: Proper calculation of eclipse.
        """
        sunFromEarth = np.linalg.norm(self.sunPos)
        satECI = self.states[-1][0:3]
        sunFromSat   = np.linalg.norm(self.sunPos - satECI)
        if (sunFromEarth < sunFromSat):
            return True
        
        else:
            return False
    def updateState(self):
        """
            Checks if we are in the NSW bounding box 
            which means that we should be imaging

            TODO: Turn the bounding box into a more complex polygon        
        """
        lat = self.GLLH[-1][0]
        long = self.GLLH[-1][1]
        newState = SAFE
        if (lat > NSW_BOUNDING[0][0] and lat < NSW_BOUNDING[1][0]):
            if(long > NSW_BOUNDING[0][1] and long < NSW_BOUNDING[1][1]):
                newState = IMAGING

        self.state = newState

    def connectGNSS(self, gnssConstellation, groundStations):
        """
            Adds the GNSS satellites from a simulation into the attributes
            of the satellite to be used in GNSS position determination.     
        """
        self.gnssConstellation = gnssConstellation
        self.GNSS.groundStations = groundStations
        return





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
    sat = Satellite("ISS.txt", "ISS")