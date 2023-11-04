import csv
import numpy as np
from orbitalTransforms import *
from skyfield.api import Star, load
from skyfield.data import hipparcos
from skyfield.framelib import itrs

"""
    starTracker

    Object implementation of a star tracker in a satellite.
    Takes a text file of stars and their location in earth's frame.

    Currently only takes one star for proof of concept
"""

STAR_MAGNITUDE_THRESHOLD= 0.8

class starTracker:
    def __init__(self, accuracy):
        self.accuracy =  np.deg2rad(accuracy)

        #Load bright stars from skyfield

        with load.open(hipparcos.URL) as f:
            df = hipparcos.load_dataframe(f)
        df = df[df['ra_degrees'].notnull()]
        df = df[df['magnitude'] <= STAR_MAGNITUDE_THRESHOLD]
        self.df = df
        self.planets = load('de421.bsp')
        self.earth = self.planets['earth']
        self.ts = load.timescale()



    def initialiseStars(self, starFile):
        """
            Returns a dictionary of stars from a config file

            Inputs: 
            stars - csv file in the format star_name, star_x, star_y, star_z

            Outputs: A dictionary where keys are the star names and the values are the coordinates of the star
            
            Note: coordinates for stars are in Earth centered earth fixed frame.
        """
        starData = {}
        # Open and read the CSV file
        with open(starFile, mode='r', newline='') as file:
            reader = csv.reader(file)
            # Iterate through each row and populate the dictionary
            for row in reader:
                if len(row) >= 4:  # Ensure there are at least 4 columns in each row
                    star_name = row[0]
                    star_x = float(row[1])
                    star_y = float(row[2]) 
                    star_z = float(row[3])
        
                    
                    # Create a dictionary entry with star name as key and coordinates as value
                    starData[star_name] = np.array([star_x, star_y, star_z])
                    
        return starData
    

    def getActualReading(self, satPos, time = None):
        """
            Uses the starfield API to get the actual positions of known stars in the ECI frame
            with respect to the satellite

            Inputs:

            sat  - the satellite position in ECI frame.
            time - datetime object of the time of observation (Currently not used)  
        
        """

        #Load bright stars from starfield (this gives around 15 observations )
        bright_stars = Star.from_dataframe(self.df)
        t = self.ts.utc(2000, 1, 1)
        astrometric = self.earth.at(t).observe(bright_stars)
        ra, dec, distance = astrometric.radec()
        x, y, z = astrometric.frame_xyz(itrs).au
        
        #Build the observation Matrix
        obs = np.zeros((len(x), 3))
        obs[:,0] = x * 1.496e+11 #Converting AU to metres
        obs[:,1] = y * 1.496e+11
        obs[:,2] = z * 1.496e+11
        #Build a matrix where each row is the position
        pos = np.zeros((len(x), 3))
        pos[:,] = satPos
        obsFromSat = obs - pos
        return obsFromSat

    def getReading(self, satPos, satAttitude):
        """
            Fakes the error involved in a star tracker reading by adding
            gaussian noise to the actual vector between the satellite and
            a star.

            Inputs:

            starName    - the name of the star to use in the comparison
            satpos      - numpy array of the satellite coordinate in ECI     
            satAttitude - the actual attitude of the satellite in heading, pitch, roll
        """  
        #Get the true vector between the satellite and the star in ECI frame:
        starsFromSat = self.getActualReading(satPos)

        #Calculate directional cosine matrix
        phi = satAttitude[0]
        theta = satAttitude[1]
        psi = satAttitude[2]
        C = np.array([[np.cos(theta)*np.cos(psi), np.cos(theta)*np.sin(psi), -np.sin(theta)],
                     [np.sin(phi)*np.sin(theta)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.sin(phi)*np.sin(theta)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.sin(phi)*np.cos(theta)],
                     [np.cos(phi)*np.sin(theta)*np.cos(psi) + np.sin(phi)*np.sin(psi), np.cos(phi)*np.sin(theta)*np.sin(psi)- np.sin(phi)*np.cos(psi), np.cos(phi)*np.cos(theta)]])
        

        #Start building the array of readings. 
        readings = np.zeros(np.shape(starsFromSat))
        for i in range(len(starsFromSat)):
            #Get star in the body frame
            starInBody = C @ starsFromSat[i]
            #Get polar coordinate in the satellite body frame
            el,az,r = CART2POLAR(starInBody)
            #Add Gausian noise based on the star tracker's accuracy
            el = np.random.normal(el, self.accuracy)
            az = np.random.normal(az, self.accuracy)
            inBodyPolar =  np.array([el,az,r]) #This is with noise
            #Convert back to cartesian coordinates in the body frame
            starInBody  =  POLAR2CART(inBodyPolar)
            readings[i,:] = starInBody
     
        return readings


 

def CART2POLAR(X):
    """
        CART2POLAR

        Inputs:
        X - Position of an object in ECEF frame

        Outputs:
        3x1 vector containing elevation, azimuth and distance

    """
    x     =   X[0]
    y     =   X[1]
    z     =   X[2]
    R     =   np.sqrt(x**2 + y**2 + z**2)
    el    =   np.arcsin(z/R)
    az    =   np.arctan2(y,x)

    return np.array([el,az, R])

def POLAR2CART(X):
    """
            POLAR2CART

            Inputs:
            X - Position of an object in ECEF frame

            Outputs:
            3x1 vector containing elevation, azimuth andx distance

    """
    el    =    X[0]
    az    =    X[1]
    r     =    X[2]

    x     =    r * np.cos(az)*np.cos(el)
    y     =    r * np.sin(az)*np.cos(el)
    z     =    r * np.sin(el)
    return np.array([x,y,z])

