import csv
import numpy as np
from orbitalTransforms import *


"""
    starTracker

    Object implementation of a star tracker in a satellite.
    Takes a text file of stars and their location in earth's frame.

    Currently only takes one star for proof of concept
"""
class starTracker:
    def __init__(self, starFile, accuracy):
        self.stars    =  self.initialiseStars(starFile)
        self.accuracy =  np.deg2rad(accuracy)

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
    
    def printStars(self):
        """
            Prints all the stars and their corresponding coordinates
            within the startracker config.       
        """
        for key, value in self.stars.items():
            print(f'{key}, {value}')

    def getActualReading(self,starName, satPos):
        """
            Gets the actual vector of a star from the satellite in ECI frame

            Inputs:

            starName    - the name of the star to use in the comparison. Assumes that it is in the text file.
            satpos      - numpy array of the satellite coordinate in ECI     
        """
        #Get the true vector between the satellite and the star:
        starFromSat = self.stars[starName] - satPos
        return starFromSat

    def getReading(self,starName, satPos, satAttitude):
        """
            Fakes the error involved in a star tracker reading by adding
            gaussian noise to the actual vector between the satellite and
            a star.

            Inputs:

            starName    - the name of the star to use in the comparison
            satpos      - numpy array of the satellite coordinate in ECI     
            satAttitude - the actual attitude of the satellite in heading, pitch, roll
        """
        #Get the true vector between the satellite and the star:
        starFromSat = self.stars[starName] - satPos

        #Calculate directional cosine matrix
        phi = satAttitude[0]
        theta = satAttitude[1]
        psi = satAttitude[2]

        C = np.array([[np.cos(theta)*np.cos(psi), np.cos(theta)*np.sin(psi), -np.sin(theta)],
                     [np.sin(phi)*np.sin(theta)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.sin(phi)*np.sin(theta)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.sin(phi)*np.cos(theta)],
                     [np.cos(phi)*np.sin(theta)*np.cos(psi) + np.sin(phi)*np.sin(psi), np.cos(phi)*np.sin(theta)*np.sin(psi)- np.sin(phi)*np.cos(psi), np.cos(phi)*np.cos(theta)]])
        
        #Get the coordinate of the star in the satellite body frame in polar form
        starInBody = C @ starFromSat
        inBodyPolar = CART2POLAR(starInBody)
        el = inBodyPolar[0]
        az = inBodyPolar[1]
        R  = inBodyPolar[2]
        
        #Add Gausian noise
        el = np.random.normal(el, self.accuracy)
        az = np.random.normal(az, self.accuracy)
        
        #Convert back to cartesian in the body frame of the satellite
        inBodyPolar =  np.array([el,az,R])
        starInBody  =  POLAR2CART(inBodyPolar)

        return starInBody

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

