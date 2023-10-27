import numpy as np
from orbitalTransforms import *

GAUSSIAN_NOISE = 0
ELEVATION_NOISE = 1


class groundStation:
    def __init__(self, lat,long, el):
        """
            groundStation class

            Requires:
            lat: Latitude coordinate of the groundstation
            long: longitude coordinate of the groundstation
            el: elevation of the groundstation above sea level.
        """

        self.GLLH = [lat,long,el]               #GLLH coordinate of the ground station
    
        #Groundstation parameters (change these if necessary)
        self.nq = 0.2 #https://cddis.nasa.gov/lw17/docs/papers/posters/01-Makram_etal_pap.pdf
        self.Et = 13.5e-3 #https://link.springer.com/article/10.1007/s12567-019-00247-x
        self.y = 532*10**(-9)
        self.h = 6.62607015*10**(-34)
        self.c = 300000000
        self.nt = 0.9
        self.nr = 0.9
        self.G = 5.16e10
        self.sigma = 46e6 #https://cddis.nasa.gov/lw17/docs/papers/posters/01-Makram_etal_pap.pdf
        self.Ta = 1 # Assume Clear
        self.Ar = 0.7 #https://www.researchgate.net/figure/The-second-generation-of-Helwan-SLR-station_fig1_279313965



    def addObservation(self, t, satellite):
        """
            TODO:
            
            Adds observation of a satellite to the observations dictionary

            Inputs;
            t         - time of observation
            satellite - satellite object being simulated        
        """


        return
        
    
    def calculatePhotons(self, R):
        """
            calculatePhotons

            Calculates the number of photons on the receiver for a measurement

            inputs: 
            R: the distance between the ground station and the satellite

            outputs:
            Number of photons received by the groundstations       
        
        """
        Npe = self.nq*(self.Et * self.y/(self.h*self.c)) * self.nt * self.G * self.sigma * (1/(4*np.pi*R**2))**2 * self.Ar * self.nr * self.Ta**2
        return Npe


