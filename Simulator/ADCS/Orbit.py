class Orbit:
        """orbit class contains information about a two line element set."""
        def __init__(self, inclination, rightAscension, eccentricity, argumentPerigee, meanAnomaly, meanMotion, epoch):

            self.epoch      =   epoch
            #Line2
            self.inclination                =   inclination       
            self.rightAscension             =   rightAscension    
            self.eccentricity               =   eccentricity      
            self.argumentPerigee            =   argumentPerigee   
            self.meanAnomaly                =   meanAnomaly       
            self.meanMotion                 =   meanMotion        
        def getTimeSinceVernal(self):
            """
                Uses the TLE epoch year and Julian day to calculate
                the amount of time that has passed since the vernal equinox
                
                Inputs:
                N/A

                Returns:
                Amount of time passed since last vernal equinox in seconds.
            """
            #Find time elapsed at vernal equinox
            return self.epoch
        def printOrbit(self):
            print("Inclination:     " + str(self.inclination ))
            print("Eccentricity:    " + str(self.eccentricity))
            print("Right Ascension: " + str(self.rightAscension))
            print("Perigee:         " + str(self.argumentPerigee))
            print("Mean Anomaly:    " + str(self.meanAnomaly))
            print("Mean Motion:     " + str(self.meanMotion))
