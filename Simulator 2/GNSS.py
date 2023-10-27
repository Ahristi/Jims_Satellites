"""
    GNSS

    Implementation of a GNSS receiver in the satellite.


"""



class GNSS:
    def __init__(self, sat):
        self.satellite = sat
        self.positionEstimates = []  #Array of the estimated positions as determined using RTK method
        self.position          = None  #Calculated position
        self.groundStations    = []  #Array of groundstation objects

    def estimatePosition(self):
        gnssConstellation = self.satellite.gnssConstellation #List of GNSS satellite objects
        actualPosition = self.satellite.states[-1][0:3]      #Actual position of the satellite in ECI frame
        
        #Just returning the actual position for now 
        self.position = actualPosition
        return actualPosition
