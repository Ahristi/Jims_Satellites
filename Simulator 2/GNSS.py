"""
    GNSS

    Implementation of a GNSS receiver in the satellite.


"""
#Satellite States
SAFE    =   0
IMAGING =   1
READOUT =   2



class GNSS:
    def __init__(self, sat):
        self.satellite = sat
        self.positionEstimates = []  #Array of the estimated positions as determined using RTK method
        self.position          = None  #Calculated position
        self.groundStations    = []  #Array of groundstation objects

    def estimatePosition(self):
        gnssConstellation = self.satellite.gnssConstellation        #List of GNSS satellite objects
        actualPosition = self.satellite.states[-1][0:3]             #Actual position of the satellite in ECI frame
        t =  self.satellite.times[-1] + self.satellite.tSinceVernal #Time since vernal equinox for ECEF conversion
        state = self.satellite.state                                #The state of the satellite 



        # Check the state of the satellite
        if state == "IMAGING":
            # Code to execute when the state is IMAGING (high precision differential GNSS)
            # For example, you can calculate and store the estimated position
            self.position = actualPosition
            self.positionEstimates.append(self.position)
        elif state == "SAFE":
            # Code to execute when the state is SAFE (lower precision regular GNSS)
            # You can do something else here
            pass  # Placeholder for SAFE state actions


        #Just returning the actual position for now 
        self.position = actualPosition                  #Return the calculated GNSS Position
        self.positionEstimates.append(self.position)    #Append the calculated GNSS Position
