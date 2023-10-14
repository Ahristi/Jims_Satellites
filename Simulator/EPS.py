"""
    EPS - Electrical Power System
   
"""


#Constants for the power used by different states as determined by the power budget
#TODO: Calculate this from the subsystem objects.
IMAGING_POWER = 21.98
READOUT_POWER = 16.48
SAFE_POWER    = 3.28


class EPS:
    def __init__(self, capacity):
        
        #Current parameters
        self.charge = capacity
        self.maxCharge = capacity

        #Datalogging
        self.charges = [self.charge] #Numpy array containing the battery charge at each time step

        return
    


    def calculateCharge(self, sat, h):
        """
            Updates the state of charge of the EPS battery pack

            Inputs:
            sat - satellite object which the EPS belongs to        
            h   - time step
        """
        eclipsing = sat.eclipses[-1] #Check if eclipse is occuring:

        #Just do a dodgy state machine for now where the satellite images during eclipse
        #And is in safe during non eclipse
        #I know this makes no sense but I just want to do it to format the plots for 
        #When the integration is done

        #TODO: Calculate the correct amount of power being used
        if (eclipsing):
            self.charge = max(self.charge - IMAGING_POWER/3600,0)
        else:
            self.charge = min(self.charge + (self.calculateSolarPower(sat)- SAFE_POWER )/3600, self.maxCharge)

        #Log the current charge
        self.charges.append(self.charge)
        

    def calculateSolarPower(self, sat):
        """
            Calculates the amount of energy harvested by the solar panels
            based on the satellite state

            Inputs:
            sat - satellite object which the EPS belongs to 
            
            Returns:

            power - power harvested in Watts      
        """

        #TODO: Actually do the calcuation based on satellite attitude
        return 28

