"""
    EPS - Electrical Power System
   
"""


#Constants for the power used by different states as determined by the power budget
#TODO: Calculate this from the subsystem objects.
IMAGING_POWER = 21.98
READOUT_POWER = 16.48
SAFE_POWER    = 3.28

SAFE    =   0
IMAGING =   1
READOUT =   2



class EPS:
    def __init__(self, capacity, sat):
        
        #Pointer to the satellite it is apart of
        self.sat = sat


        #Current parameters
        self.charge = capacity
        self.maxCharge = capacity

        #Datalogging
        self.charges = [self.charge] #Numpy array containing the battery charge at each time step

        return
    


    def calculateCharge(self, h):
        """
            Updates the state of charge of the EPS battery pack

            Inputs:    
            h   - time step
        """
        eclipsing = self.sat.eclipses[-1] #Check if eclipse is occuring:

        #Just do a dodgy state machine for now where the satellite images during eclipse
        #And is in safe during non eclipse
        #I know this makes no sense but I just want to do it to format the plots for 
        #When the integration is done

        #TODO: Calculate the correct amount of power being used
        powerDraw = self.calculatePowerDraw()

        if (eclipsing):
            self.charge = max(self.charge - powerDraw/3600,0)
        else:
            self.charge = min(self.charge + (self.calculateSolarPower()- powerDraw)/3600, self.maxCharge)

        #Log the current charge
        self.charges.append(self.charge)
        

    def calculateSolarPower(self):
        """
            Calculates the amount of energy harvested by the solar panels
            based on the satellite state

            Inputs:
            sat - satellite object which the EPS belongs to 
            
            Returns:

            power - power harvested in Watts      
        """

        #TODO: Actually do the calcuation based on satellite attitude
        return 8.3
    
    def calculatePowerDraw(self):
        """
            Calculate the amount of power draw based on the current state.
            For the moment this is hard coded but it can be made more detailed later        
        """
        if (self.sat.state == SAFE):
            return SAFE_POWER
        elif (self.sat.state == IMAGING):
            return IMAGING_POWER
        elif (self.sat.state == READOUT):
            return READOUT_POWER
        else:
            return 0
        

