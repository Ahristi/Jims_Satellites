#Class to get sat position using GNSS


class GNSS:
    def __init__(self, sat, GNSS):
        """
            ADCS object which simulates attitude determination in the satellite
            Contains different sensor objects to determine the attitude using WNLLS

            Inputs:

            starTrackerFile     -   Name of the csv file containing the star tracker config
            sat                 -   Satellite object that the ADCS belongs to
        
        """
        

        self.satellite = sat

        self.estimatedPositions = []
        

