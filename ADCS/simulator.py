import numpy as np
from ADCS import ADCS
from starTracker import *
from magnetometer import *
from sunSensor import *
from satellite import *
from orbitalTransforms import *
u = 3.986004418*10**14  

class simulator:

    def __init__(self, satellites, groundstations = None, numSatsforLink = 10):
        self.satellites = satellites
        self.groundstations = groundstations
        self.estimatePositions = False
        self.numSatsforLink = numSatsforLink

    def simulate(self, t0, t_end, h):
        for sat in self.satellites:
            sat.times.append(t0)
        while t0 < t_end:
            for sat in self.satellites:
                #Use 4th order runge kuta to solve
                k1 = h * self.motionEquation(sat.states[-1])
                k2 = h * self.motionEquation(sat.states[-1] + k1/2)
                k3 = h * self.motionEquation(sat.states[-1] + k2/2)
                k4 = h * self.motionEquation(sat.states[-1] + k3)
                
                state_new = sat.states[-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
                
                sat.times.append(t0+h)
                
                sat.states.append(state_new)

                X_ECEF = ECI2ECEF([state_new[0], state_new[1], state_new[2]],t0+h)

                GLLH = ECEF2GLLH([X_ECEF[0],X_ECEF[1],X_ECEF[2]])
            t0 += h

    def motionEquation(self, state):
        X = state[0:3]
        a  = calculateAcceleration(X)
        dS = np.concatenate((state[3:6], a), axis=0)
        return dS

def calculateAcceleration(X):
    """
        Calculates position and velocity vector from TLE data

        Inputs:
        e: Eccentric anomaly
        h: Angular Momentum
        trueAnom: True Anomaly

        Returns:

        r: Position vector of the satellite
        v: Velocity vector of the satellite
    """
    r = np.linalg.norm(X)
    a = (-u/r**3) * X 
    return a


if __name__ == "__main__":
        #Make the star tracker
    cross = starTracker("star_config.csv", 0.011)
    
    #Make the magnetometer
    J2000 = datetime(2000, 1, 1, 12)
    mg = Magnetometer(J2000,0.1)

    #Make the sun sensor
    sunsens = sunSensor(10)
    adcs = ADCS(mg,cross,sunsens)
    sat = Satellite("ISS.txt", J2000, adcs)

    roll = np.deg2rad(45)
    pitch = np.deg2rad(45)
    yaw = np.deg2rad(45)
    sat.setAttitude(np.array([roll,pitch,yaw]))
    sat.ADCS.connectToSatellite(sat)
    #Begin simulation
    sim = simulator([sat])
    sim.simulate(0,5,1)

    
    



    
    

