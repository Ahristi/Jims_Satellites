"""
    Simulator

    Simulator class that propogates multiple satellite orbits.
    Performs satellite and ground station operations

"""
import numpy as np
from Satellite import *
from orbitalTransforms import *


class Simulator:

    def __init__(self, satellites, groundstations):

        self.satellites = satellites            #Array of satellite objects to be simulated
        self.groundstations = groundstations    #Array of groundstations to be simulated


    def simulate(self, t0, t_end, h, f):
        """
            Uses a 4th order numerical propogator to simulate the satellite orbit.

            Inputs:
            t0      -   initial time
            t_end   -   end time
            h       -   time step
            f       -   orbital function for simulation
        """
        for sat in self.satellites:
            sat.times.append(t0)
        while t0 < t_end:
            for sat in self.satellites:
                k1 = h * f(sat.states[-1])
                k2 = h * f(sat.states[-1] + k1/2)
                k3 = h * f(sat.states[-1] + k2/2)
                k4 = h * f(sat.states[-1] + k3)

                state_new = sat.states[-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
                sat.states.append(state_new)
                sat.times.append(t0+h)

                currentECEF = ECI2ECEF([state_new[0], state_new[1], state_new[2]],t0+h)
                currentGLLH = ECEF2GLLH([currentECEF[0],currentECEF[1],currentECEF[2]])

                sat.ECEF.append(currentECEF)
                sat.GLLH.append(currentGLLH)
         
            t0 += h

    def showOrbit():
        """
            Plots the satellite orbit using pyvista

            Assumes that the orbit has already been propogated using simulate()

            NOTE: This only plots the first satellite orbit at the moment        
        """



if __name__ == "__main__":
    sat = Satellite("ISS.txt", "ISS")
    sim = Simulator([sat], [])
    sim.simulate(0,24*60*60, 10, motionEquation)
    print(sat.states)

