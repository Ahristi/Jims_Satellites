"""
    Simulator

    Simulator class that propogates multiple satellite orbits.
    Performs satellite and ground station operations

"""
import numpy as np
from Satellite import *
from orbitalTransforms import *


class simulator:

    def __init__(self, satellites, groundstations):

        self.satellites = satellites            #Array of satellite objects to be simulated
        self.groundstations = groundstations    #Array of groundstations to be simulated


    def simulate(self, t0, t_end, h, f):
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
                

                X_ECEF = ECI2ECEF([state_new[0], state_new[1], state_new[2]],t0+h) 
                GLLH = ECEF2GLLH([X_ECEF[0],X_ECEF[1],X_ECEF[2]])
            t0 += h
