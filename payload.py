"""
    Payload

    Python implementation of Imaging system for Jim's Satellites Cubesat

    Creates bounding box of visibility based on estimated attitude from ADCS
    and estimated position of satellite

    NOTE:   Currently not considering roll in creation of bounding box and am
            assuming that swath is perpindicular to the velocity vectorwhich is
            technically incorrect as the velocity is not always perpindicular to
            position however, with a circular orbit, its close enough
"""

import numpy as np
import math
from Simulator.orbitalTransforms import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d



class payload:
    def __init__(self, attitude, posvel, time):
        """
            Payload Object
            
            Inputs:
            attitude:   [yaw, pitch, roll] estimated attitude of satellite in ECI frame
            posvel:     [Xp, Yp, Zp, Xv, Yv, Zv] estimated position/velocity of satellite in ECI frame
            time:       time since vernal in seconds
        """
        # self.satellite = None ***************************

        self.direction = np.array([0,0,0], dtype='int64')          # vector decribing satellites pointing direction in ECO
        self.observation = np.array([0,0,0], dtype='int64')        # vector describing ECI position satellite is pointing at
        self.bound_1_ECI = np.array([0,0,0], dtype='int64')        # Left swath bound in ECI frame
        self.bound_2_ECI = np.array([0,0,0], dtype='int64')        # Right swath bound in ECI frame
        self.bound_1 = np.array([0,0,0], dtype='int64')            # Left swath bound in lla
        self.bound_2 = np.array([0,0,0], dtype='int64')            # Right swath bound in lla
        self.pos = posvel[0:3]
        self.vel = posvel[3:]
        self.attitude = attitude
        self.time = time

    # ******************** MAKE CONNECT TO SATELLITE FUNCTION HERE ***************************

    def obtain_pointing(self):
        """
            1. Determines the direction unit vector of satellite
            2. Uses direction vector to find nadir angle
            3. Finds length of direction vector from satellite to earths surface
            4. Adds direction vector to position vector to obtain pointing vector
        """

        yaw = self.attitude[0]
        pitch = self.attitude[1]
        roll = self.attitude[2]

        # Obtain direction unit vector
        self.direction = np.array([np.cos(yaw)*np.cos(pitch), np.sin(yaw)*np.cos(pitch), np.sin(pitch)])

        # Obtain Nadir angle (angle between pointing and pos vector)
        nadir = np.arccos(np.dot(self.direction, self.pos)/(np.linalg.norm(self.direction)*np.linalg.norm(self.pos)))
        #print(np.rad2deg(nadir))

        # Finding distance from sat to earth in direction of direction vector
        R_E = 6378137  # 1 is used for testing purposes
        rho = np.arcsin(R_E/np.linalg.norm(self.pos))
        epsi = np.arccos(np.sin(nadir)/np.sin(rho))         # grazing angle
        lam = np.pi/2 - epsi - nadir                        # earth central angle
        D = R_E*(np.sin(lam)/np.sin(nadir))
        self.direction = D * self.direction

        # Obtain ECI vector of point on Earth's surface the satellite is observing
        self.observation = self.pos - self.direction

        # Rotate centre observed point about velocity axis to find swath bounds
        rotation_angle = np.deg2rad(0.175171)           # this is the earth central angle for a 32km swath, calculated by an excel sheet
        axis = self.vel
        self.bound_1_ECI = np.dot(rotation_matrix(axis, rotation_angle), self.observation)
        self.bound_2_ECI = np.dot(rotation_matrix(axis, -rotation_angle), self.observation)

        # convert each bound to an ECEF then to long/lat
        self.bound_1 = ECI2ECEF(self.bound_1_ECI, self.time)
        self.bound_2 = ECI2ECEF(self.bound_2_ECI, self.time)

        self.bound_1 = ECEF2GLLH(self.bound_1)
        self.bound_2 = ECEF2GLLH(self.bound_2)



def rotation_matrix(axis, theta):
    """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d

    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                    [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                    [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


if __name__ == "__main__":

    # TESTING SCRIPT
    posvel = np.array([6878137, 6878137, 6878137, -6878137, -6878137, 2*6878137], dtype='int64')

    yaw = np.arctan2(posvel[1],posvel[0])
    pitch = np.arctan2(posvel[2], np.sqrt(posvel[0]**2 +  posvel[1]**2))
    roll = np.pi/2

    attitude = np.array([yaw + np.deg2rad(5),pitch - np.deg2rad(5),roll])

    camera = payload(attitude, posvel, 0)
    camera.obtain_pointing()

    #print(np.linalg.norm(camera.direction))
    print("position vector", camera.pos, "length =", np.linalg.norm(camera.pos))
    print("direction vector", camera.direction, "length =", np.linalg.norm(camera.direction))
    print("observation vector", camera.observation, "length =", np.linalg.norm(camera.observation))
    print("Bound 1 Lat, Long =", camera.bound_1[0:2])
    print("Bound 2 Lat, Long =", camera.bound_2[0:2])

    fig = plt.figure(1)
    ax = plt.axes(projection='3d')
    ax.quiver
    ax.plot3D([0, posvel[0]], [0, posvel[1]], [0, posvel[2]], 'red', label = "ECI position")
    ax.plot3D([0, posvel[3]], [0, posvel[4]], [0, posvel[5]], 'orange', label = "ECI Velocity")
    ax.plot3D([0, camera.direction[0]], [0, camera.direction[1]], [0, camera.direction[2]], 'blue', label = "Pointing direction")
    ax.plot3D([0, camera.observation[0]], [0, camera.observation[1]], [0, camera.observation[2]], 'green', label = "Observed centre ground point")
    ax.plot3D([0, camera.bound_1_ECI[0]], [0, camera.bound_1_ECI[1]], [0, camera.bound_1_ECI[2]], 'olive')
    ax.plot3D([0, camera.bound_2_ECI[0]], [0, camera.bound_2_ECI[1]], [0, camera.bound_2_ECI[2]], 'olive')
    ax.legend()
    plt.show()