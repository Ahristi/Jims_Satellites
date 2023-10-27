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
from orbitalTransforms import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


#Satellite States
SAFE    =   0
IMAGING =   1
READOUT =   2

class Payload:
    def __init__(self,sat):
        """
            Payload Object
            
            Inputs:
            attitude:   [yaw, pitch, roll] estimated attitude of satellite in ECI frame
            posvel:     [Xp, Yp, Zp, Xv, Yv, Zv] estimated position/velocity of satellite in ECI frame
            time:       time since vernal in seconds
        """
        self.satellite = sat
        self.direction = np.array([0,0,0], dtype='int64')          # vector decribing satellites pointing direction in ECO
        self.observation = np.array([0,0,0], dtype='int64')        # vector describing ECI position satellite is pointing at

        self.bound_1_ECI = np.array([0,0,0], dtype='int64')        # Left swath bound in ECI frame
        self.bound_1_2_ECI = np.array([0,0,0], dtype='int64')      # Second left swath bound in ECI frame
        self.bound_2_ECI = np.array([0,0,0], dtype='int64')        # Right swath bound in ECI frame
        self.bound_2_2_ECI = np.array([0,0,0], dtype='int64')      # Second right swath bound in ECI frame
        self.bound_1 = np.array([0,0,0], dtype='int64')            # Left swath bound in lla
        self.bound_1_2 = np.array([0,0,0], dtype='int64')          # Second left swath bound in lla
        self.bound_2 = np.array([0,0,0], dtype='int64')            # Right swath bound in lla
        self.bound_2_2 = np.array([0,0,0], dtype='int64')          # Second right swath bound in lla

        self.images        = []                                    #Bounding boxes of each image taken                    
        self.mappingErrors = []                                    #Magnitude of the distance between image centere and estimated centre


    def obtainPointing(self):
        """
            1. Determines the direction unit vector of satellite
            2. Uses direction vector to find nadir angle
            3. Finds length of direction vector from satellite to earths surface
            4. Adds direction vector to position vector to obtain pointing vector
        """
 
        #Dodgy Maneouver to not image when we aren't over NSW.
        if (self.satellite.state != IMAGING):
            return


        pos      = self.satellite.GNSS.position
        vel      = self.satellite.states[-1][3:6]
        attitude = self.satellite.ADCS.attitude 

        roll   =  attitude[0]
        pitch  =  attitude[1]
        yaw    =  attitude[2]        
        t      =  self.satellite.times[-1] + self.satellite.tSinceVernal

        # Obtain direction unit vector
        self.direction = np.array([np.cos(yaw)*np.cos(pitch), np.sin(yaw)*np.cos(pitch), np.sin(pitch)])

        # Obtain Nadir angle (angle between pointing and pos vector)
        nadir = np.arccos(np.dot(self.direction, pos)/(np.linalg.norm(self.direction)*np.linalg.norm(pos)))
        #print(np.rad2deg(nadir))

        # Finding distance from sat to earth in direction of direction vector
        R_E = 6378137  # 1 is used for testing purposes
        rho = np.arcsin(R_E/np.linalg.norm(pos))
        epsi = np.arccos(np.sin(nadir)/np.sin(rho))         # grazing angle
        lam = np.pi/2 - epsi - nadir                        # earth central angle
        D = R_E*(np.sin(lam)/np.sin(nadir))
        self.direction = D * self.direction

        # Obtain ECI vector of point on Earth's surface the satellite is observing
        self.observation = pos - self.direction
        self.obtainMappingError(lam)

        # Rotate centre observed point about velocity axis to find swath bounds
        rotation_angle = np.deg2rad(0.175171)               # this is the earth central angle for a 32km swath, calculated by an excel sheet
        axis = vel
        self.bound_1_ECI = np.dot(rotation_matrix(axis, rotation_angle), self.observation)
        self.bound_2_ECI = np.dot(rotation_matrix(axis, -rotation_angle), self.observation)

        rotation_angle = np.deg2rad(8.98e-5)               # this is the earth central angle for a 20m bound thickness
        axis2 = np.cross(vel, self.observation)/2e7
        self.bound_1_2_ECI = np.dot(rotation_matrix(axis2, rotation_angle), self.bound_1_ECI)
        self.bound_2_2_ECI = np.dot(rotation_matrix(axis2, rotation_angle), self.bound_2_ECI)

        # convert each bound to an ECEF then to long/lat
        self.bound_1 = ECI2ECEF(self.bound_1_ECI, t)
        self.bound_2 = ECI2ECEF(self.bound_2_ECI, t)

        self.bound_1 = ECEF2GLLH(self.bound_1)  #Tony this is the worst variable naming convention I have ever seen. I'm leaving it here so the markers laugh at you
        self.bound_2 = ECEF2GLLH(self.bound_2)

        self.bound_1_2 = ECI2ECEF(self.bound_1_2_ECI, t)
        self.bound_2_2 = ECI2ECEF(self.bound_2_2_ECI, t)

        self.bound_1_2 = ECEF2GLLH(self.bound_1_2)
        self.bound_2_2 = ECEF2GLLH(self.bound_2_2)

        new_photo = [self.bound_1, self.bound_2, self.bound_1_2, self.bound_2_2]
        self.images.append(new_photo)

    def obtainMappingError(self, eca):
        """
            Obtains the mapping err or at the current observation in magnitude of meters
            from the true observation.

            Inputs:

            eca - the earth central angle (RADIANS) between the observation and true pointing mad in the obtainPointing function.
        
        """
        R_E = 6378137

        mappingError = eca/(2*np.pi) * 2*np.pi*R_E
        self.mappingErrors.append(mappingError)
    



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

    camera = Payload(attitude, posvel, 0)
    camera.obtainPointing()

    print("position vector", camera.pos, "length =", np.linalg.norm(camera.pos))
    print("direction vector", camera.direction, "length =", np.linalg.norm(camera.direction))
    print("observation vector", camera.observation, "length =", np.linalg.norm(camera.observation))
    print("Bound 1 Lat, Long =", camera.bound_1[0:2])
    print("Bound 1_2 Lat, Long =", camera.bound_1_2[0:2])
    print("Bound 2 Lat, Long =", camera.bound_2[0:2])
    print("Bound 2_2 Lat, Long =", camera.bound_2_2[0:2])

    axis2 = np.cross(camera.vel, camera.observation)/2e7

    fig = plt.figure(1)
    ax = plt.axes(projection='3d')
    ax.quiver
    ax.plot3D([0, posvel[0]], [0, posvel[1]], [0, posvel[2]], 'red', label = "ECI position")
    ax.plot3D([0, posvel[3]], [0, posvel[4]], [0, posvel[5]], 'orange', label = "ECI Velocity")
    ax.plot3D([0, axis2[0]], [0, axis2[1]], [0, axis2[2]], 'pink', label = "axis2")
    ax.plot3D([0, camera.direction[0]], [0, camera.direction[1]], [0, camera.direction[2]], 'blue', label = "Pointing direction")
    ax.plot3D([0, camera.observation[0]], [0, camera.observation[1]], [0, camera.observation[2]], 'green', label = "Observed centre ground point")
    ax.plot3D([0, camera.bound_1_ECI[0]], [0, camera.bound_1_ECI[1]], [0, camera.bound_1_ECI[2]], 'olive')
    ax.plot3D([0, camera.bound_2_ECI[0]], [0, camera.bound_2_ECI[1]], [0, camera.bound_2_ECI[2]], 'olive')
    ax.plot3D([0, camera.bound_1_2_ECI[0]], [0, camera.bound_1_2_ECI[1]], [0, camera.bound_1_2_ECI[2]], 'olive')
    ax.plot3D([0, camera.bound_2_2_ECI[0]], [0, camera.bound_2_2_ECI[1]], [0, camera.bound_2_2_ECI[2]], 'olive')
    ax.legend()
    plt.show()