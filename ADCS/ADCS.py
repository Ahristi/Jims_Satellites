"""
    ADCS (Attitude Determination and Control System)
    
    Object oriented implementation of cubeADCS platform in software simulation.
     
    Note: Control is not implemented, the name is just used for convenience. 
"""

from satellite import Satellite
from orbitalTransforms import *
from magnetometer import Magnetometer
from starTracker import starTracker
from datetime import datetime

class ADCS:
    def __init__(self, magnetometer = None, startracker = None):
        self.magnetometer = magnetometer
        self.starTracker = startracker
        self.satellite = None

    def connectToSatellite(self, satellite):
        self.satellite = satellite

    def determineAttitude(self):
        """
            determineAttitude:

            Determines the satellite attitude based on the sensor readings
            by applying weighted non linear least squares.        
        """

        #Initial guess
        guess = np.array([-2.0,2.0,2.0])
        delta_x = np.array([1,1,1])
    
        max_iter = 3000
        tol = 1e-16
        iter_count = 0

        #For now just assume only sensors are magnetometer and startracker
        allSensors = [self.starTracker, self.magnetometer]

        #3 times length for 3 outputs and 3 rows for 3 inputs
        H = np.zeros((3 * len(allSensors), 3))

        delta_y = np.zeros((3 *  len(allSensors), 1))
        starName = "kentauras"
        actualAttitude = self.satellite.attitude

        #Get unit vector for the sensor readings
        startrackReading = self.starTracker.getReading(starName, self.satellite.X, actualAttitude)
        startrackReading = startrackReading/np.linalg.norm(startrackReading)

        magnetReading = self.magnetometer.getReading(self.satellite.time, self.satellite.X, actualAttitude)
        magnetReading = magnetReading/np.linalg.norm(magnetReading)

        #Get unit vectors of the actual readings
        starActual = self.starTracker.getActualReading(starName, self.satellite.X)
        starActual = starActual/np.linalg.norm(starActual)

        magnetActual = self.magnetometer.getActualReading(self.satellite.time, self.satellite.X)
        magnetActual = magnetActual/np.linalg.norm(magnetActual)



        while np.linalg.norm(delta_x) > tol:
            
            # Build H matrix
            psi   = guess[0]
            theta = guess[1]
            phi  = guess[2]
            

            H_starTracker   =   jacobianDCM(psi,theta,phi, startrackReading)
            H_magnet        =   jacobianDCM(psi,theta,phi, magnetReading)
      
            H[0:3, :] = H_starTracker
            H[3:6, :] = H_magnet
          
            C = np.array([[np.cos(theta)*np.cos(psi), np.cos(theta)*np.sin(psi), -np.sin(theta)],
                          [np.sin(phi)*np.sin(theta)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.sin(phi)*np.sin(theta)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.sin(phi)*np.cos(theta)],
                          [np.cos(phi)*np.sin(theta)*np.cos(psi) + np.sin(phi)*np.sin(psi), np.cos(phi)*np.sin(theta)*np.sin(psi)- np.sin(phi)*np.cos(psi), np.cos(phi)*np.cos(theta)]])
            
            #Calculate residuals   

            #print(self.magnetometer.getActualReading(self.satellite.time, self.satellite.X))
            #print(C.T @ magnetReading)
  
            delta_y[0:3, 0] = (starActual - C.T @ startrackReading)
            delta_y[3:6, 0] = (magnetActual - C.T @ magnetReading)
            
            # Use non-linear least squares to estimate error in x
            delta_x = np.linalg.inv(H.T @ H) @ H.T @ delta_y

            guess += delta_x.flatten()
            
            iter_count += 1

            if iter_count >= max_iter:
                print('Failed to Converge !!')
                break
  
        return guess


    def determineAttitudeINV(self):
        """
            determineAttitude:

            Determines the satellite attitude based on the sensor readings
            by applying weighted non linear least squares.        
        """

        #Initial guess
        guess = np.array([3.0,3.0,-1.0])
        delta_x = np.array([1,1,1])
    
        max_iter = 300
        tol = 1e-16
        iter_count = 0

        #For now just assume only sensors are magnetometer and startracker
        allSensors = [self.starTracker, self.magnetometer]

        #3 times length for 3 outputs and 3 rows for 3 inputs
        H = np.zeros((3 * len(allSensors), 3))

        delta_y = np.zeros((3 *  len(allSensors), 1))
        starName = "kentauras"
        actualAttitude = self.satellite.attitude

        #Get unit vector for the sensor readings
        startrackReading = self.starTracker.getReading(starName, self.satellite.X, actualAttitude)
        startrackReading = startrackReading/np.linalg.norm(startrackReading)



        magnetReading = self.magnetometer.getReading(self.satellite.time, self.satellite.X, actualAttitude)
        magnetReading = magnetReading/np.linalg.norm(magnetReading)

        #Get unit vectors of the actual readings
        starActual = self.starTracker.getActualReading(starName, self.satellite.X)
        starActual = starActual/np.linalg.norm(starActual)

        magnetActual = self.magnetometer.getActualReading(self.satellite.time, self.satellite.X)
        magnetActual = magnetActual/np.linalg.norm(magnetActual)



        while np.linalg.norm(delta_x) > tol:
             
            # Build H matrix
            phi  = guess[0]
            theta = guess[1]
            psi   = guess[2]

            H_starTracker   =   jacobianDCMINV(psi,theta,phi, starActual)
            H_magnet        =   jacobianDCMINV(psi,theta,phi, magnetActual)
      
            H[0:3, :] = H_starTracker
            H[3:6, :] = H_magnet
          
            C = np.array([[np.cos(theta)*np.cos(psi), np.cos(theta)*np.sin(psi), -np.sin(theta)],
                          [np.sin(phi)*np.sin(theta)*np.cos(psi)  - np.cos(phi)*np.sin(psi), np.sin(phi)*np.sin(theta)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.sin(phi)*np.cos(theta)],
                          [np.cos(phi)*np.sin(theta)*np.cos(psi) + np.sin(phi)*np.sin(psi), np.cos(phi)*np.sin(theta)*np.sin(psi)- np.sin(phi)*np.cos(psi), np.cos(phi)*np.cos(theta)]])
            
            #Calculate residuals   

            #print(self.magnetometer.getActualReading(self.satellite.time, self.satellite.X))
            #print(C.T @ magnetReading)
          
            delta_y[0:3, 0] = (startrackReading - C @ starActual)
            delta_y[3:6, 0] = (magnetReading    - C @ magnetActual)


            # Use non-linear least squares to estimate error in x
            delta_x = np.linalg.inv(H.T @ H) @ H.T @ delta_y

            guess += delta_x.flatten()
            
            iter_count += 1

            if iter_count >= max_iter:
                print('Failed to Converge !!')
                break
  
        return guess


def jacobianDCM(yaw, pitch, roll, obs):
    """ 
        jacobianDCM

        Computes the jacobian of the directional cosine matrix

        Inputs;

        yaw:    yaw (psi) euler angle of the satellite in ECI frame
        pitch:  pitch (theta) euler angle of the satellite in ECI frame 
        roll:   roll (phi) euler angle of the satellite in ECI frame
        obs:    observation made of a known vector in the body frame

    """
    #Convenience variables
    xb      =  obs[0]
    yb      =  obs[1]
    zb      =  obs[2]
    psi     =  yaw
    theta   =  pitch
    phi     =  roll  

    
    dxdpsi      =   -xb * np.cos(theta) * np.sin(psi) + yb * (-np.sin(phi) * np.sin(theta) * np. sin(psi) - np.cos(phi) * np.cos(psi)) + zb * (-np.cos(phi) * np.sin(theta) * np.sin(psi) + np.sin(theta) * np.cos(psi))
    dxdtheta    =   -xb*np.sin(theta) * np.cos(psi)   + yb * np.sin(phi) * np.cos(theta) * np.cos(psi)+ zb * np.cos(phi) * np.cos(theta) * np.cos(psi)
    dxdphi      =    yb*(np.cos(phi) * np.sin(theta) * np.cos(psi) + np.sin(phi) * np.sin(psi)) + zb* (-np.sin(phi) * np.sin(theta)*np.cos(psi) + np.cos(phi) * np.sin(psi))
    
    dydpsi      =   xb*np.cos(theta) * np.cos(psi) + yb*(np.sin(phi) * np.sin(theta) * np.cos(psi) - np.cos(phi) * np.sin(psi)) + zb * (np.cos(phi) * np.sin(theta) * np.cos(psi) + np.sin(phi) * np.sin(psi))
    dydtheta    =   -xb * np.sin(theta) * np.sin(psi) + yb * np.sin(phi) * np.cos(theta) * np.sin(psi) + zb* np.cos(phi) * np.cos(theta) * np.sin(psi)
    dydphi      =   yb*(np.cos(phi) * np.sin(theta) * np.sin(psi) - np.sin(phi) * np.cos(psi)) + zb* (-np.sin(phi) * np.sin(theta) * np.sin(psi) - np.cos(phi) * np.cos(psi))

    dzdpsi      =   0
    dzdtheta    =   -xb*np.cos(theta) - yb* np.sin(phi) * np.sin(theta) - zb * np.cos(phi) * np.sin(theta)
    dzdphi      =   yb * np.cos(phi) * np.cos(theta) - zb * np.sin(phi) * np.cos(theta)


    H = np.array([[dxdpsi, dxdtheta, dxdphi],
                 [dydpsi, dydtheta, dydphi],
                 [dzdpsi, dzdtheta, dzdphi]])
    
    return H


def jacobianDCMINV(yaw, pitch, roll, obs):
    """ 
        jacobianDCM

        Computes the jacobian of the directional cosine matrix

        Inputs;

        yaw:    yaw (psi) euler angle of the satellite in ECI frame
        pitch:  pitch (theta) euler angle of the satellite in ECI frame 
        roll:   roll (phi) euler angle of the satellite in ECI frame
        obs:    observation made of a known vector in the body frame

    """
    #Convenience variables
    x      =  obs[0]
    y      =  obs[1]
    z      =  obs[2]
    psi     =  yaw
    theta   =  pitch
    phi     =  roll  

    H11 = 0
    H12 = -np.cos(psi) * np.sin(theta) * x - np.sin(psi) * np.sin(theta) * y - np.cos(theta) * z
    H13 = -np.sin(psi) * np.cos(theta) * x + np.cos(psi) * np.cos(theta) * y

    H21 = (np.cos(psi) * np.sin(theta) *np.cos(phi) + np.sin(psi) * np.sin(phi)) * x + (np.sin(psi) * np.sin(theta) * np.cos(phi) - np.cos(psi) * np.sin(phi)) * y + np.cos(theta) * np.cos(phi) * z
    H22 = np.cos(psi) * np.cos(theta) * np.sin(phi) * x + np.sin(psi) * np.cos(theta) * np.sin(phi) * y - np.sin(theta) * np.sin(phi) * z
    H23 = (-np.sin(psi) * np.sin(theta) * np.sin(phi)  - np.cos(psi) * np.cos(phi)) * x + (np.cos(psi) * np.sin(theta) * np.sin(phi) - np.sin(psi) * np.cos(phi)) * y

    H31 = (-np.cos(psi) * np.sin(theta) * np.sin(phi) + np.sin(psi) * np.cos(phi)) * x + (-np.sin(psi) * np.sin(theta) * np.sin(phi) - np.cos(psi) * np.cos(phi)) * y - np.cos(theta) * np.sin(phi) * z
    H32 = np.cos(psi) * np.cos(theta) * np.cos(phi) * x + np.sin(psi) * np.cos(theta) * np.cos(phi) * y - np.sin(theta) * np.cos(phi) * z
    H33 = (-np.sin(psi) * np.sin(theta) * np.cos(phi) + np.cos(psi) * np.sin(phi)) * x + (np.cos(psi) * np.sin(theta) * np.cos(phi) + np.sin(psi) * np.sin(phi))*y


    H = np.array([[H11, H12, H13],
                 [H21, H22, H23],
                 [H31, H32, H33]])
    
    return H

    
    

if __name__ == "__main__": 
    cross = starTracker("star_config.csv", 0.0)
    J2000 = datetime(2000, 1, 1, 12)
    mg = Magnetometer(J2000,0)
    adcs = ADCS(mg,cross)

    sat = Satellite("ISS.txt", J2000, adcs)
    pitch = np.arctan2(sat.X[1], sat.X[2])
    yaw   = np.arctan2(sat.X[0], sat.X[1])
    roll = np.arctan2(np.sin(pitch)*np.cos(yaw), np.cos(pitch)*np.cos(yaw))

    sat.setAttitude(np.array([roll,pitch,yaw]))
    sat.ADCS.connectToSatellite(sat)
    print(np.array([roll,pitch,yaw]))
    print(sat.ADCS.determineAttitudeINV())

    

        
    

