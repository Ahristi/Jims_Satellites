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
from sunSensor import sunSensor

class ADCS:
    def __init__(self, magnetometer = None, startracker = None, sunsensor = None):
        self.magnetometer = magnetometer
        self.starTracker = startracker
        self.satellite = None
        self.sunSensor = sunsensor

    def connectToSatellite(self, satellite):
        self.satellite = satellite

    def determineAttitude(self):
        """
            determineAttitude:

            Determines the satellite attitude based on the sensor readings by applying weighted non linear least squares.        
        """

        #Initial guess
        guess = np.array([0.0,0.0,0.0])
        delta_x = np.array([10,10,10])
    
        max_iter = 400
        tol = 1e-15
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


    def determineAttitudeM(self):
        """
            determineAttitude:

            Determines the satellite attitude based on the sensor readings by applying weighted non linear least squares.        
        """

        #Initial guess
        guess = np.array([0.0,0.0,0.0])
        delta_x = np.array([10,10,10])
    
        max_iter = 400
        tol = 1e-15
        iter_count = 0

        #For now just assume only sensors are magnetometer and startracker
        allSensors = [self.starTracker, self.magnetometer, self.sunSensor]

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

        sunReading = self.sunSensor.getReading(actualAttitude)
        sunReading = sunReading/np.linalg.norm(sunReading)


        #Get unit vectors of the actual readings
        starActual = self.starTracker.getActualReading(starName, self.satellite.X)
        starActual = starActual/np.linalg.norm(starActual)

        magnetActual = self.magnetometer.getActualReading(self.satellite.time, self.satellite.X)
        magnetActual = magnetActual/np.linalg.norm(magnetActual)

        sunActual = self.sunSensor.getActualReading()
        sunActual = sunActual/np.linalg.norm(sunActual)


        while np.linalg.norm(delta_x) > tol:
    


            # Build H matrix
            phi  = guess[0]
            theta = guess[1]
            psi   = guess[2]

            H_starTracker   =   jacobianDCM(psi,theta,phi, starActual)
            H_magnet        =   jacobianDCM(psi,theta,phi, magnetActual)
            H_sun        =   jacobianDCM(psi,theta,phi, sunActual)
      
            H[0:3, :] = H_starTracker
            H[3:6, :] = H_magnet
            H[6:9, :] = H_sun
          
            C = np.array([[np.cos(theta)*np.cos(psi), np.cos(theta)*np.sin(psi), -np.sin(theta)],
                          [np.sin(phi)*np.sin(theta)*np.cos(psi)  - np.cos(phi)*np.sin(psi), np.sin(phi)*np.sin(theta)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.sin(phi)*np.cos(theta)],
                          [np.cos(phi)*np.sin(theta)*np.cos(psi) + np.sin(phi)*np.sin(psi), np.cos(phi)*np.sin(theta)*np.sin(psi)- np.sin(phi)*np.cos(psi), np.cos(phi)*np.cos(theta)]])
            
            #Calculate residuals   
          
            delta_y[0:3, 0] =   (startrackReading - C @ starActual)
            delta_y[3:6, 0] =   (magnetReading    - C @ magnetActual)
            delta_y[6:9, 0] =   (sunReading    - C @ sunActual)
   
         
            #Hard code the weighting matrix for now
            W = np.array([[100/self.starTracker.accuracy**2,0,0,0,0,0,0,0,0],
                        [0,100/self.starTracker.accuracy**2,0,0,0,0,0,0,0],
                        [0,0,100/self.starTracker.accuracy**2,0,0,0,0,0,0],
                        [0,0,0,1/self.magnetometer.accuracy**2,0,0,0,0,0],
                        [0,0,0,0,1/self.magnetometer.accuracy**2,0,0,0,0],
                        [0,0,0,0,0,1/self.magnetometer.accuracy**2,0,0,0],
                        [0,0,0,0,0,0,1/self.sunSensor.accuracy**2,0,0],
                        [0,0,0,0,0,0,0,1/self.sunSensor.accuracy**2,0],
                        [0,0,0,0,0,0,0,0,1/self.sunSensor.accuracy**2]])
    
           #W = np.eye(9)
        

            # Use non-linear least squares to estimate error in x
            delta_x = np.linalg.inv(H.T @ W @  H) @ H.T @ W @ delta_y

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
    print(np.array([roll,pitch,yaw]))

    #Get the angle between the two attitudes for pointing error

    unit = np.array([1,1,1])
    C = directionalCosine(roll,pitch,yaw)
    actualUnit = C @ unit


    estimate = sat.ADCS.determineAttitudeM()
    estimateRoll = estimate[0]
    estimatePitch = estimate[1]
    estimateYaw = estimate[2]

    C = directionalCosine(estimateRoll,estimatePitch,estimateYaw)
    estimateUnit = C @ unit

    




    

        
    

