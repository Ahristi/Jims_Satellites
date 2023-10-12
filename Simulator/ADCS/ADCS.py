"""
    ADCS (Attitude Determination and Control System)
    
    Object oriented implementation of cubeADCS platform in software simulation.
     
    Note: Control is not implemented, the name is just used for convenience. 
"""

from satellite import Satellite
from orbitalTransforms import *
from Simulator.magnetometer import Magnetometer
from Simulator.starTracker import starTracker
from datetime import datetime
from Simulator.sunSensor import sunSensor

NUM_STAR_TRACKER = 2 


class ADCS:
    def __init__(self, magnetometer = None, startracker = None, sunsensor = None):
        
        self.magnetometer = magnetometer        
        self.starTracker1 = startracker         
        self.starTracker2 = startracker
        self.sunSensor = sunsensor

        self.satellite = None

    def connectToSatellite(self, satellite):
        self.satellite = satellite


    def determineFineAttitude(self):
        """
            determineAttitude:

            Determines the satellite attitude based on the sensor readings by applying weighted non linear least squares.        
        """

        #Initial guess
        guess = np.array([0,0,0])
        delta_x = np.array([10,10,10])
    
        max_iter = 1000
        tol = 1e-15
        iter_count = 0

        #3 times length for 3 outputs and 3 rows for 3 inputs
        H = np.zeros((3 * NUM_STAR_TRACKER, 3))

        delta_y = np.ones((3 *  NUM_STAR_TRACKER, 1))*1000


        starName = "kentauras" #This is hard coded but we can add multiple stars later
        


        #Get readings and convert to unit vector
        pos = self.satellite.states[-1][0:3]          #Current pos is the most recent one
        actualAttitude = self.satellite.attitudes[-1] #Current attitude is the most recent one

        startrack1Reading = self.starTracker.getReading(starName, pos, actualAttitude)
        startrack1Reading = startrack1Reading/np.linalg.norm(startrack1Reading)

        startrack2Reading = self.starTracker.getReading(starName, pos, actualAttitude)
        startrack2Reading = startrack2Reading/np.linalg.norm(startrack2Reading)

        #Get unit vector of the actual reading
        starActual = self.starTracker.getActualReading(starName, self.satellite.X)
        starActual = starActual/np.linalg.norm(starActual)


        while np.linalg.norm(delta_x) > tol:


            # Build H matrix
            phi  = guess[0]
            theta = guess[1]
            psi   = guess[2]

            H_starTracker   =   jacobianDCM(psi,theta,phi, starActual)

            H[0:3, :] = H_starTracker
            H[3:6, :] = H_starTracker
          
            C = directionalCosine(phi,theta,psi)
            
            
            #Calculate residuals   
            delta_y[0:3, 0] =   (startrack1Reading - C @ starActual)  
            delta_y[3:6, 0] =   (startrack2Reading    - C @ starActual) 

            #Both star trackers are the same so weight them the same
            W = np.eye(9)
            #TODO: Maybe make one star tracker better than the other to pretend one was manufactured a bit better or is in a better part of the satellite

            # Use non-linear least squares to estimate error in x
            delta_x = np.linalg.inv(H.T @ W @  H) @ H.T @ W @ delta_y
            guess += delta_x.flatten()
      
            #Deal with angle wrap around
            guess[0] = guess[0] % (2*np.pi)
            guess[2] = guess[2] % (2*np.pi)
            iter_count += 1

            if iter_count >= max_iter:
                break
  
        return guess



    def determineAttitude(self):
        """
            determineAttitude:

            Determines the satellite attitude based on the sensor readings by applying weighted non linear least squares.        
        """

        #Initial guess
        guess = np.array([1.5,-0.8,3])
        delta_x = np.array([10,10,10])
    
        max_iter = 1000
        tol = 1e-15
        iter_count = 0

        #For now just assume only sensors are magnetometer and startracker
        allSensors = [self.starTracker, self.magnetometer, self.sunSensor]

        #3 times length for 3 outputs and 3 rows for 3 inputs
        H = np.zeros((3 * len(allSensors), 3))

        delta_y = np.ones((3 *  len(allSensors), 1))*1000
        starName = "kentauras"
        actualAttitude = self.satellite.attitude

        #Get unit vector for the sensor readings
        startrackReading = self.starTracker.getReading(starName, self.satellite.X, actualAttitude)
        startrackReading = startrackReading/np.linalg.norm(startrackReading)

        magnetReading = self.magnetometer.getReading(self.satellite.time, self.satellite.X, actualAttitude)
        magnetReading = magnetReading/np.linalg.norm(magnetReading)

        sunReading = self.sunSensor.getReading(actualAttitude)
        sunReading =sunReading/np.linalg.norm(sunReading)


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
            H_sun           =   jacobianDCM(psi,theta,phi, sunActual)
  
            H[0:3, :] = H_starTracker
            H[3:6, :] = H_starTracker
          
            C = directionalCosine(phi,theta,psi)
            
            
            #Calculate residuals   
            delta_y[0:3, 0] =   (startrackReading - C @ starActual)  
            delta_y[3:6, 0] =   (magnetReading    - C @ magnetActual) 
            delta_y[6:9, 0] =   (sunReading    - C @ sunActual) 
         
            """
                #Hard code the weighting matrix for now
                W = np.array([[self.starTracker.accuracy**2/self.starTracker.accuracy,0,0,0,0,0,0,0,0],
                            [0,self.starTracker.accuracy**2/self.starTracker.accuracy,0,0,0,0,0,0,0],
                            [0,0,self.starTracker.accuracy**2/self.starTracker.accuracy,0,0,0,0,0,0],
                            [0,0,0,self.starTracker.accuracy**2/self.magnetometer.accuracy,0,0,0,0,0],
                            [0,0,0,0,self.starTracker.accuracy**2/self.magnetometer.accuracy,0,0,0,0],
                            [0,0,0,0,0,self.starTracker.accuracy**2/self.magnetometer.accuracy,0,0,0],
                            [0,0,0,0,0,0,self.starTracker.accuracy**2/self.sunSensor.accuracy,0,0],
                            [0,0,0,0,0,0,0,self.starTracker.accuracy**2/self.sunSensor.accuracy,0],
                            [0,0,0,0,0,0,0,0,self.starTracker.accuracy**2/self.sunSensor.accuracy]])
            """
            #Hard code the weighting matrix for now
            W = np.array([[1000,0,0,0,0,0,0,0,0],
                        [0,1000,0,0,0,0,0,0,0],
                        [0,0,1000,0,0,0,0,0,0],
                        [0,0,0,0.01,0,0,0,0,0],
                        [0,0,0,0,0.01,0,0,0,0],
                        [0,0,0,0,0,0.01,0,0,0],
                        [0,0,0,0,0,0,0.01,0,0],
                        [0,0,0,0,0,0,0,0.01,0],
                        [0,0,0,0,0,0,0,0,0.01]])
            W = np.eye(9)
      
            # Use non-linear least squares to estimate error in x
            delta_x = np.linalg.inv(H.T @ W @  H) @ H.T @ W @ delta_y
            guess += delta_x.flatten()
      

            guess[0] = guess[0] % (2*np.pi)
      
            guess[2] = guess[2] % (2*np.pi)
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
    #I can't work out the weights correctly so just model the star tracker using the magnetometer and sun sensor


    #Make the star tracker
    cross = starTracker("star_config.csv", 0.0027)
    
    #Make the magnetometer
    J2000 = datetime(2000, 1, 1, 12)
    mg = Magnetometer(J2000,0.0027)

    #Make the sun sensor
    sunsens = sunSensor(0.0027)

    adcs = ADCS(mg,cross,sunsens)

    sat = Satellite("ISS.txt", J2000, adcs)

    sat.setAttitude()
    sat.ADCS.connectToSatellite(sat)
    print(sat.attitude)

    #Get the angle between the two attitudes for pointing error

    unit = np.array([1,1,1])
    C = directionalCosine(sat.attitude[0],sat.attitude[1],sat.attitude[2])
    actualUnit = C @ unit
    sigma = 0
 
    for i in range(1000):
        estimate = sat.ADCS.determineFineAttitude()
        estimateRoll = estimate[0] 
        estimatePitch = estimate[1] 
        estimateYaw = estimate[2]   
        C = directionalCosine(estimateRoll,estimatePitch,estimateYaw)
        estimateUnit = C @ unit
        sigma+=angleBetweenVectors(estimateUnit,actualUnit)
        print(angleBetweenVectors(estimateUnit,actualUnit))
    print(sigma/1000)
    

