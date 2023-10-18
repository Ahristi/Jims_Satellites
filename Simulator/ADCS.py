"""
    ADCS (Attitude Determination and Control System)
    
    Object oriented implementation of cubeADCS platform in software simulation.
     
    Note: Control is not implemented, the name is just used for convenience. 
"""

from orbitalTransforms import *
from magnetometer import Magnetometer
from starTracker import starTracker
from datetime import datetime
from sunSensor import sunSensor

NUM_STAR_TRACKER = 2 
MAGNETOMETER_ACCURACY = 0.1     #Degrees
STARTRACKER_ACCURACY  = 0.0027  #Degrees
SUNSENSOR_ACCURACY    = 0.1     #Degrees

#Satellite states 
SAFE    =   0
IMAGING =   1
READOUT =   2


class ADCS:
    def __init__(self, starTrackerFile, sat):
        """
            ADCS object which simulates attitude determination in the satellite
            Contains different sensor objects to determine the attitude using WNLLS

            Inputs:

            starTrackerFile     -   Name of the csv file containing the star tracker config
            sat                 -   Satellite object that the ADCS belongs to
        
        """
        
        self.magnetometer = Magnetometer(datetime(2000, 1, 1, 12), MAGNETOMETER_ACCURACY) #Just hard coding epoch
        self.starTracker1 = starTracker(starTrackerFile, STARTRACKER_ACCURACY)      
        self.starTracker2 = starTracker(starTrackerFile, STARTRACKER_ACCURACY)     
        self.sunSensor    = sunSensor(SUNSENSOR_ACCURACY)

        self.satellite = sat
        self.attitude = np.array([0,0,0])

        #Datalogging
        self.estimatedAttitudes = []
        self.determineAttitude() #Get initial attitude

    def connectToSatellite(self, satellite):
        """
            Legacy code that connects a satellite object to its ADCS. 
            This is done in the constructor now but keep just in case for testing.
        """
        self.satellite = satellite

    def determineAttitude(self):
        if (self.satellite.state == IMAGING):
            self.determineCoarseAttitude()
        elif (self.satellite.state == SAFE):
            self.determineCoarseAttitude()
        else:
            self.determineCoarseAttitude()

        self.estimatedAttitudes.append(self.attitude)

    def determineFineAttitude(self):
        """
            determineAttitude:

            Determines the satellite attitude based on the two star tracker readings by applying weighted non linear least squares.        
        """

        #Initial guess
        guess = np.array([0.0,0.0,0.0])
        delta_x = np.array([10,10,10])
    
        max_iter = 1000
        tol = 1e-8
        iter_count = 0

        #3 times length for 3 outputs and 3 rows for 3 inputs
        H = np.zeros((3 * NUM_STAR_TRACKER, 3))
        delta_y = np.ones((3 *  NUM_STAR_TRACKER, 1))*1000

        X = self.satellite.states[-1][0:3]


        starName = "kentauras" #This is hard coded but we can add multiple stars later
        
        #Get readings and convert to unit vector
        pos = self.satellite.states[-1][0:3]          #Current pos is the most recent one
        actualAttitude = self.satellite.attitudes[-1] #Current attitude is the most recent one

        startrack1Reading = self.starTracker1.getReading(starName, pos, actualAttitude)
        startrack1Reading = startrack1Reading/np.linalg.norm(startrack1Reading)

        startrack2Reading = self.starTracker2.getReading(starName, pos, actualAttitude)
        startrack2Reading = startrack2Reading/np.linalg.norm(startrack2Reading)

        #Get unit vector of the actual reading
        starActual = self.starTracker1.getActualReading(starName, X)
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
            W = np.eye(6)
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


    def getDesiredAttitude(self):
        """
            Gets the desired attitude of the satellite. 
            For the time being attitude change is assumed to be instant.

            #TODO: Consider satellite state
        """
        x = self.satellite.states[-1][0]
        y = self.satellite.states[-1][1]
        z = self.satellite.states[-1][2]

        yaw = np.arctan2(y,x)
        pitch = np.arctan2(z, np.sqrt(x**2 +  y**2))
        roll = 0
        desiredAttitude = np.array([roll,pitch,yaw])
        return desiredAttitude
    
    def getQuaternions(self, attitude):
        """
            Convert the satellite attitude to quaternions
        
        """
        roll,pitch,yaw = attitude

        cr = np.cos(roll * 0.5)
        sr = np.sin(roll * 0.5)
        cp = np.cos(pitch * 0.5)
        sp = np.sin(pitch * 0.5)
        cy = np.cos(yaw * 0.5)
        sy = np.sin(yaw * 0.5)
 
        qw = cr * cp * cy + sr * sp * sy
        qx = sr * cp * cy - cr * sp * sy
        qy = cr * sp * cy + sr * cp * sy
        qz = cr * cp * sy - sr * sp * cy


        q = np.array([qw,qx,qy,qz])
        return q
    


    def determineCoarseAttitude(self):
        """
            determineCoarseAttitude:

            Determines the satellite attitude based on coarse sensor readings by applying weighted non linear least squares.        
        """

        #Initial guess
        guess = np.array([1.0,1.0,1.0,1.0])
        delta_x = np.array([10,10,10,10])
    
        max_iter = 1000
        tol = 1e-8
        iter_count = 0

        #For now just assume only sensors are magnetometer and startracker
        allSensors = [self.magnetometer, self.sunSensor]

        #3 times length for 3 outputs and 3 rows for 3 inputs
        H = np.zeros((4 * len(allSensors), 4))

        delta_y = np.ones((4 *  len(allSensors), 1))*1000

        actualAttitude = self.satellite.attitude
        X = self.satellite.states[-1][0:3]
        
        #Get readings from our sensors
        magnetReadingQ = np.zeros(4)
        magnetReading = self.magnetometer.getReading(self.satellite.time, X, actualAttitude)
        magnetReading = magnetReading/np.linalg.norm(magnetReading)
        magnetReadingQ[1:4] = magnetReading

        sunReadingQ = np.zeros(4)
        sunReading = self.sunSensor.getReading(actualAttitude)
        sunReading =sunReading/np.linalg.norm(sunReading)
        sunReadingQ[1:4] = sunReading

        #Get the actual values which are known references
        magnetActualQ = np.zeros(4)
        magnetActual = self.magnetometer.getActualReading(self.satellite.time, X)
        magnetActual = magnetActual/np.linalg.norm(magnetActual)
        magnetActualQ[1:4] = magnetActual

        sunActualQ = np.zeros(4)
        sunActual = self.sunSensor.getActualReading()
        sunActual = sunActual/np.linalg.norm(sunActual)
        sunActualQ[1:4] = sunActual


        while np.linalg.norm(delta_x) > tol:
            H_magnet        =   quarternionJacobian(guess, magnetReadingQ)
    
            H_sun           =   quarternionJacobian(guess, sunReadingQ)
            H[0:4, :] = H_magnet
            H[4:8, :] = H_sun
          
            q = guess
            q /= np.linalg.norm(q)               
            q_ = np.array([q[0],-q[1],-q[2],-q[3]])

            #Calculate residuals    
            delta_y[0:4, 0] =   (magnetActualQ    - quaternionMultiply(quaternionMultiply(q,magnetReadingQ),q_)) 
            delta_y[4:8, 0] =   (sunActualQ    - quaternionMultiply(quaternionMultiply(q,sunReadingQ),q_)) 
      
            #Don't worry about weightings for the time being.
            W = np.eye(8)   
      
            # Use non-linear least squares to estimate error in x
            delta_x = np.linalg.inv(H.T @ W @  H) @ H.T @ W @ delta_y
            guess += delta_x.flatten()
      

            iter_count += 1

            if iter_count >= max_iter:
                #print('Failed to Converge !!')
                break
        e = quaternion2Euler(guess)

        self.attitude = e



def restrictDomain(eulerAngles):
    """
        Restricts the domain of the euler angles as follows:
        roll: -pi to pi
        pitch: -pi/2 to pi/2
        yaw:  -pi to pi

        Inputs:
        eulerAngles - the attitude guess determined by NLLS
    
    """
    guess = eulerAngles



    #Pitch
    if (guess[1] > np.pi):
        guess[1] =np.pi - guess[1]
        #check for rare double wrap
        if (guess[1] < -np.pi):
            guess[1] = guess[1]+ np.pi/2

    #Yaw
    if (guess[2] > np.pi):
        guess[2] = -np.pi + (guess[2] -np.pi)


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

