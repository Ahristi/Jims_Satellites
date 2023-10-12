import numpy as np


#Global Variables
u = 3.986004418*10**14  
R = 6378137
r = 3064192.46769972



"""
     ECI2ECEF

     Inputs:
     X - Position of an object in ECI frame
     t - Time since vernal equinox

     Outputs:
     3x1 vector representing the objects position in ECEF frame

"""
def ECI2ECEF(X, t):
    w     =    7.292115*10**(-5)
    C     =    np.array([[np.cos(w*t), -np.sin(w*t), 0], 
                         [np.sin(w*t), np.cos(w*t), 0],
                         [0,0,1]])
    ECEF = np.matmul(C,X)
    return ECEF


"""
     ECEF2GLLH
     Inputs:
     X_ECEF - Position of the satellite in ECEF frame.

     Outputs:
     3x1 vector representing the position of the satellite in geodetic LLH 
"""
def ECEF2GLLH(X):
    x       =   X[0]
    y       =   X[1]
    z       =   X[2]

    e2      =   0.00669437999
    a       =   6378137
    b       =   a*np.sqrt(1-e2)
    r       =   np.sqrt(x**2 + y**2)

    l       =     e2/2
    m       =   (r/a)**2
    n       =   ((1-e2)*z/b)**2
    i       =   -(2*l**2 + m + n)/2
    k       =   l**2 * (l**2 - m - n)
    q       =   (m+n-4*l**2)**3 / 216 + m*n*l**2
    D       =   np.sqrt((2*q - m*n*l**2)*m*n*l**2) 
    beta    =   i/3 - (q+D)**(1/3) - (q-D)**(1/3)
    t       =   np.sqrt(np.sqrt(beta**2 - k ) - (beta+i)/2) - np.sign(m-n) * np.sqrt((beta-i)/2)
    r0      =   r/(t+l)
    z0      =   (1-e2)*z/(t-l)
    lat     =   np.arctan(z0/((1-e2)*r0))
    long    =   2*np.arctan((np.sqrt(x**2+y**2)-x)/ y)
    h       =   np.sign(t-1+l)*np.sqrt((r-r0)**2 + (z-z0)**2)
    return np.array([np.rad2deg(lat), np.rad2deg(long), h])

"""
     GLLH2ECEF

     Inputs: X in geodetic lat long height

     Outputs: object position in ECEF frame

"""
def GLLH2ECEF(X):
    lat     =   np.deg2rad(X[0])
    long    =   np.deg2rad(X[1])
    h       =   X[2]
    a       =   6378137
    e2      =   0.00669437999 # eccentricity of ellipsoid for WGS-84
    R       =   a/np.sqrt(1-e2 * np.sin(lat)**2)

    x       =   (R+h)*np.cos(long)*np.cos(lat)
    y       =   (R+h)*np.cos(lat)*np.sin(long)
    z       =   (R+h - e2 * R)*np.sin(lat)
    
    return np.array([x,y,z])



def PFF2ECI(X,V,O,w,i):
    """
        Converts satellite state from 2D perifocal frame to ECI

        Inputs:
        X: Position vector of the satellite in perifocal frame
        V: Velocity vector of the satellite in perifocal frame
        O: Right ascension of ascending node
        w: Argument of perigee
        i: Inclination
    """
    C = np.array(
            [   
                [-np.sin(O)*np.cos(i)*np.sin(w) + np.cos(O)*np.cos(w), -np.sin(O)*np.cos(i)*np.cos(w) - np.cos(O)*np.sin(w), np.sin(O)*np.sin(i)],
                [np.cos(O)*np.cos(i)*np.sin(w) + np.sin(O)*np.cos(w), np.cos(O)*np.cos(i)*np.cos(w) - np.sin(O)*np.sin(w), np.cos(O)*np.sin(i)],
                [np.sin(i)*np.sin(w),                          np.sin(i)*np.cos(w),                      np.cos(i)      ]
            ])
    X_ECEF=  np.matmul(C,X)
    V_ECEF=  np.matmul(C,V)
    return X_ECEF, V_ECEF


"""
     ECEF2POLAR

     Inputs:
     X - Position of an object in ECEF frame

     Outputs:
     3x1 vector containing elevation, azimuth and distance

"""

def CART2POLAR(X):
    x     =   X[0]
    y     =   X[1]
    z     =   X[2]
    R     =   np.sqrt(x**2 + y**2 + z**2)

    el    =   np.arcsin(z/R)
    az    =   np.arctan2(y,x)

    return el,az, R


"""
     POLAR2CART

     Inputs:
     X - Position of an object in ECEF frame

     Outputs:
     3x1 vector containing elevation, azimuth andx distance

"""
def POLAR2CART(X):
     el    =    X[0]
     az    =    X[1]
     r     =    X[2]


     x     =    r * np.cos(az)*np.cos(el)
     y     =    r * np.sin(az)*np.cos(el)
     z     =    r * np.sin(el)

     return np.array([x,y,z])


def ECEF2LGDV(X, ref):
    lat = np.deg2rad(ref[0])
    long = np.deg2rad(ref[1])
    R = ref[2]
    C = np.array(
        [
            [-np.sin(long), np.cos(long), 0],
            [-np.cos(long)*np.sin(lat), -np.sin(long)*np.sin(lat), np.cos(lat)],
            [np.cos(lat)*np.cos(long), np.cos(lat)*np.sin(long), np.sin(lat)]
        ]
    )

    ref_ECI = GLLH2ECEF(ref)
    X_ENU = np.matmul(C,X-ref_ECI)
    return X_ENU




def LGDV2ECEF(X_NEU, ref):
    lat = np.deg2rad(ref[0])
    long = np.deg2rad(ref[1])
    C = np.array(
        [
            [-np.sin(long), -np.sin(lat)*np.sin(long), np.cos(lat)*np.cos(long)],
            [np.cos(long), -np.sin(lat)*np.sin(long), np.cos(lat)*np.sin(long)],
            [0, np.cos(lat), np.sin(lat)]
        ]
    )
    C = np.array(
        [
            [-np.sin(long), np.cos(long), 0],
            [-np.cos(long)*np.sin(lat), -np.sin(long)*np.sin(lat), np.cos(lat)],
            [np.cos(lat)*np.cos(long), np.cos(lat)*np.sin(long), np.sin(lat)]
        ]
    )
    ref_ECEF = GLLH2ECEF(ref)
    X_ECEF = np.matmul(np.linalg.inv(C),X_NEU) + np.array(ref_ECEF)
    return X_ECEF

def ECEF2ECI(X_ECEF,t):
    w     =    7.292115*10**(-5)
    C     =    np.array([[np.cos(w*t), -np.sin(w*t), 0], 
                         [np.sin(w*t), np.cos(w*t), 0],
                         [0,0,1]])
    ECEF = np.matmul(np.linalg.inv(C),X_ECEF)
    return ECEF

def angleBetweenVectors(vector1, vector2):
    """
        Gets the angle between two vectors

        Inputs:
        vector1: numpy column vector
        vector2: numpy column vector

        Outputs:

        theta: the angle between the two vectors in degrees

    
    """
    dotProduct = np.dot(np.transpose(vector1), vector2) #Transpose because np has a cry when you try to dot product vectors
    normVector1 = np.linalg.norm(vector1)
    normVector2 = np.linalg.norm(vector2)
    
    # Ensure the denominator is not zero to avoid division by zero errors
    if normVector1 == 0 or normVector2 == 0:
        return None
    
    cos_theta = dotProduct / (normVector1 * normVector2)
    # Use arccosine to calculate the angle in radians
    theta = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    
    return np.rad2deg(np.linalg.norm(theta)) #Return a float instead of a numpy vector containing a single float

def directionalCosine(roll,pitch,yaw):
    """
        Function to return the directional cosine matrix
    """
    theta = pitch
    phi   = roll
    psi   = yaw 
    
    C = np.array([[np.cos(theta)*np.cos(psi), np.cos(theta)*np.sin(psi), -np.sin(theta)],
                    [np.sin(phi)*np.sin(theta)*np.cos(psi)  - np.cos(phi)*np.sin(psi), np.sin(phi)*np.sin(theta)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.sin(phi)*np.cos(theta)],
                    [np.cos(phi)*np.sin(theta)*np.cos(psi) + np.sin(phi)*np.sin(psi), np.cos(phi)*np.sin(theta)*np.sin(psi)- np.sin(phi)*np.cos(psi), np.cos(phi)*np.cos(theta)]])
    
    return C




def motionEquation(state):
    """
        Computes the derivative of the satellite state vector
        which contains the velocity and acceleration of the satellite.

        Inputs:

        state: 6x1 numpy vector of the following format
                    x
                    y
                    z
                    dx/dt
                    dy/dt
                    dz/dt

        Returns:
        dS: Derivate of the state vector in the following format
                    dx/dt
                    dy/dt
                    dz/dt
                    d^2x/dt^2
                    d^2y/dt^2
                    d^2z/dt^2
    """
    X = state[0:3]
    a  = calculateAcceleration(X)
    dS = np.concatenate((state[3:6],a), axis=0)
    return dS


def motionEquationJ2(state):
    """
        Computes the derivative of the satellite state vector
        accounting for J2 pertubations.

        Inputs:

        state: 6x1 numpy vector of the following format
                    x
                    y
                    z
                    dx/dt
                    dy/dt
                    dz/dt

        Returns:
        dS: Derivate of the state vector in the following format
                    dx/dt
                    dy/dt
                    dz/dt
                    d^2x/dt^2
                    d^2y/dt^2
                    d^2z/dt^2
    """
    X = state[0:3]  
    x = state[0]
    y = state[1]
    z = state[2]

    J2 = 1.08263*10**-3
    r = np.linalg.norm(X)
    a  = calculateAcceleration(X)
    p = 1.5*J2*u*R**2/r**4 * np.array([[(x/r)*(5*z**2/r**2-1)], [(y/r)*(5*z**2/r**2-1)], [(z/r)*(5*z**2/r**2-3)]])
    dS = np.concatenate((state[3:6],a), axis=0)
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


# Define the equation of motion for a Keplerian orbit
def keplerOrbit(params, t):
    """
        Determine the initial state of a satellite based on the 
        orbital parameters and time.

        Inputs:
        params       - Orbital parameters of the sattelite in this order:
                                e, i, w, RAA, M0, meanMotion
        t            - the time passed since epoch
    
    """
    i, RAA, e, w, M0, meanMotion, tsinceV = params
   
    #Get number of revolutions per second
    meanMotionSeconds   =   meanMotion/(24*60*60)
    #Calculate period in seconds
    T   =   1/meanMotionSeconds
    n = 2*np.pi/T
    M_e = M0 + n*t

    #Calculate Eccentric Anomaly
    E   =   calculateEccentricAnomaly(M_e,e)
    #Calculate True Anomaly

    trueAnom    =   2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2))

    #Calculate perigee
    a   =   (((T*np.sqrt(u))/(2*np.pi))**(2/3))
        
    #Calculate the distance from the focus
    r   =   (a*(1-e**2)/(1+e*np.cos(trueAnom)))

    #Calculate the Angular Momentum
    h   =    np.sqrt(a*u*(1-e**2))
    
    #Calculate radius of perigee
    r_p =  (h**2/u)* (1/(1+e))
    
    #Calculate semi-major axis
    r_a =  2*a - r_p

    #Get initial state of the satellite in perifocal frame
    state = calculateState(e,h, trueAnom)
    
    #Convert perifocal to ECI
    X, V = PFF2ECI(state[0],state[1],RAA,w,i)


    return X,V  # Return as a NumPy array





def calculateState(e,h,trueAnom):
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
    r = (h**2 / u)*(1/(1+e*np.cos(trueAnom))) * np.array([np.cos(trueAnom), np.sin(trueAnom), 0])
    v = (u/h)* np.array([-np.sin(trueAnom), e + np.cos(trueAnom), 0])
    return r,v


def calculateEccentricAnomaly(M_e,e):
    """
        Calcualtes the eccentric anomaly of a satellite

        Inputs:
        M_e: Mean Anomaly of a satellite
        e:   eccentricity of satellite's orbit

        Returns:
        E: eccentric anomaly

    """
    #Assume that eccentric anomaly is equal to mean anomaly
    E = M_e
    while (np.abs(M_e + e*np.sin(E) - E) > 0.00000001):
        E = M_e + e*np.sin(E)
    return E


if __name__ == "__main__":
            
    print(ECEF2GLLH(np.array([[-1052724.43586093],[5229852.40712444],[1966215.61789402]])))
    X_POLAR = np.array([[-5.10793795e+01],[-6.61054982e+00],[ 5.98658251e+08]])
    X_NEU = POLAR2CART(X_POLAR)
    groundstation = [[-26.6025329910276941], [118.52774285958469], [0]]    
    
    X_NEU =    [[-5932283.11672236],
        [ 3057896.77453712],
        [-5131477.2339245 ]]
    X_ECEF = LGDV2ECEF(X_NEU,groundstation)             
    print(X_NEU)
    X_NEU = ECEF2LGDV(X_ECEF, groundstation)
    print(X_NEU)



