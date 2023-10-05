import numpy as np


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



def PFF2ECEF(X,V,O,w,i):
    """
        Converts satellite state from 2D perifocal frame to ECEF

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

     az    =    np.deg2rad(az)
     el    =    np.deg2rad(el)

     x     =    r * np.cos(az)*np.cos(el)
     y     =    r * np.sin(az)*np.cos(el)
     z     =    r * np.sin(el)

     return np.array([x,y,z])


def ECEF2LGDV(X, ref):
    lat = np.deg2rad(ref[0][0])
    long = np.deg2rad(ref[1][0])
    R = ref[2][0]
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
    lat = np.deg2rad(ref[0][0])
    long = np.deg2rad(ref[1][0])
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



