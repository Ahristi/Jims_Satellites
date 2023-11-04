import numpy as np
import matplotlib.pyplot as plt
import csv

R_E = 6378137

EastError = []
NorthError = []


def trueMappingError():
    # CSV FORMAT : X_eci_true, Y_eci_true, Z_eci_true, X_eci_observed, Y_eci_observed, Z_eci_observed
    with open('OBSERVATIONS.csv', mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0

        for row in csv_reader:
            xt = float(row[0])
            yt = float(row[1])
            zt = float(row[2])
            xo = float(row[3])
            yo = float(row[4])
            zo = float(row[5])


            # Find difference in azimuth angles of both vectors (for longitutde error)
            az1 = np.arctan2(yt, xt)
            az2 = np.arctan2(yo, xo)
            east_angle = az2 - az1
            # Convert angle to error in metres
            east_err =  east_angle/(2*np.pi) * 2*np.pi*R_E
            EastError.append(east_err)

            # Find difference in elevation angles of both vectors (for lattitude error)
            el1 = np.arcsin(zt/np.linalg.norm(np.array([xt, yt, zt])))
            el2 = np.arcsin(zo/np.linalg.norm(np.array([xo, yo, zo])))
            north_angle = el2 - el1
            # Convert angle to error in metres
            north_err =  north_angle/(2*np.pi) * 2*np.pi*R_E
            NorthError.append(north_err)

            line_count += 1


    SimNum = 10000
    sd = 13.12 

    # Simulate errors/noise in the East and North directions
    EastError_estimate = np.random.normal(0, sd, SimNum)
    NorthError_estimate = np.random.normal(0, sd, SimNum)



    # Plotting
    plt.figure(figsize=(8,6))
    #plt.scatter(EastError_estimate, NorthError_estimate, s=5, alpha=0.5, color='blue', label = 'Predicted Error Spread')
    plt.scatter(EastError, NorthError, s=5, alpha=0.5, color='red', label = 'True Error from Simulation') 
    plt.axhline(0, color='black', linewidth=0.5, linestyle='-')
    plt.axvline(0, color='black', linewidth=0.5, linestyle='-')
    plt.xlabel('East Error (m)')
    plt.ylabel('North Error (m)')
    plt.title('Nadir Mapping Error Using Sensor Fusion')
    #plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlim([-1000, 1000])
    plt.ylim([-1000, 1000])
    plt.show()

if __name__ == "__main__":
    trueMappingError()