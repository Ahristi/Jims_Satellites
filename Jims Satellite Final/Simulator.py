"""
    Simulator

    Simulator class that propogates multiple satellite orbits.
    Performs satellite and ground station operations

"""
import numpy as np
from Satellite import *
from navSatellite import *
from orbitalTransforms import *
from groundStation import *
import pyvista as pv
from pyvista import examples
import matplotlib.pyplot as plt
from matplotlib.image import imread
from tqdm import tqdm

SUN_R = 149597870700 #1AU
SUN_W = -2*np.pi/(24*60*60) #Angular velocity of the sun


class Simulator:

    def __init__(self, satellites, gnssConstellation = [], groundstations = []):

        self.satellites        = satellites        #Array of satellite objects to be simulated
        self.gnssConstellation = gnssConstellation #Array of GNSS satellites to be simulated
        self.groundstations    = groundstations    #Array of groundstations to be simulated
        self.sunAngle          =   0               #Mean anomaly of the sun in ECI frame


    def simulate(self, t0, t_end, h, f):
        """
            Uses a 4th order numerical propogator to simulate the satellite orbit.

            Inputs:
            t0      -   initial time
            t_end   -   end time
            h       -   time step
            f       -   orbital function for simulation
        """
        #Make the progress bar
        pbar = tqdm(total=100) 
        pbar.set_description("Simulating ")

        #Connect the GNSS constellation with Jim's Constellation and the groundstations
        for sat in self.satellites:
            sat.times.append(t0)
            sat.connectGNSS(self.gnssConstellation, self.groundstations)
        for gSat in self.gnssConstellation:
            gSat.times.append(t0)


        #Begin propogating all satellite orbits
        while t0 < t_end:
            #Propogate the sun
            self.sunAngle += SUN_W*h
            #Simulate GNSS satellites
            for gSat in self.gnssConstellation:

                k1 = h * f(gSat.states[-1])
                k2 = h * f(gSat.states[-1] + k1/2)
                k3 = h * f(gSat.states[-1] + k2/2)
                k4 = h * f(gSat.states[-1] + k3)

                state_new = gSat.states[-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
                gSat.states.append(state_new)
                gSat.times.append(t0+h)
                currentECEF = ECI2ECEF([state_new[0], state_new[1], state_new[2]], gSat.tSinceVernal + t0+h)
                currentGLLH = ECEF2GLLH([currentECEF[0],currentECEF[1],currentECEF[2]])
                gSat.ECEF.append(currentECEF)
                gSat.GLLH.append(currentGLLH)
                gSat.sunPos = self.calculateSunPos()

                #Service the satellite's routines
                gSat.tick()
            #Simulate Jim's Constellation
            for sat in self.satellites:
                k1 = h * f(sat.states[-1])
                k2 = h * f(sat.states[-1] + k1/2)
                k3 = h * f(sat.states[-1] + k2/2)
                k4 = h * f(sat.states[-1] + k3)
                state_new = sat.states[-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
                sat.states.append(state_new)
                sat.times.append(t0+h)
                currentECEF = ECI2ECEF([state_new[0], state_new[1], state_new[2]], sat.tSinceVernal + t0+h)
                currentGLLH = ECEF2GLLH([currentECEF[0],currentECEF[1],currentECEF[2]])
                sat.ECEF.append(currentECEF)
                sat.GLLH.append(currentGLLH)
                sat.sunPos = self.calculateSunPos()

                #Service the satellite's routines
                sat.tick()
                
            #Update progress bar
            pbar.update(100*h/t_end)
            t0 += h

    def showOrbit(self):
        """
            Plots the satellite orbit using pyvista

            Assumes that the orbit has already been propogated using simulate()

            NOTE: This only plots the first satellite orbit at the moment        
        """

        title = "Simulated Satellite Orbit"
        sat = self.satellites[0]
        x = []
        y = []
        z = []
        #This is dodgy I know 
        for state in sat.states:
            x.append(state[0])
            y.append(state[1])
            z.append(state[2])

        mesh = examples.planets.load_earth(radius = 6378137)
        texture = examples.load_globe_texture()
        pl = pv.Plotter()
        scatter_points = np.column_stack((x, y, z))
        pl.add_title(title, font_size=18, color='cyan', font=None, shadow=False)
        scatter1 = pl.add_points(scatter_points, color="deeppink", point_size=2)
      
        image_path = examples.planets.download_stars_sky_background(
            load=False
        )
        pl.add_background_image(image_path)
        _ = pl.add_mesh(mesh, texture=texture)
        pl.show()


    def showAttitudes(self):
        """
            Plots the satellite attitudes from a simulation

            Assumes that the orbit has already been propogated using simulate()

            NOTE: This only plots the first satellite attitudes at the moment        
        """
        sat = self.satellites[0]
        allAttitudes = np.array(sat.attitudes)
        estimatedAttitudes = np.array(sat.ADCS.estimatedAttitudes)
        points = self.getStateTimes(sat.times,sat.imaging)

        # Create subplots with 3 rows and 1 column
        fig, axes = plt.subplots(3, 1, figsize=(8, 10))

        # Plot roll
        axes[0].plot(sat.times, allAttitudes[:,0], color='blue')
        axes[0].plot(sat.times, estimatedAttitudes[:,0], color='red')
        axes[0].axvspan(points[0], points[1], color='orange', alpha=0.5, label = "Imaging")
        axes[0].legend()
        axes[0].set_title('Roll')
        axes[0].set_xlabel('Time (s)')
        axes[0].set_ylabel('Roll (Rad)')



        # Plot pitch
        axes[1].plot(sat.times, allAttitudes[:,1], color='blue')
        axes[1].plot(sat.times, estimatedAttitudes[:,1], color='red')
        axes[1].axvspan(points[0], points[1], color='orange', alpha=0.5, label = "Imaging")
        axes[1].legend()  
        axes[1].set_title('Pitch')
        axes[1].set_xlabel('Time (s)')
        axes[1].set_ylabel('Pitch (Rad)')
     

        # Plot yaw
        axes[2].plot(sat.times, allAttitudes[:,2], color='blue')
        axes[2].plot(sat.times, estimatedAttitudes[:,2], color='red')
        axes[2].axvspan(points[0], points[1], color='orange', alpha=0.5, label = "Imaging")
        axes[2].legend()     
        axes[2].set_title('Yaw')
        axes[2].set_xlabel('Time (s)')
        axes[2].set_ylabel('Yaw (Rad)')
    

        # Adjust spacing between the subplots
        plt.tight_layout()

        # Display the plots
        plt.show()

    def showCharges(self):
        """
            Plots the state of charge of the satellite's battery throughout the orbit

            Assumes that the orbit has already been propogated using simulate()

            NOTE: This only plots the first satellite charge at the moment    
        """
        sat = self.satellites[0]
        fig, ax = plt.subplots()
        ax.plot(sat.times, sat.EPS.charges, color = "r", label = "Charge")

        #Show Eclipses
        points = self.getStateTimes(sat.times, sat.eclipses)
        for i in range(0, len(points), 2):
            if (i+1 < len(points)):
                if (i == 0):
                    ax.axvspan(points[i], points[i+1], color='blue', alpha=0.5, label = "Eclipse")
                else:
                    ax.axvspan(points[i], points[i+1], color='blue', alpha=0.5)

        ax.set_ylim(bottom=0) 
        ax.set_ylim(top=23) 
        ax.legend()     
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Battery Charge (Whr)")
        ax.set_title("Battery charge over mission")
        plt.show()

    def showGroundTrack(self):
        """
            Plots the ground track of one satellites orbit

            Assumes that the orbit has already been propogated using simulate()

            NOTE: This only plots the first satellite ground track
        
        """
        #Plot the ground track in a new window 
        sat = self.satellites[0]
        GLLH = np.array(sat.GLLH)

        fig, ax2 = plt.subplots()
        imagePath = 'MapOfEarth.jpg'
        image = imread(imagePath)
        ax2.grid(alpha=0.5)
        ax2.imshow(image, extent = [-180,180,-90,90])
        ax2.scatter(GLLH[:,1], GLLH[:,0], s = 0.01, color ="red", zorder = 2)
        ax2.set_title("Satellite Ground Track")
        ax2.set_xlabel("Longitude (deg)")
        ax2.set_ylabel("Latitude (deg)")


    def calculateSunPos(self):
        """
            Calculates the position of the sun in ECI 
            frame based on the current sun angle. 

            I probably need to check this because the sun
            might orbit the wrong way lol.   
        """
        x = SUN_R*np.cos(self.sunAngle)
        y = SUN_R*np.sin(self.sunAngle)
        z = 0

        return np.array([x,y,z])


    def getStateTimes(self, times, eclipses):
        """
            Get the start and end times of a boolean array.
            Used to determine the start and end times of eclipses and imaging.      
        """
        points = []
        inEclipse = False #Bool for if we have already detect this eclipse or not


        for i in range(len(eclipses)):
            if eclipses[i]: 
                if not inEclipse:
                    points.append(times[i])
                    inEclipse = True
            elif inEclipse:
                points.append(times[i])
                inEclipse = False
        return points
    def showGroundTrack(self):
        """
            Plots the ground track of one satellites orbit

            Assumes that the orbit has already been propogated using simulate()

            NOTE: This only plots the first satellite ground track
        
        """
        #Plot the ground track in a new window 
        sat = self.satellites[0]
        GLLH = np.array(sat.GLLH)

        fig, ax2 = plt.subplots()
        imagePath = 'MapOfEarth.jpg'
        image = imread(imagePath)
        ax2.grid(alpha=0.5)
        ax2.imshow(image, extent = [-180,180,-90,90])
        ax2.scatter(GLLH[:,1], GLLH[:,0], s = 0.01, color ="red", zorder = 2)
        ax2.set_title("Satellite Ground Track")
        ax2.set_xlabel("Longitude (deg)")
        ax2.set_ylabel("Latitude (deg)")

    def showPositions(self):
        print("Showing Positions:")
        sat = self.satellites[0]  # Assuming you want to plot the positions of the first satellite
        sat_positions = np.array(sat.GNSS.positionEstimates)  
        x = sat_positions[:, 0]  # Extract X coordinates
        y = sat_positions[:, 1]  # Extract Y coordinates
        z = sat_positions[:, 2]  # Extract Z coordinates

        # Create a 3D scatter plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z, c='b', marker='o', label='Satellite Positions')

        # Label the axes
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')

        # Add a legend
        ax.legend()

        # Show the plot
        plt.show()

    

     
if __name__ == "__main__":

    print("Generating Satellites...")
    sat1 = Satellite("Satellites/sat1.txt", "SAT1")
    nav1 = navSatellite("navSatellites/GNSS1.txt", "GNSS1")
    nav2 = navSatellite("navSatellites/GNSS2.txt", "GNSS2")
    nav3 = navSatellite("navSatellites/GNSS3.txt", "GNSS3")
    nav4 = navSatellite("navSatellites/GNSS4.txt", "GNSS4")
    gs = groundStation(-32.9986, 148.2621, 415)
    print("Satellites created")
    sim = Simulator([sat1], [nav1,nav2,nav3,nav4], [gs])
    sim.simulate(0,6*60*60,5, motionEquation)
    sim.showPositions()
    sim.showGroundTrack()
    sim.showOrbit() 
    sim.showAttitudes()
    sim.showCharges()




