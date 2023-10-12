"""
    Simulator

    Simulator class that propogates multiple satellite orbits.
    Performs satellite and ground station operations

"""
import numpy as np
from Satellite import *
from orbitalTransforms import *
import pyvista as pv
from pyvista import examples
import matplotlib.pyplot as plt

class Simulator:

    def __init__(self, satellites, groundstations):

        self.satellites = satellites            #Array of satellite objects to be simulated
        self.groundstations = groundstations    #Array of groundstations to be simulated


    def simulate(self, t0, t_end, h, f):
        """
            Uses a 4th order numerical propogator to simulate the satellite orbit.

            Inputs:
            t0      -   initial time
            t_end   -   end time
            h       -   time step
            f       -   orbital function for simulation
        """
        for sat in self.satellites:
            sat.times.append(t0)
        while t0 < t_end:
            for sat in self.satellites:
                k1 = h * f(sat.states[-1])
                k2 = h * f(sat.states[-1] + k1/2)
                k3 = h * f(sat.states[-1] + k2/2)
                k4 = h * f(sat.states[-1] + k3)


                #Propogate the satellite orbit
                state_new = sat.states[-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
                sat.states.append(state_new)
                sat.times.append(t0+h)
                currentECEF = ECI2ECEF([state_new[0], state_new[1], state_new[2]],t0+h)
                currentGLLH = ECEF2GLLH([currentECEF[0],currentECEF[1],currentECEF[2]])
                sat.ECEF.append(currentECEF)
                sat.GLLH.append(currentGLLH)

                #Service the satellite's routines
                sat.tick()
         
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

        # Create subplots with 3 rows and 1 column
        fig, axes = plt.subplots(3, 1, figsize=(8, 10))

        # Plot the first graph on the top
        axes[0].plot(sat.times, allAttitudes[:,0], color='blue')
        axes[0].set_title('Plot 1')
        axes[0].legend()

        # Plot the second graph in the middle
        axes[1].plot(sat.times, allAttitudes[:,1], color='green')
        axes[1].set_title('Plot 2')
        axes[1].legend()

        # Plot the third graph at the bottom
        axes[2].plot(sat.times, allAttitudes[:,2], color='red')
        axes[2].set_title('Plot 3')
        axes[2].legend()

        # Adjust spacing between the subplots
        plt.tight_layout()

        # Display the plots
        plt.show()

if __name__ == "__main__":
    sat = Satellite("ISS.txt", "ISS")
    sim = Simulator([sat], [])
    sim.simulate(0,24*60*60, 10, motionEquation)
    sim.showOrbit()
    sim.showAttitudes()



