import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from statistics import mean
import copy 
import pandas as pd
import pyarrow
import math
import scipy.io
from scipy.optimize import fsolve

theta_cps = 65 # heading of the CPS 
theta_wave = 220 # arbitraty wave heading which matches average wave heading from plotted data ('going to' convention)
theta_vert = 30
wave_velocity_horizontal = 0.5 #aritrary 0.5 m/s

#plotting rotation for each plot
elev = [0, 90]
azim = -90

#for headings from true north in clockwise east components are sin(theta), north components are cos(theta)
original_wave_directional_vec = np.array([math.sin(math.radians(theta_wave)), 
                                          math.cos(math.radians(theta_wave)), 
                                          0])

#defining the vector of the wave
original_wave_vector = wave_velocity_horizontal*original_wave_directional_vec

#calculate the magnitude of the normal and axial horizontal velocity components
Uw_norm = math.sin(math.radians(theta_cps - (theta_wave))) * wave_velocity_horizontal
Uw_ax = math.cos(math.radians(theta_cps - (theta_wave))) * wave_velocity_horizontal

#setting up the normal direction vector (why - 90 works??? double check this!!!)
Uw_norm_directional_vec = np.array([math.sin(math.radians(theta_cps - 90)), 
                                    math.cos(math.radians(theta_cps - 90)), 
                                    0])

#setting up the axial direction vector
Uw_ax_directional_vec = np.array([math.sin(math.radians(theta_cps)),
                                  math.cos(math.radians(theta_cps)), 
                                  0])

Uw_ax_vector = Uw_ax*Uw_ax_directional_vec

Uw_norm_vector = Uw_norm*Uw_norm_directional_vec

#check that the two velocity components add to 0.5
check_vel = math.sqrt((Uw_norm)**2 + (Uw_ax)**2)
if check_vel == wave_velocity_horizontal: 
    print('vel checked successfully')
else: 
    print(f'wave vel = {wave_velocity_horizontal}, check = {check_vel} calculation incorrect!')

#setting up the CPS direction vector (convention pos x = east, pos y = north, pos z = up vertically)
cps_direction_vector = np.array([-1*math.cos(math.radians(theta_vert)) * math.sin(math.radians(theta_cps)), 
                                 -1*math.cos(math.radians(theta_vert)) * math.cos(math.radians(theta_cps)), 
                                 math.sin(math.radians(theta_vert))])

#defining the start and end points of if the direction vector of the CPS were to be plotted
#this is so I can get the touch down point of the CPS and extract the 'starting point' for the direction vector to define the beam vector (so that 'up' is positive) 
CPS_points = pd.DataFrame()
CPS_points['x'] = [0, math.cos(math.radians(theta_vert))*math.sin(math.radians(theta_cps))]
CPS_points['y'] = [0, math.cos(math.radians(theta_vert))*math.cos(math.radians(theta_cps))]
CPS_points['z'] = [math.sin(math.radians(theta_vert)), 0]

CPS_starting_points = pd.DataFrame()
CPS_starting_points.loc[0, 'x'] = CPS_points.loc[1, 'x']
CPS_starting_points.loc[0, 'y'] = CPS_points.loc[1, 'y']
CPS_starting_points.loc[0, 'z'] = CPS_points.loc[1, 'z']

CPS_midpoint = [CPS_points['x'].mean(), CPS_points['y'].mean(), CPS_points['z'].mean()]

# # lets check that the beam direction vector is correct
# fig = plt.figure()
# ax = fig.add_subplot(projection = '3d')
# ax.plot3D([-1,1], [0,0], [0,0], color = 'black') #x axis
# ax.plot3D([0,0], [-1,1], [0,0], color = 'black') #y axis
# ax.plot3D([0,0], [0,0], [0,1], color = 'black') #z axis
# # plotting the CPS 'beam'
# ax.plot3D(CPS_points['x'], CPS_points['y'], CPS_points['z'], color = 'red', label = 'CPS')
# #plotting a unit vector of the cps direction vector
# ax.quiver(CPS_starting_points.loc[0, 'x'], CPS_starting_points.loc[0, 'y'], CPS_starting_points.loc[0, 'z'], cps_direction_vector[0], cps_direction_vector[1], cps_direction_vector[2], color = 'blue', label = 'unit vector?')
# # ax.quiver(0, 0, 0.5, 0.757443, 0.419857, 0, color = 'green', label = 'hard coded coords vec')
# ax.legend(loc = 'lower left', bbox_to_anchor = (1, 1))
# plt.title('Testing the direction vector direction')
# ax.view_init(elev = elev, azim = azim)
# plt.show()

for i in range(2):
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    ax.plot3D([-1,1], [0,0], [0,0], color = 'black') #x axis
    ax.plot3D([0,0], [-1,1], [0,0], color = 'black') #y axis
    ax.plot3D([0,0], [0,0], [0,1], color = 'black') #z axis
    # plotting the CPS 'beam'
    ax.quiver(CPS_starting_points.loc[0, 'x'], CPS_starting_points.loc[0, 'y'], CPS_starting_points.loc[0, 'z'], cps_direction_vector[0], cps_direction_vector[1], cps_direction_vector[2], color = 'red', label = 'CPS')
    
    #plotting the initial wave velocity
    ax.quiver(CPS_midpoint[0], CPS_midpoint[1], CPS_midpoint[2], original_wave_vector[0], original_wave_vector[1], original_wave_vector[2], color = 'blue', label = 'initial wave velocity')
    
    # plotting the axial velocity
    ax.quiver(CPS_midpoint[0], CPS_midpoint[1], CPS_midpoint[2], Uw_ax_vector[0], Uw_ax_vector[1], Uw_ax_vector[2], color = 'green', label = 'Uw_ax')
    
    # plotting the normal velocity
    ax.quiver(CPS_midpoint[0], CPS_midpoint[1], CPS_midpoint[2], Uw_norm_vector[0], Uw_norm_vector[1], Uw_norm_vector[2], color = 'orange', label = 'Uw_norm')
    
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')
    ax.set_title(f'CPS and velocities plot elev: {elev[i]}, azim: {azim}')
    ax.legend(loc = 'lower left', bbox_to_anchor = (1, 1))
    ax.view_init(elev = elev[i], azim = azim )
    plt.tight_layout()
    plt.show()

#####################
# now need to find the flow along the beam, this needs to be done using the dot product
#####################

# projecting the axial (horizontally parallel flow) onto the CPS to get the magnitude of the velocity along the CPS
vel_along_CPS_magnitude = np.dot(Uw_ax_vector, cps_direction_vector)

#Create the vector of the flow along the CPS for plotting
vel_along_CPS_vector = vel_along_CPS_magnitude * cps_direction_vector

# # calculating the flow along the CPS using derived formula
# vel_ax_CPS = Uw_ax*math.cos(math.radians(theta_vert))

# vel_ax_vector = vel_ax_CPS*cps_direction_vector
print(np.cross(vel_along_CPS_vector, cps_direction_vector))

for i in range(2):
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    ax.plot3D([-1,1], [0,0], [0,0], color = 'black') #x axis
    ax.plot3D([0,0], [-1,1], [0,0], color = 'black') #y axis
    ax.plot3D([0,0], [0,0], [0,1], color = 'black') #z axis
    # plotting the CPS 'beam'
    ax.quiver(CPS_starting_points.loc[0, 'x'], CPS_starting_points.loc[0, 'y'], CPS_starting_points.loc[0, 'z'], cps_direction_vector[0], cps_direction_vector[1], cps_direction_vector[2], color = 'red', label = 'CPS')
    #plotting the initial wave velocity
    ax.quiver(CPS_midpoint[0], CPS_midpoint[1], CPS_midpoint[2], original_wave_vector[0], original_wave_vector[1], original_wave_vector[2], color = 'blue', label = 'initial wave velocity')
    # plotting the axial velocity
    ax.quiver(CPS_midpoint[0], CPS_midpoint[1], CPS_midpoint[2], Uw_ax_vector[0], Uw_ax_vector[1], Uw_ax_vector[2], color = 'green', label = 'Uw_ax')
    # plotting the normal velocity
    ax.quiver(CPS_midpoint[0], CPS_midpoint[1], CPS_midpoint[2], Uw_norm_vector[0], Uw_norm_vector[1], Uw_norm_vector[2], color = 'orange', label = 'Uw_norm')
    
    #plotting the flow along the cps
    ax.quiver(CPS_midpoint[0], CPS_midpoint[1], CPS_midpoint[2], vel_along_CPS_vector[0], vel_along_CPS_vector[1], vel_along_CPS_vector[2], color = 'cyan', label = 'vertical axial flow component')
    
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')
    ax.set_title(f'checking flow along the CPS elev:{elev[i]}, azim:{azim}')
    ax.set_box_aspect([1,1,1])
    ax.legend(loc = 'lower left', bbox_to_anchor = (1, 1))
    ax.view_init(elev = elev[i], azim = azim)
    # ax.set_box_aspect([0.5,0.5,0.5])
    plt.tight_layout()
    plt.show()

#### now I need to plot and find the vertical velocity component of the flow along the beam
V_vert = vel_along_CPS_magnitude*math.sin(math.radians(theta_vert))
vertical_direction_vector = np.array([0,0,1])
V_vert_vector = V_vert*vertical_direction_vector

for i in range(2):
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    ax.plot3D([-1,1], [0,0], [0,0], color = 'black') #x axis
    ax.plot3D([0,0], [-1,1], [0,0], color = 'black') #y axis
    ax.plot3D([0,0], [0,0], [0,1], color = 'black') #z axis
    # plotting the CPS 'beam'
    ax.quiver(CPS_starting_points.loc[0, 'x'], CPS_starting_points.loc[0, 'y'], CPS_starting_points.loc[0, 'z'], cps_direction_vector[0], cps_direction_vector[1], cps_direction_vector[2], color = 'red', label = 'CPS')
    #plotting the initial wave velocity
    ax.quiver(CPS_midpoint[0], CPS_midpoint[1], CPS_midpoint[2], original_wave_vector[0], original_wave_vector[1], original_wave_vector[2], color = 'blue', label = 'initial wave velocity')
    # plotting the axial velocity
    ax.quiver(CPS_midpoint[0], CPS_midpoint[1], CPS_midpoint[2], Uw_ax_vector[0], Uw_ax_vector[1], Uw_ax_vector[2], color = 'green', label = 'Uw_ax')
    # plotting the normal velocity
    ax.quiver(CPS_midpoint[0], CPS_midpoint[1], CPS_midpoint[2], Uw_norm_vector[0], Uw_norm_vector[1], Uw_norm_vector[2], color = 'orange', label = 'Uw_norm')
    
    #plotting the flow along the cps
    ax.quiver(CPS_midpoint[0], CPS_midpoint[1], CPS_midpoint[2], vel_along_CPS_vector[0], vel_along_CPS_vector[1], vel_along_CPS_vector[2], color = 'cyan', label = 'vertical axial flow component')
    
    #plotting the vertical flow
    ax.quiver(CPS_midpoint[0], CPS_midpoint[1], CPS_midpoint[2], V_vert_vector[0], V_vert_vector[1], V_vert_vector[2], color = 'purple', label = 'vertical velocity component')
    
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')
    ax.set_title(f'checking flow along the CPS elev:{elev[i]}, azim:{azim}')
    ax.set_box_aspect([1,1,1])
    ax.legend(loc = 'lower left', bbox_to_anchor = (1, 1))
    ax.view_init(elev = elev[i], azim = azim)
    # ax.set_box_aspect([0.5,0.5,0.5])
    plt.tight_layout()
    plt.show()