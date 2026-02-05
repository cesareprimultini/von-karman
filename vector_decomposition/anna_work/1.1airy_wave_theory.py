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

### note there may be issues with segments 3 and 4 as they are technically below the seabed
depth = 53.5
g = 9.81
start_easting = 213980
start_northing = 2696820
TDP = 15
CPS_heading = 65
num_segments = 5

met_input_filepath = r"P:\2198-CIP Feng Miao I (FM1) OWF\03 calcs\CTR 009 FLS flow amp\Stress and Fatigue Model\compiled_metocean_hindcast_2198_DRAFTNOTFINAL.csv"
catenary_profile_filepath = r"P:\2198-CIP Feng Miao I (FM1) OWF\03 calcs\CTR 009 FLS flow amp\Stress and Fatigue Model\Segments\catenary_profile_data.xlsx"

met_data = pd.read_csv(met_input_filepath)
met_data= met_data.dropna(axis = 1, how = 'all')

catenary_profile_data = pd.read_excel(catenary_profile_filepath, sheet_name = 'Sheet1')


#### THIS WAS ALL BASED ON OLD METHOD OF DEFINING SHALLOW AND DEEP WAVES RATHER THAN USING THE MORE SOPHISTIATED VERSION OF THE EQUATION
# #lets find if its deep or shallow
# met_data['depth toggle'] = depth/(g*((met_data['Tp'])**2))
# met_data.loc[met_data['depth toggle'] < 0.0025, 'shallow/deep'] = 'shallow'
# met_data.loc[met_data['depth toggle'] > 0.08, 'shallow/deep'] = 'deep'
# met_data.loc[(met_data['depth toggle'] > 0.0025) & (met_data['depth toggle'] < 0.08), 'shallow/deep'] = 'intermediate'

# # there appears to be some intermediate cases, lets plot the residual from the deep threshold
# met_data['residual deep'] = np.where(met_data['shallow/deep'] == 'intermediate', 0.08 - met_data['depth toggle'], np.nan)

# #lets plot the residuals
# residuals = pd.DataFrame()
# residuals['residual deep'] = met_data['residual deep']
# residuals = residuals.dropna(axis = 0, how = 'any')
# residuals = residuals.reset_index(drop = True)

# fig, ax = plt.subplots()
# plt.plot()
# ax.scatter(x = residuals.index, y = residuals['residual deep'], color = 'b', label = 'data points residual')
# ax.set_xlabel('index')
# ax.set_ylabel('residual')
# ax.set_title('deep data residuals')
# ax.axhline(y = (0.08 - 0.0025), color = 'red', label = 'shallow residual')
# plt.legend()
# plt.show()

# # #lets separate the df into 'deep' , 'shallow' and 'intermediate' 
# # df_deep = met_data[met_data[]] 

#### FROM HERE ONWARDS IS THE MORE SOPHISTICATED VERSION
#testing joel's airy wave vs mine
def blimit(met_df):
    k_lst = []
    L_lst = []
    for row in range(len(met_df)):
        T = met_df.loc[row, 'Tp']
        d = depth
        def dispersion(k):
            return (2*np.pi/T)**2-g*k*np.tanh(k*d)
        k = fsolve(dispersion,0)[0]
        print(f"k = {k}")
        L = 2*np.pi/k
        k_lst.append(k)
        L_lst.append(L)
    met_df['k'] = k_lst
    met_df['L'] = L_lst

blimit(met_data)

# #testing my method to validate both
# def find_L(met_df): 
#     k_lst = []
#     L_lst = []
#     for row in range(len(met_df)): 
#         T = met_df.loc[row, 'Tp']
#         d = depth
#         def dis_L(L): 
#             return (L - (g/(2*np.pi)) * (T**2) * np.tanh((2*np.pi*d)/L))
#         L = fsolve(dis_L,0)[0]
#         print(f'L = {L}')
#         k = (2*np.pi)/L
#         k_lst.append(k)
#         L_lst.append(L)
#     met_df['k_A'] = k_lst
#     met_df['L_A'] = L_lst    
    
# find_L(met_data)

# the data currently has negative values as it is supposed to have a scour hole? 
# lets check if there are any negatives and then if there are negatives then we need to move it above the seabed to get velocities
min_point = min(catenary_profile_data["Z' NB (m)"].to_list())
if min_point < 0: 
    catenary_profile_data["Z' NB (m)"] = catenary_profile_data["Z' NB (m)"] + abs(min_point)
#lets create the segments and then define their easting and northings 
#need to get rid of any values greater than 15m in the x axis (TDP is set as 15 m)
start_point_x = catenary_profile_data.loc[0, "X' BM (m)"]
start_point_z = catenary_profile_data.loc[0, "Z' NB (m)"]
condition = catenary_profile_data["X' BM (m)"] > TDP
catenary_profile_data = catenary_profile_data.drop(catenary_profile_data[condition].index)
end_index = catenary_profile_data.index[-1]
end_point = catenary_profile_data.loc[end_index, "X' BM (m)"]

fig, ax = plt.subplots()
ax.plot(catenary_profile_data["X' BM (m)"], catenary_profile_data["Z' NB (m)"], label = 'CPS')
ax.axhline(y = 0, color = 'yellow', label = 'seafloor')
ax.scatter(start_point_x, start_point_z, color = 'red', label = 'J tube exit')
ax.set_xlabel('horizontal displacement (m)')
ax.set_ylabel('vertical displacement (m)')
ax.set_title('Catenary profile')
plt.legend()
plt.show()

#have plotted and checked data
CPS_separated = {}
length_of_segments_x = (end_point - start_point_x)/num_segments
for i in range(num_segments): 
    temp_df = catenary_profile_data.copy()
    if i == 0: 
        start = 0
        fin = length_of_segments_x
    if i > 0: 
        start = length_of_segments_x*i
        fin = length_of_segments_x*(i+1)
    temp_df = temp_df[(temp_df["X' BM (m)"] >= start) & (temp_df["X' BM (m)"] <= fin)]
    CPS_separated[f"segment {i+1}"] = temp_df
      
    
segment_colours = ['red', 'orange', 'yellow', 'green', 'blue']
fig, ax = plt.subplots()
for i, (key, df) in enumerate(CPS_separated.items()):
    ax.plot(df["X' BM (m)"], df["Z' NB (m)"], color = segment_colours[i], linewidth = '4')
ax.axhline(y = 0, color = 'gray', label = 'seafloor')
ax.plot(catenary_profile_data["X' BM (m)"], catenary_profile_data["Z' NB (m)"], label = 'CPS', color = 'black', linewidth = 1)
ax.scatter(start_point_x, start_point_z, color = 'black', label = 'J tube exit', s = 100)
ax.set_xlabel('horizontal displacement (m)')
ax.set_ylabel('vertical displacement (m)')
ax.set_title('Catenary profile')
plt.legend()
plt.show()

CPS_segments = pd.DataFrame()
#need to now make the relevant data points and put into ABSTAB style sheet with all the relevant parameters (depth, heading, start and end eastings and northings ect)
for i, (key, df) in enumerate(CPS_separated.items()):
    df = df.reset_index(drop = True)
    #adding the x and z start and end point
    CPS_segments.loc[i,'X start'] = df.loc[0, "X' BM (m)"]
    CPS_segments.loc[i,'Z start'] = df.loc[0, "Z' NB (m)"]
    end_index = df.index[-1]
    CPS_segments.loc[i,'X end'] = df.loc[end_index, "X' BM (m)"]
    CPS_segments.loc[i,'Z end'] = df.loc[end_index, "Z' NB (m)"]
    #first find the equivalent cable length looking from above (top view)
    delta_z = CPS_segments.loc[i, 'Z end'] - CPS_segments.loc[i, 'Z start']
    delta_x = CPS_segments.loc[i, 'X end'] - CPS_segments.loc[i, 'X start']
    elevation = np.degrees(np.arctan2(delta_z, delta_x)) #calculates delta_z/delta_x
    #converting it to be measured vertical going clockwise
    # then - 90 to make it from the horizontal
    # elevation = ((90 - elevation) % 360) - 90
    #given elevation as negative (anticlockwise from z axis) we want to convert to depth below horizontal (how far is it dipping from the horizontal of the first point?)
    #converting to from vertical/north
    # elevation = (elevation + 360) % 360
    CPS_segments.loc[i, 'elevation'] = elevation
    
    if i == 0: 
        CPS_segments.loc[i,'easting start'] = start_easting
        CPS_segments.loc[i,'northing start'] = start_northing
        end_easting = start_easting + np.sin(np.radians(CPS_heading))*length_of_segments_x
        end_northing = start_northing + np.cos(np.radians(CPS_heading))*length_of_segments_x
        CPS_segments.loc[i,'easting end'] = end_easting
        CPS_segments.loc[i,'northing end'] = end_northing 
    else: 
        CPS_segments.loc[i, 'easting start'] = end_easting
        CPS_segments.loc[i, 'northing start'] = end_northing
        end_easting = end_easting + np.sin(np.radians(CPS_heading))*length_of_segments_x
        end_northing = end_easting + np.cos(np.radians(CPS_heading))*length_of_segments_x
        CPS_segments.loc[i,'easting end'] = end_easting
        CPS_segments.loc[i,'northing end'] = end_northing
        
    end_index = df.index[-1]
    CPS_segments.loc[i,'depth (midpoint)'] = depth - ((df.loc[0, "Z' NB (m)"] + df.loc[end_index, "Z' NB (m)"])/2)
    # CPS_segments['horiz disp'] = df["X' BM (m)"]
    CPS_segments.loc[i,'cable heading'] = CPS_heading
    segment_end_index = df.index[-1]
    # x1 = df.loc[0, "X' BM (m)"]
    # x2 = df.loc[segment_end_index, "X' BM (m)"]
    # z1 = df.loc[0, "Z' NB (m)"]
    # z2 =  df.loc[segment_end_index, "Z' NB (m)"]
    # print(f"x1: {x1}, x2: {x2}, z1: {z1}, z2: {z2}")
    CPS_segments.loc[i, 'cable length'] = math.sqrt(((delta_x)**2) + ((delta_z)**2))
    
## lets plot the found segments to double check that they make sense
plt.figure()
segment_colours = ['red', 'orange', 'yellow', 'limegreen', 'royalblue']
new_colours = ['maroon', 'chocolate', 'gold', 'darkgreen', 'midnightblue']
for i, (key, df) in enumerate(CPS_separated.items()): 
    plt.plot(df["X' BM (m)"], df["Z' NB (m)"], color = segment_colours[i], linewidth = '4')
    plt.plot([CPS_segments.loc[i, "X start"], CPS_segments.loc[i, 'X end']], [CPS_segments.loc[i, 'Z start'], CPS_segments.loc[i, 'Z end']], color = new_colours[i], linewidth = 2)
plt.xlabel('horizontal displacement')
plt.ylabel('distance from seafloor')
plt.title('Checking derived straight line segments')

# # Now that k and L have been found the vertical and horizontal partical velocity at the seabad can be found
# #first 'A' needs to be found, which is the wave phase angle (need to consider 90 and 180 as one will have a max vertical and a max horizontal velocity)
airy_wave_dict = {}
for i in range(len(CPS_segments)):
    temp_met_data = met_data.copy()
    A = np.radians(0)
    temp_met_data['0 deg horiz Uw'] = ((np.pi * temp_met_data['Hs'] * np.cosh(temp_met_data['k']*((CPS_segments.loc[i, 'depth (midpoint)']*-1) + depth))) / (temp_met_data['Tp']*np.sinh(temp_met_data['k']*depth))) * (round(np.cos(A)))
    temp_met_data['0 deg vert Uw'] = ((np.pi * temp_met_data['Hs'] * np.sinh(temp_met_data['k'] * ((CPS_segments.loc[i, 'depth (midpoint)']*-1) + depth))) / (temp_met_data['Tp'] * np.sinh(temp_met_data['k'] * depth))) * (round(np.sin(A)))
    A = np.radians(90)
    temp_met_data['90 deg horiz Uw'] = ((np.pi * temp_met_data['Hs'] * np.cosh(temp_met_data['k']*((CPS_segments.loc[i, 'depth (midpoint)']*-1) + depth))) / (temp_met_data['Tp']*np.sinh(temp_met_data['k']*depth))) * (round(np.cos(A)))
    temp_met_data['90 deg vert Uw'] = ((np.pi * temp_met_data['Hs'] * np.sinh(temp_met_data['k'] * ((CPS_segments.loc[i, 'depth (midpoint)']*-1) + depth))) / (temp_met_data['Tp'] * np.sinh(temp_met_data['k'] * depth))) * (round(np.sin(A)))
    airy_wave_dict[f"segment {i}"] = temp_met_data
    
# met_data['horiz_Uw'] = ((np.pi*met_data['Hs']*np.cosh(met_data['k']*(depth + depth))) / (met_data['Tp']*np.sinh(met_data['k']*depth))) * np.cos(met_data['Angular freq'])
# met_data['vert__Uw'] = ((np.pi*met_data['Hs']*np.sinh(met_data['k']*(depth + depth))) / (met_data['Tp']*np.sinh(met_data['k']*depth))) * np.sin(met_data['Angular freq'])

#plot the relationship between hs and tp 
for i, (key, df) in enumerate(airy_wave_dict.items()): 
    if i == 0: 
        plt.figure()
        plt.scatter(df['Tp'], df['Hs'], color = 'g')
        plt.ylabel('Hs (m)')
        plt.xlabel('Tp (sec)')
        plt.title('Met Data Relationship between Hs and Tp')
        plt.show()
    if i != 0: 
        break


#lets test the relationship between tp and velocity
test_df = pd.DataFrame()
test_df['Tp'] = list(range(1,20,1))
test_df['Hs'] = 3.5
test_df['k'] = 0.001
test_df['d'] = 50
test_df['z'] = -45
A = np.radians(0)
test_df['0 deg horiz Uw'] = ((np.pi * test_df['Hs'] * np.cosh(test_df['k']*((test_df['z']*-1) + test_df['d']))) / (test_df['Tp']*np.sinh(test_df['k']*test_df['d']))) * (round(np.cos(A)))
test_df['0 deg vert Uw'] = ((np.pi * test_df['Hs'] * np.sinh(test_df['k'] * ((test_df['z']*-1) + test_df['d']))) / (test_df['Tp'] * np.sinh(test_df['k'] * test_df['d']))) * (round(np.sin(A)))
A = np.radians(90)
test_df['90 deg horiz Uw'] = ((np.pi * test_df['Hs'] * np.cosh(test_df['k']*((test_df['z']*-1) + test_df['d']))) / (test_df['Tp']*np.sinh(test_df['k']*test_df['d']))) * (round(np.cos(A)))
test_df['90 deg vert Uw'] = ((np.pi * test_df['Hs'] * np.sinh(test_df['k'] * ((test_df['z']*-1) + test_df['d']))) / (test_df['Tp'] * np.sinh(test_df['k'] * test_df['d']))) * (round(np.sin(A)))

plt.figure()
plt.scatter(test_df['Tp'], test_df['0 deg horiz Uw'], color = 'red', label = 'horizontal velocity')
plt.scatter(test_df['Tp'], test_df['0 deg vert Uw'], color = 'b', label = 'Vertical velocity')
plt.xlabel('Tp (sec)')
plt.ylabel('Velocity (m/s)')
plt.title('Testing the relationship between Tp and velocity for 0 degrees')
plt.legend()
plt.show()

plt.figure()
plt.scatter(test_df['Tp'], test_df['90 deg horiz Uw'], color = 'red', label = 'horizontal velocity')
plt.scatter(test_df['Tp'], test_df['90 deg vert Uw'], color = 'b', label = 'Vertical velocity')
plt.xlabel('Tp (sec)')
plt.ylabel('Velocity (m/s)')
plt.title('Testing the relationship between Tp and velocity for 90 degrees')
plt.legend()
plt.show()

# now lets plot the horizontal and vertical components against the tp and hs values to check that they are making sense
for i, (key, df) in enumerate(airy_wave_dict.items()):
    fig, ax = plt.subplots(2,2)
    #first lets plot velocity vs Hs for 0 deg 
    ax[0,0].scatter(df['Hs'], df['0 deg horiz Uw'], color = 'red', label = 'Horizontal Velocity')
    ax[0,0].scatter(df['Hs'], df['0 deg vert Uw'], color = 'b', label = 'Vertical Velocity')
    ax[0,0].set_xlabel('Hs (m)')
    ax[0,0].set_ylabel('Velocity (m/s)')
    ax[0,0].set_title('Hs vs Uw phase angle = 0 deg')
    # ax[0,0].set_legend()
    # plot velocity vs Hs for 90 deg
    ax[0,1].scatter(df['Hs'], df['90 deg horiz Uw'], color = 'red', label = 'Horizontal Velocity')
    ax[0,1].scatter(df['Hs'], df['90 deg vert Uw'], color = 'b', label = 'Vertical Velocity')
    ax[0,1].set_xlabel('Hs (m)')
    ax[0,1].set_ylabel('Velocity (m/s)')
    ax[0,1].set_title('Hs vs Uw phase angle = 90 deg')
    # ax[0,1].set_legend()
    #plot velocity vs Tp for 0 deg 
    ax[1,0].scatter(df['Tp'], df['0 deg horiz Uw'], color = 'red', label = 'Horizontal Velocity')
    ax[1,0].scatter(df['Tp'], df['0 deg vert Uw'], color= 'b', label = 'Vertical Velocity')
    ax[1,0].set_xlabel('Tp (sec)')
    ax[1,0].set_ylabel('Velocity (m/s)')
    ax[1,0].set_title('Tp vs Uw phase angle = 0 deg')
    #plot velocity vs Tp for 90 deg
    ax[1,1].scatter(df['Tp'], df['90 deg horiz Uw'], color = 'red', label = 'Horizontal Velocity')
    ax[1,1].scatter(df['Tp'], df['90 deg vert Uw'], color = 'b', label = 'Vertical Velocity')
    ax[1,1].set_xlabel('Tp (sec)')
    ax[1,1].set_ylabel('Velocity (m/s)')
    ax[1,1].set_title('Tp vs Uw, phase angle = 90 deg')
    plt.suptitle(f'segment {i}')
    plt.legend()
    plt.tight_layout()
    plt.show()