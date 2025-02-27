"""
GENDATA FOR WAVE INTERATIONS PROJECT (05_05 , Seperation by 1.5 ridges, Fr 0.5)
25 February 2025
Klema, Matthew
Edited by: Flower, Sydney
"""

import numpy as np 
import xmitgcm
import matplotlib.pyplot as plt
import math

# -- Data Format Specifications -----------------------------------------

ieee = 'F'

# -- Constant & Calculated Parameters -----------------------------------

H0 = 10                          # full depth of the model (m)
h0 = 5                           # height of topography (from full depth H0) (m)
ht_hd = h0 / H0                  # depth ratio
L = 13                           # horizontal distance from start to location of topographic peak (m)
Ls = np.sqrt(L**2 + h0**2)       # slope length (m)
T = 0.5                          # wave period (hr)
om = 2 * np.pi / T / 3600        # tidal/forcing frequency (s^-1)
N0 = 1e-2                        # initial buoyancy frequency (s^-1)
u0 = 0.0149                     # open-boundary inflow velocity (m/s)
g = 9.81                         # gravitational acceleration (m/s^2)
f = 0                            # Coriolis frequency (s^-1)
s = np.sqrt((om**2 - f**2) / (N0**2 - om**2))  # characteristic wave slope
gamma = (h0 / L) / s             # relative slope steepness
lambda_x = (2 * H0) / s          # horizontal wavelength (m)
Fr = (u0 * np.pi * s) / (H0 * om)  # Froude number      
Fr_ex = (u0 * T * 3600) / (np.pi * Ls)  # Excursion Froude number

# -- Set Computational Domain/Grid Spacing ------------------------------

ny = 1                            # number of grid points in the lateral direction
nx = 2800                         # number of grid points in the horizontal direction
nz = 300                          # number of grid points in the vertical direction

Lx = 3.73 * lambda_x  +100          # total length of horizontal domain (m)
dy = Lx / nx                       # length of the domain in lateral direction (m)
  
dz = np.zeros(nz) + H0 / nz       # grid spacing in the vertical (m)
z = np.cumsum(dz)                 # vector of z values, bottom up (m)

dx = np.zeros(nx)                 # grid spacing in the horizontal (m)
for i in range(nx):
    dx[i] = Lx / (nx + 1)
x = np.cumsum(dx)                 # vector of x values, Left-to-Right (m)

# Writes dx & dz to binary files for input to the model
xmitgcm.utils.write_to_binary(dx.flatten(ieee), 'delXvar.bin', dtype='float64')
xmitgcm.utils.write_to_binary(dz.flatten(ieee), 'delZvar.bin', dtype='float64')

    
# -- Temperature Profile ------------------------------------------------

talpha = 2.0e-4                   # linear EOS thermal expansion coefficient (1/degree-C)
N2 = N0**2                        # square of the buoyancy frequency (s^-2)

# calculates the initial temperature profile
Tz = N2 / (g * talpha)
Tref = (Tz * z - np.mean(Tz * z)) * -1

# initialize temperature profile over array for entire domain
t = np.zeros((nx, ny, nz))
for k in range(nz):
    t[:, :, k] += Tref[k]

# Writes temp profile and grid to binary files for input to the model
xmitgcm.utils.write_to_binary(Tref.flatten(ieee), 'Tref.bin',dtype='float64')
xmitgcm.utils.write_to_binary(t.flatten(ieee), 'T0.bin',dtype='float64')


# -- Topography ---------------------------------------------------------

# generates topographic displacement heights above bed
h = ((h0 / 2) * (1 + np.cos((2 * np.pi * (x - 1.5 * lambda_x)) / (2 * L) - np.pi / 2)))/2


# Adjusting the height of the second curve to be half the height of the first
# h_second_curve = ((h0 / 2) * (1 + np.cos((2 * np.pi * (x - 0.5 * lambda_x)) / (2 * L) - np.pi / 2)))
h_second_curve = (h0 / 2) * (1 + np.cos((2 * np.pi * ((x - 37.5) - 0.5 * lambda_x)) / (2 * L) - np.pi / 2))

# converts to depth from surface to top of topography (m)
H = -H0 + h

# eliminates all topographic features except desired feature
shift = 3.679
# Topography Generation ---------------------------------------------------------

# Combine topographies: first curve full height, second curve half height
for i in range(len(x)):
    if x[i] < (152):
        H[i] = -H0
for i in range(len(x)):
    if x[i] > (178.81) and x[i] < (214.53):
         H[i] = -H0
for i in range(len(x)):
    if x[i] > (214.53) and x[i] < (240):
        H[i] = -H0  + h_second_curve[i]  
for i in range(len(x)):
    if x[i] > (240):
        H[i] = -H0 
  

# Writes topo to binary file used for input to the model
xmitgcm.utils.write_to_binary(H.flatten(ieee), 'topo.bin', dtype='float64')

# Writes topo to binary file used for input to the model
xmitgcm.utils.write_to_binary(H.flatten(ieee),'topo.bin',dtype='float64')
       
     
#---------Plots-----
plt.figure(1, figsize=(12, 10))

# subplot 1: position vs horizontal grid spacing
plt.subplot(2, 2, 1)
plt.plot(x, dx, 'b', label='$\Delta x$', linewidth=2)
plt.plot(z, dz, 'g', label='$\Delta z$', linewidth=2)
plt.xlim([0, 120])
plt.ylim([0, 0.5])
plt.grid(True)
plt.xlabel('x(m), z(m)', fontsize=18)
plt.ylabel('dx(m), dz(m)', fontsize=18)
plt.title('Grid Spacing', fontsize=18)
plt.legend(fontsize=16)

# subplot 2: plot of temperature profile
plt.subplot(2, 2, 2)
plt.plot(Tref, -z, 'r', linewidth=2)
plt.grid(True)
plt.xlabel('Temperature ($\pm$)', fontsize=18)
plt.ylabel('Depth (m)', fontsize=18)
plt.title('Temperature Profile', fontsize=18)

# subplot 3: plot of topograpfor i in range(len(x)):
plt.subplot2grid((2, 2), (1, 0), colspan=2)
plt.plot(x, H, 'k', linewidth=2)
plt.grid(True)
plt.xlabel('x (m)', fontsize=18)
plt.ylabel('Depth (m)', fontsize=18)
plt.title('Topography', fontsize=18)
plt.xlim([0, Lx])
plt.ylim([-H0, 0])
plt.xticks(np.arange(0, 301, 10))

plt.tight_layout()
plt.show()

# -- Printed Parameters -------------------------------------------------

ULx = Lx / u0                               # time for wave to traverse domain (s)
CFL_cond = 0.2                              # set CFL condition for stability
deltaTa = math.floor((CFL_cond/2) * dx[0] / u0)  # time step calculated for CFL criterion (s)
num_T1 = ULx / deltaTa                      # number of time steps needed for one wave to cross domain
num_T5 = 5 * num_T1                         # number of time steps needed for one wave to cross domain 5 times

viscAz = 1e-5              # model viscosity (vertical)
viscAh = 1e-5             # model viscosity (horizontal)
diffKhT = 1e-5             # model diffusion (horizontal)
diffKzT = 1e-5           # model diffusion (horizontal)

deltaTvH = math.floor(0.075 * (dx[0]**2 / viscAh))    # max time step based on horizontal viscosity 
deltaTvZ = math.floor(0.15 * (dz[0]**2 / viscAz))     # max time step based on vertical viscosity 
    
print('\t\t[General Simulation]')
print('------------------------------------------------------------------')
print(f'The topography-to-depth ratio is {ht_hd:.2f}')
print(f'The initial Froude number of the simulation is {Fr:.3f}')
print(f'The excursion Froude number of the simulation is {Fr_ex:.2f}')
print(f'The total length of the domain is {Lx:.0f} meters')

gamma = round(gamma, 2)               # rounds gamma to the nearest hundredth
if gamma == 1.0 or gamma == 1.00:
    print(f'The relative slope steepness is {gamma:.2f}, slope is critical')
elif gamma > 1:
    print(f'The relative slope steepness is {gamma:.2f}, slope is supercritical')
else:
    print(f'The relative slope steepness is {gamma:.2f}, slope is subcritical')
print(f'If wave has a velocity of {u0:.3f} m/s ({u0*100:.2f} cm/s) and the domain is {Lx:.0f} m\n it will take {ULx:.0f} seconds ({ULx/3600:.2f} hours) to traverse the domain.\n')

print('\t\t[Modify in "obcs_calc.F"]')
print('------------------------------------------------------------------')
print(f'The buoyancy frequency of the simulation, N, is {N0:.1E} s^-1 (N^2 = {N0**2:.2E} s^-2)')
print(f"The model depth is {H0:.0f} meters so twice the model depth 'DH={2*H0:.1f} meters'")
print(f'Tidal/forcing frequency (omega) is {om:.5f} rad/s^-1 (obTimeScale = {(2*math.pi)/om:.0f}s/{(2*math.pi)/om/3600:.2f}hr)')
print(f'The open-boundary inflow velocity, Uinflow, is {u0*100:.2f} cm/s ({u0:.3f} m/s)\n')

print('\t\t[Modify in "SIZE.h"]')
print('------------------------------------------------------------------')
print(f'The number of grid points in the horizontal, nx={nx:.0f} (dx={dx[0]:.2f}m)')
print(f'The number of grid points in the vertical, nz={nz:.0f} (dz={dz[0]:.2f}m)\n')

print('\t\t[Modify in "data"]')
print('------------------------------------------------------------------')
print(f"Variable 'sref' should be changed to a value for every cell, 'sref={nz:.0f}*35'")
print(f"The single grid cell size in the y-direction = 'dy = {dy:.2f}m'\n")

# -- README File --------------------------------------------------------

with open('README.txt', 'w') as fid:
    fid.write('\t\t[General Simulation]\n')
    fid.write('------------------------------------------------------------------\n')
    fid.write(f'The topography-to-depth ratio is {ht_hd:.2f}\n')
    fid.write(f'The initial Froude number of the simulation is {Fr:.3f}\n')
    fid.write(f'The excursion Froude number of the simulation is {Fr_ex:.2f}\n')
    fid.write(f'The total length of the domain is {Lx:.0f} meters\n')

    gamma = round(gamma, 2)               # rounds gamma to the nearest hundredth
    if gamma == 1.0 or gamma == 1.00:
        fid.write(f'The relative slope steepness is {gamma:.2f}, slope is critical\n')
    elif gamma > 1:
        fid.write(f'The relative slope steepness is {gamma:.2f}, slope is supercritical\n')
    else:
        fid.write(f'The relative slope steepness is {gamma:.2f}, slope is subcritical\n')
    fid.write(f'If wave has a velocity of {u0:.3f} m/s ({u0*100:.2f} cm/s) and the domain is {Lx:.0f} m\n it will take {ULx:.0f} seconds ({ULx/3600:.2f} hours) to traverse the domain.\n\n')

    fid.write('\t\t[Modify in "obcs_calc.F"]\n')
    fid.write('------------------------------------------------------------------\n')
    fid.write(f'The buoyancy frequency of the simulation, N, is {N0:.1E} s^-1 (N^2 = {N0**2:.2E} s^-2)\n')
    fid.write(f"The model depth is {H0:.0f} meters so twice the model depth 'DH={2*H0:.1f} meters'\n")
    fid.write(f'Tidal/forcing frequency (omega) is {om:.5f} rad/s^-1 (obTimeScale = {(2*math.pi)/om:.0f}s/{(2*math.pi)/om/3600:.2f}hr)\n')
    fid.write(f'The open-boundary inflow velocity, Uinflow, is {u0*100:.2f} cm/s ({u0:.3f} m/s)\n\n')
    fid.write(f'endTime= 22500')
  
    fid.write('\t\t[Modify in "SIZE.h"]\n')
    fid.write('------------------------------------------------------------------\n')
    fid.write(f'The number of grid points in the horizontal, nx={nx:.0f} (dx={dx[0]:.2f}m)\n')
    fid.write(f'The number of grid points in the vertical, nz={nz:.0f} (dz={dz[0]:.2f}m)\n\n')

    fid.write('\t\t[Modify in "data"]\n')
    fid.write('------------------------------------------------------------------\n')
    fid.write(f"Variable 'sref' should be changed to a value for every cell, 'sref={nz:.0f}*35'\n")
    fid.write(f"The single grid cell size in the y-direction = 'dely = {dy:.2f}m'\n\n")

    fid.write('\t\t[Values set in "data"]\n')
    fid.write('------------------------------------------------------------------\n')
    fid.write(f"Simulation vertical viscosity set at value 'viscAz = {viscAz:.2e}'\n")
    fid.write(f"Simulation horizontal viscosity set at value 'viscAh = {viscAh:.2e}'\n")
    fid.write(f"Simulation horizontal diffusion set at value 'diffKhT = {diffKhT:.2e}'\n")
    fid.write(f"Simulation vertical diffusion set at value 'diffKzT = {diffKzT:.2e}'\n")
    fid.write("Simulation time step set at value 'deltaT = 0.5s'\n")
