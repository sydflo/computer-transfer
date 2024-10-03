"""
Created on Tue Jul 16 14:26:34 2024

@author: drmasha
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import MITgcmutils as mitgcm
from MITgcmutils import mds
import matplotlib as mpl


#%% Data Extraction
data_path = ('/home/drmasha/SPRING_RESERCH/RESULTS/DoubleRidge/025_05_/0.84/run')

XC = mds.rdmds('/home/drmasha/SPRING_RESERCH/RESULTS/DoubleRidge/025_05_/0.84/run/XC')
RC = np.squeeze(mds.rdmds('/home/drmasha/SPRING_RESERCH/RESULTS/DoubleRidge/025_05_/0.84/run/RC'))  # loads cell data

##To ensure XC and RC are !D arrays
XC=XC.flatten()
RC=RC.flatten(momentum_w = np.squeeze(mds.rdmds('/home/drmasha/SPRING_RESERCH/RESULTS/DoubleRidge/025_05_/0.84/run/momentumvars_w', file))
)
uo = 0.0251  # deceleration of open boundary velocity for simulation
file = 224400  # loads data at numbered step set
start = 100
final = 2800

# state-variable extraction at step set by 'file' variable
statevars = np.squeeze(mds.rdmds('/home/drmasha/SPRING_RESERCH/RESULTS/DoubleRidge/025_05_/0.84/run/statevars', file))

U = statevars[0, :, :].T  # zonal velocity
U[U == 0] = np.NaN
V = statevars[1, :, :].T  # meridional velocity
V[V == 0] = np.NaN
W = statevars[2, :, :].T  # vertical velocity
W[W == 0] = np.NaN
T = statevars[3, :, :].T  # temperature data
T[T == 0] = np.NaN

PHIHYD = statevars[5, :, :].T # hydrostatic pressure potential anomaly
PHIHYD[PHIHYD == 0] = np.NaN
PHINH = statevars[5, :, :].T  # non-hydrostatic pressure potential anomaly
PHINH[PHINH == 0] = np.NaN
dRHOdR = statevars[6, :, :].T  # stratification
dRHOdR[dRHOdR == 0] = np.NaN
RHOanom = statevars[6, :, :] .T # density anomaly (rho-rho_c)
RHOanom[RHOanom == 0] = np.NaN

# extract 2D data from the state - variable file at step set by 'file'
statevars2D = np.squeeze(mds.rdmds('/home/drmasha/SPRING_RESERCH/RESULTS/DoubleRidge/025_05_/0.84/run/statevars2d', file))
ETAN = statevars2D[0, :].T # surface height anomaly
PHIBOT = statevars2D[1, :].T # bottom pressure potential anomaly

# data extraction for zonal momentum quantities
momentum_u = np.squeeze(mds.rdmds('/home/drmasha/SPRING_RESERCH/RESULTS/DoubleRidge/025_05_/0.84/run/momentumvars_u', file))
TOTAL = momentum_u[0, :, :].T
TOTAL[TOTAL == 0] = np.NaN  # total zonal momentum tendency (m/s/day)
Um_diss = momentum_u[1, :, :].T
Um_diss[Um_diss == 0] = np.NaN  # momentum tendency from dissipation (m/s^2)
Um_adv = momentum_u[2, :, :].T
Um_adv[Um_adv == 0] = np.NaN  # momentum tendency from advection (m/s^2)
Um_dPI = momentum_u[3, :, :].T
# momentum tendency from pressure gradient (m/s^2)
Um_dPI[Um_dPI == 0] = np.NaN
Um_Ext = momentum_u[4, :, :].T
Um_Ext[Um_Ext == 0] = np.NaN  # momentum tendency from external forcing (m/s^2)
momKE = momentum_u[5, :, :].T
momKE[momKE == 0] = np.NaN  # kinetic energy from momentum eq (m^2/s^2)
AB_gU = momentum_u[6, :, :].T
AB_gU[AB_gU == 0] = np.NaN  # momentum tendency from Adams- Bashforth(m/s^2)

# data extraction for vertical momentum quantities
momentum_w = np.squeeze(mds.rdmds('/home/drmasha/SPRING_RESERCH/RESULTS/DoubleRidge/025_05_/0.84/run/momentumvars_w', file))
Wm_diss = momentum_w[3, :, :].T
Wm_diss[Wm_diss == 0] = np.NaN  # momentum tendency from dissipation(m/s^2)
Wm_adv = momentum_w[2, :, :].T
Wm_adv[Wm_adv == 0] = np.NaN  # momentum tendency from advection(m/s^2)
AB_gW = momentum_w[0, :, :].T
AB_gW[AB_gW == 0] = np.NaN  # momentum tendency from Adams-Bashforth (m/s^2)
ADVr_Th = momentum_w[1, :, :].T
# vertical advective flux advective flux of potential temp. (degC-m^3/s)
ADVr_Th[ADVr_Th == 0] = np.NaN

del momentum_u, momentum_w, statevars, statevars2D
# check if the array is not empty before normalizing
if np.isnan(T).all():
    raise ValueError('Temperature array contains only NaN values.')

#norm_T = mpl.colors.Normalize(vmin=np.nanmin(T), vmax=np.nanmax(T), clip=False)

X,Y=np.meshgrid(XC,RC)
print(T.shape)


# %% Figure 1 - Temperature over full domain at specified time

fig1, ax1 = plt.subplots(figsize=(14, 3.5))
#extent = [XC[0], XC[2800], -RC[0], -RC[200]]
mesh1 = ax1.pcolormesh(X,-Y, T, cmap='jet',shading='auto')
#contour=ax1.contour(XC,-RC,T,levels=10,color='k',linewidths=1)
ax1.set_xlabel('Distance (m)', fontsize=18)
ax1.set_ylabel('Depth (m)', fontsize=18)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.text(56, 0.5, '$\\mathbf{(a)}$', fontsize=24,
        verticalalignment='bottom', horizontalalignment='right')
c1 = plt.colorbar(mesh1, ax=ax1)
c1.set_label('$\\pm T$ $(^\\circ C)$', fontsize=18)
plt.set_cmap('jet')
fig1.set_size_inches(14, 3.5)
plt.xticks(np.arange(0, 200, 10))
plt.show()

fig1.savefig('Fig_1.esp',format='esp',dpi=400)
#fig1.savefig('~/Documents/dissertation/documents/images/Fig_crit_05_01.esp',format='esp',dpi=600) #600

# %% Figure 2 - Zonal velocity over full domain at specified time
fig2, ax2 = plt.subplots(figsize=(16, 4))
X,Y=np.meshgrid(XC[start:final],RC)
mesh2= ax2.pcolormesh(X[start:final,start:final],-Y[start:final,start:final],U[start:final, start:final], cmap='viridis',shading='auto')
ax2.set_xlabel('Distance(m)', fontsize=18)
ax2.set_ylabel('Depth(m)', fontsize=18)
c2 = plt.colorbar(mesh2, ax=ax2)
c2.set_label('$\\pm U_0$ $(m/s)$', fontsize=18)
mesh2.set_clim([-0.02, 0.02])
plt.show()
fig2.savefig('Fig_2.eps', format='eps', dpi=600)


# %% Figure 3 - Quiver plot at specified time step
fig3, ax3 = plt.subplots(figsize=(16, 5))
X,Y=np.meshgrid(XC[start:final],RC)
mesh3 = ax3.pcolormesh(X[start:final,start:final], -Y[start:final,start:final], T[start:final, start:final], cmap='jet',shading='nearest')
ax3.quiver(XC[start:final], -RC, W[:, start:final], U[:, start:final], scale=2)
fig3.set_size_inches(12, 3.75)
fig3.savefig('Fig_3.eps', format='eps', dpi=600)
plt.show()

# %% Figure 6 - Wave snapshots
t_w = 1800      # wave period (s)
dt = 1          # time step (s)
U0 = 0.0251      # velocity amplitude of forcing
start = 100     # starting vertical profile to not display entire domain
final = 2800    # ending vertical profile to not display entire domain
files = [3600, 7200, 10800, 14400, 18000, 21600]

fig6, axes6 = plt.subplots(3, 2, figsize=(16, 10))

for i, file in enumerate(files):
    ax6 = axes6.flat[i]
    t_ratio = (file * dt) / t_w

    statevars = np.squeeze(mds.rdmds('/home/drmasha/SPRING_RESERCH/RESULTS/DoubleRidge/025_05_/0.1/01 RESULTS/statevars', file))

    T = statevars[:, :, 3].T
    T[T == 0] = np.nan
    T_masked = np.ma.masked_where(np.isnan(T), T)

    mesh6 = ax6.pcolormesh(X[start:final,start:final], -Y[start:final,start:final], T_masked[start:final,start:final], cmap='jet',shading='nearest')
    
    ax6.set_xticks([])
    ax6.set_yticks([])

    label = f'{t_ratio:.2f}'
    ax6.text(50, 1.5, f'$(i)$ $t/T= ${label}', fontsize=15, transform=ax6.transAxes, va='top', ha='right')
