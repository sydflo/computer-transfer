		[General Simulation]
------------------------------------------------------------------
The topography-to-depth ratio is 0.50
The initial Froude number of the simulation is 0.500
The excursion Froude number of the simulation is 0.61
The total length of the domain is 300 meters
The relative slope steepness is 1.03, slope is supercritical
If wave has a velocity of 0.015 m/s (1.49 cm/s) and the domain is 300 m
 it will take 20152 seconds (5.60 hours) to traverse the domain.

		[Modify in "obcs_calc.F"]
------------------------------------------------------------------
The buoyancy frequency of the simulation, N, is 1.0E-02 s^-1 (N^2 = 1.00E-04 s^-2)
The model depth is 10 meters so twice the model depth 'DH=20.0 meters'
Tidal/forcing frequency (omega) is 0.00349 rad/s^-1 (obTimeScale = 1800s/0.50hr)
The open-boundary inflow velocity, Uinflow, is 1.49 cm/s (0.015 m/s)

endTime= 22500		[Modify in "SIZE.h"]
------------------------------------------------------------------
The number of grid points in the horizontal, nx=2800 (dx=0.11m)
The number of grid points in the vertical, nz=300 (dz=0.03m)

		[Modify in "data"]
------------------------------------------------------------------
Variable 'sref' should be changed to a value for every cell, 'sref=300*35'
The single grid cell size in the y-direction = 'dely = 0.11m'

		[Values set in "data"]
------------------------------------------------------------------
Simulation vertical viscosity set at value 'viscAz = 1.00e-05'
Simulation horizontal viscosity set at value 'viscAh = 1.00e-05'
Simulation horizontal diffusion set at value 'diffKhT = 1.00e-05'
Simulation vertical diffusion set at value 'diffKzT = 1.00e-05'
Simulation time step set at value 'deltaT = 0.5s'
