# 3D-NS-FSI
3D Navier-Stokes, Fluid Structure Interaction solver, using a sharp interface Immersed Boundary Method

The goal for this code is to simulate complex Fluid-Structure Interaction (FSI) flow problems, by solving the Navier-Stokes equations for compressible flow.
The N-S eqs are solved numerically using a 2nd order Finite Difference Method (FDM) and marching in time. 

## Dependencies
- libconfig++	- For reading config files
- SFML		    - Was meant to be used for output. Now only used for timing. Will probably remove this dependency.
- Eigen		    - For solving 8x8 linear system for trilinear interpolation.

## Input
All input is given via the ConfigFile. All variables from this file are imported using libconfig++, into a class named ConfigSettings. 
If new inputs are introduced, they should be added to ConfigFile, added as member in the ConfigSettings class, and added to the method (of ConfigSettings) that reads each setting from ConfigFile.

