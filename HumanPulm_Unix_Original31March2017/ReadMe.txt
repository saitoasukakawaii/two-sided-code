This directory contains a combined c++ and Fortran90 code that simulates the results 
presented in Qureshi et al.(2014, Biomech Model Mechanobiol). The code can be executed
by running ./2SidedPulm command in a Mac terminal. (Check that an appropriate compiler
is available on your machine and the correct path of compiler is written in the Makefile).
The results can be produced both using an MRI measured Flow profile at the inlet of MPA
or a Pressure profile as an inflow boundary condition. Appropriate changes would be required
in sor06.h, arteries.h and in arteries.C, to switch between two inflow boundary conditions.
(Use cmd(F) to look for Q0 and Ps in arteries.C and arteries.h for change)

Files:

2SidedPulm: A C-Shell script to create and run an executable “sor06”
Makefile: Compiles c++ and fortran bits of the code together
sor06.h: Header file to control and initiate all basic parameter values
sor06.C: A main solver that calls every thing and then printout the data
arteries.h: Header file that declares class “Tube” and all its objects and functions,
            which are the defined and computed in arteries.C
arteries.C: Main computational code that computes all the variables associated 
            with every artery and vein. Calls tools.C, tools.h, new_match.f90, arteries.h
	    ludcump.h and nr3.h etc

gvPA_4096.dat: MRI measured inflow
p0_in4096.dat: Computed/simulated pressure at the inlet of MPA (x=0)
zero.dat: A file of zeros to handle outflow static pressure
Above data files are used as an inflow and out flow boundary conditions. Data in these files may
be interpolated to increase the number of time steps in the computational scheme. However, the
number of time steps should be an integer power of 2 (To meet the requirements of FFT in tools.f90)

new_match.f90: Fortran code that computes a 2 by 2 admittance matrix for each connected structured
tree in frequency domain.

Example.m and gnuplot.m:
Matlab files to plot the simulated data in the form of 3D and 1D plots.

