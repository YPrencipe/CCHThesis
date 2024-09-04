README FILE COMPOUND COAXIAL HELICOPTER MATLAB FILES
------------------------------------------------------

Thank you for reading the README file for the code of the MATLAB files for the Compound Coaxial Helicopter thesis.

This README serves to guide a future student, continuing the work of the thesis, through the various files used to construct the results.


Main Files
------------
- coaxial_heli_parameters.m
	- Defines the helicopter parameters used in the model

- f_xk6.m
	- Main 6dof function used for trim and simulation --> calculates the inflow, flapping dynamics,
	forces and moments, and equations of motion
	- f_xk.m is a 3DoF version of the EOM's

- trim_6dof.m
	- Trimming using Newton-Rhapson and linearisation using central-difference
	- Linearisation starts at line 302
	- try_3_trim_3dof.m is a 3DoF version of the trimming

- nonLinSim_6dof.m
	- Nonlinear simulation using either longitudinal cyclic or elevator

- EMFController6dofSimulationBobAccellDecell.m
	- Simulates 3 potential manoeuvres: 3-2-1-1, Bob up/down w accel decell, or step input
	- Implements the EMF and PID controller architecture mathematically, 
	- Control allocation (line 464-546)
	- Actuator constraints (line 550-625)
	- Note that the gains are tuned for bob up/down w accel decell, so other manoeuvres might be unstable, 	especially when changing planform parameters

- CCH_PID_tuned_v1.slx
	- Simulink file of the control system using EMF and PID using 6dof linear model
	!!! RUN initLinCCH6dof.m in order to load the correct parameters, otherwise the Simulink simulation 	   	    won't work !!!

- reportPlotting.m
	- uses various stored .mat files for neat report plotting
	- 




Contact Details
----------------
Name:			Ynias J.E. Prencipe 
Prefered e-mail: 	ynias.prencipe1999@gmail.com
Phone Number:		+32 487 53 72 48

