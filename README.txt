This README serves to guide the reader through the various files used to construct the results
for the Compound Coaxial Helicopter Modelling, Control and HQA Thesis in cooperation with NLR.

Main Files
------------
- coaxial_3dof_trim.m
	- Description: Trim using the Newton-Rhapson algorithm for the coaxial helicopter, no pusher prop or elevator added yet
	- Dependencies:
		- coaxial_heli_parameters.m: defines the helicopter parameters
		- f_xk.m: Function to determine the state vector using the small increments 
			of the trim vector during the Newton-Rhapson trim algorithm

- f_xk.m
	- Description: Function to determine the state vector using the small increments of the trim vector during the 
	Newton-Rhapson trim algorithm
	- Dependencies:
		- coaxial_heli_parameters.m: defines the helicopter parameters

- coaxial_heli_parameters.m
	- Description: Defines the helicopter parameters
	- Dependencies:
		NONE

Files Used for Manual Calculations
-----------------------------------


Data Files:
--------------



================================================================================================

Contact Details
----------------
Name:			Ynias J.E. Prencipe 
Student Number:		4777158
Student e-mail:		y.j.e.prencipe@student.tudelft.nl
Prefered e-mail: 	ynias.prencipe1999@gmail.com
Phone Number:		+32 487 53 72 48
