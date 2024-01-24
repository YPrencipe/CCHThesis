% Filename: compound_heli_parameters
% Course: Thesis
% Supervisor: M.D. Pavel
% Author: Ynias Prencipe 
% Student Number: 4777158
% Date of Delivery:

g = 9.80665;
rho = 1.225;

mass = 3010;            % Su et al. 2023
Iyy = 25000;            % guess
CDS = 1.9;              % Cd * S - guess from graph in report heli course assignment
Omega = 38.94;             % [rad/s] Su et al. 2023
R = 5.4;               % Su et al. 2023
diam = 2*R;
area = pi/4*diam^2;
N_u = 5;
N_l = 5;

c = 0.257;              % blade chord [m]
c_l_a = rad2deg(0.1);            % NACA 0012
         % Su et al. 2023
solidity = 0.127;
sigma = N_u*c*R/(pi*R^2);
                        
lock = 5.41;            % Su et al. 2023
Vmax = 120;

h_l = 0.89;
h_u = h_l + 0.77;

W = mass*g;




