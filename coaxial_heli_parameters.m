%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename:     coaxial_heli_parameters
% Project:      MSc Thesis
% Supervisor:   M.D. Pavel
% Author:       Ynias Prencipe 
% Student Nr.:  4777158
% 
% Description:  Define helicopter parameters used in the model
% 
% Issues to be solved:  - 
% 
% To Add:       - 
%
% Comments:     - The current values are from an H145D3 helicopter for easy
%               verification using an old helicopter course assignment
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Atmospheric Data
g = 9.80665;
rho = 1.225;

% General Planform Data
mass = 3010;                % Su et al. 2023
W = mass*g;
Iyy = 25000;                % educated guess
CDS = 1.9;                  % Cd * S - guess from graph in report heli course assignment
Omega = 38.94;              % [rad/s] Su et al. 2023
R = 5.4;                    % Su et al. 2023
diam = 2*R;
area = pi/4*diam^2;
N_u = 5;
N_l = 5;
sigma = N_u*c*R/(pi*R^2);
lock = 5.41;                % Su et al. 2023
Vmax = 120;                 % [m/s]
h_l = 0.89;                 % Shaft spacing upper rotor - Su et al. 2023   
h_u = h_l + 0.77;           % Shaft spacing lower rotor - Su et al. 2023

% Main rotor Airfoil Data
c = 0.257;                  % blade chord [m]
c_l_a = rad2deg(0.1);       % NACA 0012








