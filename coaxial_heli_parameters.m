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

% Schedules
V_vals = [0.1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, ...
    70, 75, 80, 85, 90, 95, 100];
theta_f_trim = deg2rad([4.8, 4.8, 4.5, 4.3, 4, 3.8, 3.7, 3.5, 3.3, 3, 2.5, 2, 1.5, 1, ...
    0.5, 0 ,-0.5, -1 ,-1.5 ,-2.5, -2.5]);

% Atmospheric Data
g = 9.80665;
rho = 1.225;

% General Planform Data
mass = 3010;                % Su et al. 2023
W = mass*g;
Iyy = 25000;                % educated guess
CDS = 1.9;                  % Cd * S - guess from graph in report heli course assignment
% Omega = 35;              % [rad/s] Su et al. 2023
R = 5.49;                    % Su et al. 2023
diam = 2*R;
area = pi/4*diam^2;
N_u = 3;
N_l = 3;
lock = 6.57;                % Su et al. 2023
Vmax = 120;                 % [m/s]
h_l = 0.89;                 % Shaft spacing upper rotor - Su et al. 2023   
h_u = h_l + 0.77;           % Shaft spacing lower rotor - Su et al. 2023
hinge_offset_ratio = 0.1;
nu_b2 = 1 + 3/2*hinge_offset_ratio/(1-hinge_offset_ratio);

% Rotor RPM Scheduling

% Main rotor Airfoil Data
c = 0.257;                  % blade chord [m]
c_l_a = rad2deg(0.1);       % NACA 0012
% sigma = N_u*c*R/(pi*R^2);
sigma = 0.153;
twist_r = deg2rad(-10);
i_theta = deg2rad(-0);      % rotor shaft tilt

% Propeller Data
c_l_a_p = 5.7; C_l_p = 5.7;
c_d_p = 0.19;
omega_p = 207;  % 2300 rpm
R_p = 1.4;
r_p = 0:0.01:R_p;
N_p = 6;
c_p = 0.2;      % from manual calculation using provided solidity, number of blades, and prop radius
% sigma_p = N_p*c*R_p/(pi*R_p^2);
sigma_p = 0.142;
Kp = 0.4;       % Downwash factor from main rotors
l_p = 7.66;     % x location 
h_p = 0;        % z location
d_p = 0;        % y location
twist_p = deg2rad(-30);

% Horizontal Tail
S_h = 5.57;
l_h = 6.8;
C_l_h = 5.7;






