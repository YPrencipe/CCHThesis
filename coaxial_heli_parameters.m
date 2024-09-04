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
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation Parameters
tau = 0.1;

% Schedules
V_vals = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, ...
    70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125];
V_vals_kts = V_vals.*1.944;
theta_f_trim = deg2rad([4.8, 4.8, 4.5, 4.3, 4, 3.8, 3.7, 3.5, 3.3, 3, 2.5, 2, 1.5, 1, ...
    0.5, 0 ,-0.5, -1 ,-1.5 ,-2.5, -2.5]);

% Atmospheric Data
g = 9.80665;
rho = 1.225;

% CG Positions
x_cg = -0.0164;
y_cg = 0;
z_cg = 0;

% General Planform Data
mass = 5500;                % Su et al. 2023
W = mass*g;
Ixx = 6800;
Ixz = 5000;
Iyy = 40000;                % educated guess 
Izz = 12000;
CDS = 1.9;                  % Cd * S - guess from graph in report heli course assignment
% Omega = 40;              % [rad/s] Su et al. 2023
Omega = 35;
R = 5.49;                    % Su et al. 2023
diam = 2*R;
area = pi/4*diam^2;
N_u = 4;
N_l = 4;
lock = 6.57;                % Su et al. 2023
Vmax = 120;                 % [m/s]
h_l = 0.89;                 % Shaft spacing upper rotor - Su et al. 2023   
h_u = h_l + 0.77;           % Shaft spacing lower rotor - Su et al. 2023
hinge_offset_ratio = 0.04;
% fictive_hinge_offset = 0.04;
fictive_hinge_offset = 0.445;
nu_b2 = 1 + 3/2*fictive_hinge_offset/(1-fictive_hinge_offset);
m_bl = 50;
N_b = N_u+N_l;
v_h_u = sqrt(mass/2*g/(2*rho*area));
v_h_l = sqrt(mass/2*g/(2*rho*area));
x_MR = -0.016;      % - means it is behind the c.g.
y_MR = 0;
CD_r = 0.025;


% Rotor RPM Scheduling

% Main rotor Airfoil Data
c = 0.257;                  % blade chord [m]
c_l_a = rad2deg(0.08);       % NACA 0012
% sigma = N_u*c*R/(pi*R^2);
sigma = 0.153;
twist_r = deg2rad(-10);
gamma_s = deg2rad(-0);      % rotor shaft tilt
K_b = 220500;   % [Nm/rad]
% K_b = 900000;   % [Nm/rad]
I_b = 1000;

% Propeller Data
c_l_a_p = 5.7; C_l_p = 5.7;
c_d_p = 0.19;
omega_p = 207;  % 2300 rpm
% R_p = 1.1;
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
k_u_pr = 0.9;
k_l_pr = 0.9;


% Horizontal Tail
S_ht = 1.197*1.7;  % if for example its 3, the helicopter becomes stable!
l_h = 6.8;
C_l_ht = rad2deg(0.04);
alfa_ht_0 = deg2rad(0);  % built-in horizontal tail incidence angle

% Elevator
S_e = 0.3 * S_ht;
C_l_e = rad2deg(0.04);

% Vertical Tail
S_vt = 1.197/2.3;
% S_vt = 2*2;
l_v = 6.8;
h_v = 0.5;
C_l_vt = rad2deg(0.04);
beta_vt_0 = deg2rad(0);  % built-in horizontal tail incidence angle

% Rudder
S_r = 0.1 * S_vt;
% S_r = 0.3 * S_vt;
C_l_r = rad2deg(0.04);

% Fuselage
Vol_fus = 12.2*pi*2^2;    % l*A cross section [circular]
K_fus = 0.7;
alfa_fus_m_0 = deg2rad(-1);


% Control Allocation
epsilon = 0.0000001; %small value






