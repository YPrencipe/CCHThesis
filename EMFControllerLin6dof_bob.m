%% LINEAR Simulation

clear;clc;close all;

% Load Newest Trim States and Helicopter Parameters
load('trim_saved_6dof.mat'); load('TrimLinSave_6dof.mat'); 
load('V_transit_save.mat'); 
load('ref_3211'); 
load('theta_cmd_3211'); load('phi_cmd_3211'); load('r_cmd_3211');
load('V_x_cmd'); load('V_z_cmd');
load('t'); 
load('X_theta_x.mat'); load('Z_theta_z.mat');
coaxial_heli_parameters;

%% Simulation Parameters
tEnd = 14;
simN = 400;
u(1)=0;
v(1)=0;
w(1)=0;

% TRIM SPEED
[~, idx_trimspeed] = min(abs(V_vals - u(1)));


% U = zeros(6,simN);
U(1,1) = trim_saved_6dof(1,idx_trimspeed);  % Coll trim
U(2,1) = trim_saved_6dof(2,idx_trimspeed);  % DiffColl trim
U(3,1) = trim_saved_6dof(3,idx_trimspeed);  % LonCyc trim
U(4,1) = trim_saved_6dof(4,idx_trimspeed);  % LatCyc trim
U(5,1) = trim_saved_6dof(5,idx_trimspeed);  % roll trim
U(6,1) = trim_saved_6dof(6,idx_trimspeed);  % PropColl trim
U(7,1) = trim_saved_6dof(7,idx_trimspeed);  % quasidynamic inflow state trim upper rotor
U(8,1) = trim_saved_6dof(8,idx_trimspeed);  % quasidynamic inflow state trim lower rotor
U(9,1) = trim_saved_6dof(9,idx_trimspeed);  % quasidynamic inflow state trim pusher propeller
t(1)=0;
delta_e(1) = 0;
delta_r(1) = 0;

p(1)=0;
q(1)=0;
r(1)=0;
phi(1)=0;
theta(1)=0;
psi(1)=0;

% Time Paramaters
dt = (tEnd-t(1))/simN;

% Reset t
t = 0;

%% Simulation
for i = 1:simN

    
    % Speed definitions
    V(i) = sqrt(u(i)^2 + v(i)^2 + w(i)^2);

    % Select Derivatives for Correct Speed
    [~, idx_speed] = min(abs(V_vals - V(i)));

    % Keep certain inputs constant for t
    U(1,i) = 0;
    U(2,i) = 0;    % DiffColl
    U(3,i) = 0;
    U(4,i) = 0;
    U(5,i) = 0;    % Roll
    U(6,i) = 0;
    delta_e(i) = 0;
    delta_r(i) = 0;


    theta_cmd(i) = ref_3211(i);

    % EMF parameters
    M_theta = tf([4^2], [1, 2*0.707*4, 4^2]);
    M_theta_d = c2d(M_theta, dt);
    % Convert to state-space representation
    [num, den] = tfdata(M_theta_d, 'v');
    theta_c(i) = filter(num, den, theta_cmd(i));

    % State Update
    % [xdot(:,i)] = f_xk6(vel(:,i), U(:,i), p(i), q(i), r(i), theta(i), delta_e(i), delta_r(i), theta_cdiff(i));
    A_lin = TrimLinSave_6dof{idx_speed,1};
    B_lin = TrimLinSave_6dof{idx_speed,2};
    C_lin = TrimLinSave_6dof{idx_speed,3};
    D_lin = TrimLinSave_6dof{idx_speed,4};
    linSS = ss(A_lin, B_lin, C_lin, D_lin);
    linSSd = c2d(linSS, dt);
    A_lin_d = linSSd.A;
    B_lin_d = linSSd.B;
    C_lin_d = linSSd.C;
    D_lin_d = linSSd.D;
    x = [u(i); v(i); w(i); p(i); q(i); r(i)];
    % u_c(:,i) = [U(1,i); U(2,i); U(3,i); U(4,i); U(5,i); U(6,i); delta_e(i); delta_r(i)];
    u_c = 0;
    u_trim(:,i) = [trim_saved_6dof(1:end-3,idx_speed); 0; 0]; % last 2 rows are delta_e and delta_r equal to 0 in trim
    % u_trim = 0;
    u_real = u_c + u_trim;

    x = A_lin_d*x + B_lin_d*u_real(:,i);

    % Euler Integration
    u(i+1) = x(1);
    v(i+1) = x(2);
    w(i+1) = x(3);
    p(i+1) = x(4);
    q(i+1) = x(5);
    r(i+1) = x(6);

    phi(i+1) = phi(i)+dt*p(i);
    theta(i+1) = theta(i)+dt*q(i);
    psi(i+1) = psi(i)+dt*r(i);

    t(i+1) = t(i)+dt;

end

%% PLOTTING
figure(81)
subplot(3,2,1)
plot(t, u); grid on; legend('u');
subplot(3,2,3)
plot(t, v); grid on; legend('v');
subplot(3,2,5)
plot(t, w); grid on; legend('w');
subplot(3,2,2)
plot(t, rad2deg(phi)); grid on; legend('phi');
subplot(3,2,4)
plot(t, rad2deg(theta)); grid on; legend('theta');
subplot(3,2,6)
plot(t, rad2deg(psi)); grid on; legend('psi');



% INPUTS
figure(82)
subplot(4,2,1)
plot(t(1:end-1), rad2deg(u_real(1,:))); legend('Coll'); grid on;
subplot(4,2,3)
plot(t(1:end-1), rad2deg(u_real(2,:))); legend('DiffColl'); grid on;

subplot(4,2,5)
plot(t(1:end-1), rad2deg(u_real(3,:))); legend('LonCyc'); grid on;
subplot(4,2,7)
plot(t(1:end-1), rad2deg(u_real(4,:))); legend('LatCyc'); grid on;

subplot(4,2,2)
plot(t(1:end-1), rad2deg(u_real(6,:))); legend('PropColl');grid on;

subplot(4,2,4)
plot(t(1:end-1), rad2deg(u_real(7,:))); legend('Elev'); grid on;
subplot(4,2,6)
plot(t(1:end-1), rad2deg(u_real(8,:))); legend('Rudder'); grid on;


