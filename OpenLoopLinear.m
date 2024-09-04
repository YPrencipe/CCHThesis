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
tEnd = 8;
simN = 800;
u(1)=30;
v(1)=0;
w(1)=0;

% TRIM SPEED
[~, idx_trimspeed] = min(abs(V_vals - u(1)));
U(:,1) = trim_saved_6dof(:,idx_trimspeed);

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

u_c(:,1) = [0;0;0;0;0;0;0;0];

%% Simulation
for i = 1:simN

    
    % Speed definitions
    V(i) = sqrt(u(i)^2 + v(i)^2 + w(i)^2);
    % Select Derivatives for Correct Speed
    [~, idx_speed] = min(abs(V_vals - V(i)));

    % State Update
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

    % u_c(:,i) = [U(1,i); U(2,i); U(3,i); U(4,i); U(5,i); U(6,i); delta_e(i); delta_r(i)];
    if t(i) < 3
        u_c(:,i) = 0;
    end
    if t(i) > 3
        u_c(:,i) = 0;
        u_c(3,i) = deg2rad(1);
    end
    if t(i) > 3.5
        u_c(:,i) = 0;
    end    

    u_trim(:,i) = [U(1,1); U(2,1);U(3,1);U(4,1);U(5,1);U(6,1); 0; 0]; % last 2 rows are delta_e and delta_r equal to 0 in trim
    u_real = u_c;% + u_trim;

    x = [u(i); v(i); w(i); p(i); q(i); r(i)];
    % x = A_lin_d*x + B_lin_d*u_real(:,i);
    x = A_lin*x + B_lin*u_real(:,i);

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


