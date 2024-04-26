%% Nonlinear Simulation

clear;clc;close all;

% Load Newest Trim States and Helicopter Parameters
load('trim_saved.mat');load('qpitchhq1.mat');load('qpitchhq2.mat');
load('TrimLinSave.mat'); load('V_transit_save.mat'); load('theta_ref');
load('t'); load('theta_cmd');
coaxial_heli_parameters;

%% Simulation Parameters

% TRIM SPEED
u(1) = 0;
v(1) = 0;
w(1) = 0;
[~, idx_trimspeed] = min(abs(V_vals - u(1)));

simN = 400;
% U = zeros(6,simN);
U(1,1) = deg2rad(12);  % collective at hover trim
U(2,1) = 0;
U(3,1) = deg2rad(2);  % longitudinal cyclic at hover trim
U(4,1) = deg2rad(0);
U(5,1) = 0;
U(6,1) = deg2rad(0);  % prop collective at hover trim
U(7,1) = 0;
U(8,1) = 0;
U(9,1) = 0;
t(1)=0;

p(1)=0;
q(1)=0;
r(1)=0;
phi(1)=0;
theta(1)=0;
psi(1) = 0;

% Time Paramaters
tEnd = 6;
dt = (tEnd-t(1))/simN;

% Reset t
t = 0;

%% Simulation
for i = 1:simN

    
    % Speed definitions
    V(i) = sqrt(u(i)^2 + v(i)^2 + w(i)^2);
    vel(:,i) = [u(i); v(i); w(i); 0; 0; 0];

    % Constant Rotor and Prop Collective
    U(1,i) = U(1,1);
    U(3,i) = U(3,1);
    U(4,i) = U(4,1);
    U(5,i) = U(5,1);
    U(6,i) = U(6,1);

    % % Select Longitudinal Control Derivatives for Correct Speed
    % [~, idx_speed] = min(abs(V_vals - V(i)));
    % M_theta_c(i) = TrimLinSave{idx_speed,2}(3,2);
    % M_delta_e(i) = TrimLinSave{idx_speed,2}(3,4);

    % FREE RESPONSE after 1 deg lon cyc input
    delta_e(i) = 0;
    delta_r(i) = 0;
    if t(i)<1
        U(2,i) = U(2,1);
    end
    if t(i)>1
        U(2,i) = deg2rad(3);
    end
    if t(i)>2 
        U(2,i) = U(2,1);
    end
    

    % State Update
    [xdot(:,i), M_components(:,i)] = f_xk6(vel(:,i), U(:,i), p(i), q(i), r(i), theta(i), delta_e(i), delta_r(i));

    % Euler Integration
    u(i+1) = u(i)+dt*xdot(1,i);
    v(i+1) = v(i)+dt*xdot(2,i);
    w(i+1) = w(i)+dt*xdot(3,i);
    p(i+1) = p(i)+dt*xdot(4,i);
    q(i+1) = q(i)+dt*xdot(5,i);
    r(i+1) = r(i)+dt*xdot(6,i);
    U(7,i+1) = U(7,i)+dt*xdot(7,i);
    U(8,i+1) = U(8,i)+dt*xdot(8,i);

    phi(i+1) = phi(i)+dt*p(i);
    theta(i+1) = theta(i)+dt*q(i);
    psi(i+1) = psi(i)+dt*r(i);
    
    t(i+1) = t(i)+dt;

end

%% PLOTTING
figure(20)
subplot(4,1,1)
plot(t(1:end-1), rad2deg(U(3,1:end-1))); 
grid on; legend('LonCyc [deg]')

subplot(4,1,2)
plot(t, u, t, v, t, w); grid on; legend('u [m/s]', 'v [m/s]', 'w [m/s]')

subplot(4,1,3)
plot(t, p, t, q, t, r); grid on; legend('p [rad/s]', 'q [rad/s]', 'r [rad/s]')

subplot(4,1,4)
plot(t, rad2deg(theta), t, rad2deg(phi), t, rad2deg(psi), 'LineWidth', 1.5); grid on;
legend('\phi_f [deg]', '\theta_f [deg]', '\psi_f [deg]')

figure(21)
subplot(3,1,1)
plot(t, rad2deg(U(2,:))); grid on; legend('Differential Collective [deg]')
subplot(3,1,2)
plot(t, r); grid on; legend('r [rad/s]')
subplot(3,1,3)
plot(t, rad2deg(psi)); grid on; legend('\psi [deg]')

