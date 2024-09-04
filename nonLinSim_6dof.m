%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename:     nonLinSim_6dof
% Project:      MSc Thesis
% Supervisor:   M.D. Pavel
% Author:       Ynias Prencipe 
% Student Nr.:  4777158
% 
% Description:  Nonlinear simulation file used for control and manoeuvres
%               etc
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nonlinear Simulation

clear;clc; %close all;

% Load Newest Trim States and Helicopter Parameters
load('trim_saved.mat');load('qpitchhq1.mat');load('qpitchhq2.mat');
load('TrimLinSave.mat'); load('V_transit_save.mat'); load('theta_ref');
load('t'); load('theta_cmd'); load('trim_saved_6dof.mat')
coaxial_heli_parameters;

%% Simulation Parameters

% TRIM SPEED
u(1) = input('Initial speed [m/s]: ');
v(1) = 0;
w(1) = 0;
[~, idx_trimspeed] = min(abs(V_vals - u(1)));

simN = 400;
U(:,1) = trim_saved_6dof(:,idx_trimspeed);
t(1)=0;

p(1)=0;
q(1)=0;
r(1)=0;
phi(1)=0;
theta(1)=0;
psi(1) = 0;

% Time Paramaters
% tEnd = -0.19*u(1) + 24;
tEnd = 24;
dt = (tEnd-t(1))/simN;

% Reset t
t = 0;

% Disturbancechoice
distChoice = input('Choose Disturbance Choice \n1 = LonCyc, 2 = Elev : ');     % 1 is lonCyc 1 deg, 2 is elev 1 deg
if distChoice == 1
    nameInp = 'cyc';
else
    nameInp = 'elev';
end

%% Simulation
for i = 1:simN
    % Speed Definitions
    V(i) = sqrt(u(i)^2 + v(i)^2 + w(i)^2);
    vel(:,i) = [u(i); v(i); w(i); 0; 0; 0];

    % Constant Trim Variables from Initial Condition
    U(1,i) = U(1,1);
    U(3,i) = U(3,1);
    U(4,i) = U(4,1);
    U(5,i) = U(5,1);
    U(6,i) = U(6,1);

    % FREE RESPONSE after 1 deg lonCyc/Elev input
    delta_r(i) = 0;
    theta_cdiff(i) = 0;
    if distChoice == 1
        if t(i)<3
            U(3,i) = U(3,1);
            delta_e(i) = 0;
        end
        if t(i)>3
            U(3,i) = U(3,1) + deg2rad(1);
            delta_e(i) = 0;
        end
        if t(i)>4.5 
            U(3,i) = U(3,1);
            delta_e(i) = 0;
        end
    else
        if t(i)<3
            U(3,i) = U(3,1);
            delta_e(i) = 0;
        end
        if t(i)>3
            U(3,i) = U(3,1);
            delta_e(i) = deg2rad(1);
        end
        if t(i)>3.5 
            U(3,i) = U(3,1);
            delta_e(i) = 0;
        end
    end

    % State Update
    xdot(:,i) = f_xk6(vel(:,i), U(:,i), p(i), q(i), r(i), theta(i), delta_e(i), delta_r(i), theta_cdiff(i));

    % Euler Integration
    u(i+1) = u(i)+dt*xdot(1,i);
    v(i+1) = v(i)+dt*xdot(2,i);
    w(i+1) = w(i)+dt*xdot(3,i);
    p(i+1) = p(i)+dt*xdot(4,i);
    q(i+1) = q(i)+dt*xdot(5,i);
    r(i+1) = r(i)+dt*xdot(6,i);
    U(7,i+1) = U(7,i)+dt*xdot(7,i);
    U(8,i+1) = U(8,i)+dt*xdot(8,i);

    psidot(i) = (q(i)*sin(phi(i)) + r(i)*cos(phi(i))) / cos(theta(i));
    thetadot(i) = q(i)*cos(phi(i)) - r(i)*sin(phi(i));
    phidot(i) = p(i) + psidot(i)*sin(theta(i));

    phi(i+1) = phi(i)+dt*phidot(i);
    theta(i+1) = theta(i)+dt*thetadot(i);
    psi(i+1) = psi(i)+dt*psidot(i);
    
    t(i+1) = t(i)+dt;

end

resultsNonLinOpenLoop = [rad2deg(U(1,1:end-1));rad2deg(U(3,1:end-1));rad2deg(delta_e); ...
    u(1:end-1);v(1:end-1);-w(1:end-1);...
    p(1:end-1);q(1:end-1);r(1:end-1);...
    rad2deg(phi(1:end-1));rad2deg(theta(1:end-1));rad2deg(psi(1:end-1));...
    t(1:end-1)];

%%
nameFile = [num2str(u(1)) 'ms_' nameInp '.mat']

filePath = fullfile('nonLinSimSaves', nameFile);
save(filePath, 'resultsNonLinOpenLoop');

%% PLOTTING
f = figure(20);
t = resultsNonLinOpenLoop(13,:);

subplot(4,1,1)
plot(t, resultsNonLinOpenLoop(1,:),  t, resultsNonLinOpenLoop(2,:), ...
    t, resultsNonLinOpenLoop(3,:), 'LineWidth', 1.5); 
grid on; legend('Coll [deg]','LonCyc [deg]','Elev','Interpreter', 'latex', 'FontSize', 10)
ylim([-1;17])
xlim([0;12])

subplot(4,1,2)
plot(t, resultsNonLinOpenLoop(4,:), t, resultsNonLinOpenLoop(5,:), ...
    t, resultsNonLinOpenLoop(6,:),'LineWidth', 1.5); grid on; 
legend('u [m/s]', 'v [m/s]', 'w [m/s]', 'Interpreter', 'latex', 'FontSize', 10)
xlim([0;12])

subplot(4,1,3)
plot(t, rad2deg(resultsNonLinOpenLoop(7,:)), t, rad2deg(resultsNonLinOpenLoop(8,:)), ...
    t, rad2deg(resultsNonLinOpenLoop(9,:)),'LineWidth', 1.5); grid on; 
legend('p [deg/s]', 'q [deg/s]', 'r [deg/s]', 'Interpreter', 'latex', 'FontSize', 10)
xlim([0;12])

subplot(4,1,4)
plot(t, resultsNonLinOpenLoop(10,:), t, resultsNonLinOpenLoop(11,:), ...
    t,resultsNonLinOpenLoop(12,:), 'LineWidth', 1.5); grid on;
legend('$\phi_f$ [deg]', '$\theta_f$ [deg]', '$\psi_f$ [deg]', 'Interpreter', 'latex', 'FontSize', 10)
xlabel('Simulation runtime [s]', 'Interpreter', 'latex', 'FontSize', 13)
xlim([0;12])

if distChoice == 1
    distStr = '$\theta_{1s}$ = $\theta_{1s_{trim}}$ + 1 deg, $\delta_e$ = 0';
else
    distStr = '$\theta_{1s}$ = \theta_{1s_{trim}}, \delta_e = 1 deg';
end

if u(1) == 0
    sgtitle(['Nonlinear open-loop response at hover'], 'Interpreter', 'latex')
else
    sgtitle(['Nonlinear open-loop response at ' num2str(u(1)) ' m/s'], 'Interpreter', 'latex')
end
f.Position = [500 200 570 650];

