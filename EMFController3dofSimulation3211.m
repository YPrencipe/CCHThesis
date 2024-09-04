%% Nonlinear Simulation

clear;clc;close all;

% Load Newest Trim States and Helicopter Parameters
load('trim_saved.mat');load('qpitchhq1.mat');load('qpitchhq2.mat');
load('TrimLinSave.mat'); load('V_transit_save.mat'); 
load('ref_3211'); load('V_x_cmd'); load('V_z_cmd');
load('t'); 
load('X_theta_x.mat'); load('Z_theta_z.mat');
coaxial_heli_parameters;

%% Simulation Parameters

ManoeuvreSelection = 2; % 1. 3-2-1-1 pitch
                        % 2. Bob up/down w accel decell

if ManoeuvreSelection == 1
    tEnd = 14;
    simN = 400;
    u(1)=0;
    w(1)=0;
end
if ManoeuvreSelection == 2
    tEnd = 70;
    simN = 400;
    u(1)=40;
    w(1)=0;
    theta_cmd(1) = 0;
end

% TRIM SPEED
[~, idx_trimspeed] = min(abs(V_vals - u(1)));


% U = zeros(6,simN);
U(1,1) = trim_saved(1,idx_trimspeed);  % collective at hover trim
U(2,1) = trim_saved(2,idx_trimspeed);  % longitudinal cyclic at hover trim
U(3,1) = trim_saved(3,idx_trimspeed);  % prop collective at hover trim
U(4,1) = trim_saved(4,idx_trimspeed);
U(5,1) = trim_saved(5,idx_trimspeed);
U(6,1) = trim_saved(6,idx_trimspeed);
t(1)=0;

q(1)=0;
theta(1)=0;

% Time Paramaters
dt = (tEnd-t(1))/simN;

% CA ON or OFF
CA = 1;
V_l_transit = V_transit_save(1);
V_u_transit = V_transit_save(2)-20;

% Control Parameters Pitch
if CA == true
    if ManoeuvreSelection == 1
        Kp_theta = 150;
        Ki_theta = 1;
        Kp_q = -11; % THE GAINS ARE PROPORTIONATE TO THE B MATRIX FOR CA
        r_theta = 1.4;
    end

    if ManoeuvreSelection == 2
        Kp_theta = 45;
        Ki_theta = 0.1;
        Kp_q = -6; % THE GAINS ARE PROPORTIONATE TO THE B MATRIX FOR CA
        r_theta = 1.4;
    end

    % CA Outer OFF
    % Kp_x = -deg2rad(4);
    % Ki_x = -deg2rad(0.1);
    % CA Outer ON
    Kp_x = -0.7;
    Ki_x = -0.05;

    Kp_z = 1.2;
    Ki_z = 0.05;
end
if CA == false
    if ManoeuvreSelection == 1
        Kp_theta = -0.9;
        Ki_theta = -0.9;
        Kp_q = 0.6;
        r_theta = 1.4;
    end
    if ManoeuvreSelection == 2
        Kp_theta = -0.3;
        Ki_theta = -0.1;
        Kp_q = 0.07;
        r_theta = 1;
    end

    Kp_x = -deg2rad(4);
    Ki_x = -deg2rad(0.1);

    Kp_z = 0;
    Ki_z = 0;
end

% Initialize list to store Q values
Q_values = [];
r_values = [];

% Iterate over different values of r
% for r = 0.2:0.1:4

% EMF Control Scaling Parameter
r_x = 1;
r_z = 1;
integral_theta = 0;
integral_x = 0;
integral_z = 0;

% Reset t
t = 0;

%% Simulation
for i = 1:simN

    
    % Speed definitions
    V(i) = sqrt(u(i)^2+w(i)^2);
    vel(:,i) = [u(i); w(i); 0; 0; 0];

    % Constant Rotor and Prop Collective
    % U(1,i) = U(1,1);
    

    % Select Longitudinal Control Derivatives for Correct Speed
    [~, idx_speed] = min(abs(V_vals - V(i)));
    M_theta_c = TrimLinSave{idx_speed,2}(3,2);
    M_delta_e = TrimLinSave{idx_speed,2}(3,4);

    % X_theta_x = X_theta_x(idx_speed);
    X_theta_p = TrimLinSave{idx_speed,2}(1,3);

    Z_theta_0 = TrimLinSave{idx_speed,2}(2,1);
    % Z_theta_z = Z_theta_z(idx_speed);
    


    % Inner Loop ACAH Controllers
    error_theta(i) = theta_cmd(i) - theta(i);
    integral_theta = integral_theta + error_theta(i) * dt;
    gamma_theta(i) = r_theta * (Kp_theta*error_theta(i) + Ki_theta * integral_theta + Kp_q*q(i) );

    % Outer Loop Controllers
    error_x(i) = V_x_cmd(i) - u(i);
    integral_x = integral_x + error_x(i) * dt;
    gamma_x(i) = r_x * (Kp_x*error_x(i) + Ki_x*integral_x);

    error_z = V_z_cmd(i) - w(i);
    integral_z = integral_z + error_z * dt;
    gamma_z(i) = r_z * (Kp_z*error_z + Ki_z*integral_z);

    % If not using Control Allocation
    if CA == false
        % Pitch
        U(2,i) = gamma_theta(i);
        delta_e(i) = 0;

        % LonVel
        theta_x_cmd(i) = gamma_x(i);
        theta_p(i) = U(3,1);

        % Height
        theta_z_cmd(i) = 0;
        U(1,i) = U(1,1);

    end

    % Define Reference 3-2-1-1 Signal
    if ManoeuvreSelection == 1
        if t(i) < 1
            theta_ref(i) = 0;
        elseif t(i) >= 1 && t(i) < 4 
            theta_ref(i) = deg2rad(1);
        elseif t(i) >= 4 && t(i) < 6 
            theta_ref(i) = deg2rad(-1);
        elseif t(i) >= 6 && t(i) < 7 
            theta_ref(i) = deg2rad(1);
        elseif t(i) >= 7 && t(i) < 8 
            theta_ref(i) = deg2rad(-1);
        else t(i) >= 8;
            theta_ref(i) = 0;
        end
    end

    % Define Bob-Up Bob-Down with Accel/Decell Manoeuvre
    if t(i) < 5
        V_x_ref(i) = u(1);
    elseif t(i) >= 5 && t(i) < 25
        V_x_ref(i) = 30;
    elseif t(i) >= 25 && t(i) < 40
        V_x_ref(i) = 40;
    elseif t(i) >= 40
        V_x_ref(i) = 30;
    end

    if t(i) < 10
        V_z_ref(i) = w(1);
    elseif t(i) >= 10 && t(i) < 20
        V_z_ref(i) = -2;
    elseif t(i) >= 20 && t(i) < 45
        V_z_ref(i) = 0;
    elseif t(i) >= 45 && t(i) < 55
        V_z_ref(i) = 2;
    elseif t(i) >= 55
        V_z_ref(i) = 0;
    end

    % CA Inner Loop
    if V(i)<V_l_transit
        rho_theta_c(i) = 1;
        rho_delta_e(i) = epsilon;
    elseif V(i)>V_l_transit && V(i)<V_u_transit
        rho_theta_c(i) = 1 - 1/(V_u_transit-V_l_transit)*(V(i)-V_l_transit);
        rho_delta_e(i) = 1/(V_u_transit-V_l_transit)*(V(i)-V_l_transit);
    else
        rho_theta_c(i) = epsilon;
        rho_delta_e(i) = 1;
    end

    % CA Outer Loop
    if V(i)<V_l_transit
        rho_theta_x_cmd(i) = 1;
        rho_theta_p(i) = epsilon;
        rho_theta_0(i) = 1;
        rho_theta_z_cmd(i) = epsilon;
    elseif V(i)>V_l_transit && V(i)<V_u_transit
        rho_theta_x_cmd(i) = 1 - 1/(V_u_transit-V_l_transit)*(V(i)-V_l_transit);
        rho_theta_p(i) = 1/(V_u_transit-V_l_transit)*(V(i)-V_l_transit);
        rho_theta_0(i) = 1 - 1/(V_u_transit-V_l_transit)*(V(i)-V_l_transit);
        rho_theta_z_cmd(i) = 1/(V_u_transit-V_l_transit)*(V(i)-V_l_transit);
    else
        rho_theta_x_cmd(i) = epsilon;
        rho_theta_p(i) = 1;
        rho_theta_0(i) = epsilon;
        rho_theta_z_cmd(i) = 1;
    end


    if CA == true

        [~, idx_trimspeed] = min(abs(V_vals - u(1)));

        % LonVel CA
        % theta_x_cmd(i) = gamma_x(i);
        % U(3,i) = U(3,1);
        B_x = [1, -X_theta_p];
        % B_x = [-1, -1];
        W_x = [1/rho_theta_x_cmd(i), 0; 0, 1/rho_theta_p(i)];
        U_x(:,i) = inv(W_x)*transpose(B_x)* ...
            inv(B_x*inv(W_x)*transpose(B_x))*gamma_x(i);
        theta_x_cmd(i) = U_x(1,i);
        U(3,i) = U_x(2,i) + trim_saved(3,idx_trimspeed);

        % Height CA
        B_z = [Z_theta_0/2, 1];
        W_z = [1/rho_theta_0(i), 0; 0, 1/rho_theta_z_cmd(i)];
        U_z(:,i) = inv(W_z)*transpose(B_z)* ...
            inv(B_z*inv(W_z)*transpose(B_z))*gamma_z(i);
        U(1,i) = U_z(1,i) + trim_saved(1,idx_trimspeed);
        theta_z_cmd(i) = U_z(2,i);
        % U(1,i) = U(1,1);
        % theta_z_cmd(i) = 0;

        % Pitch CA
        B_theta = [M_theta_c, M_delta_e];
        W_theta = [1/rho_theta_c(i), 0; 0, 1/rho_delta_e(i)];
        U_theta(:,i) = inv(W_theta)*transpose(B_theta)* ...
            inv(B_theta*inv(W_theta)*transpose(B_theta))*gamma_theta(i);
        U(2,i) = U_theta(1,i) + trim_saved(2,idx_trimspeed);
        delta_e(i) = U_theta(2,i); 
    end

    

    % Actuator Rate Constraints
    if i>1
        % LonCyc
        cyclic_rate(i) = (U(2,i)-U(2,i-1))/dt;
        if abs(cyclic_rate(i)) > 28.8
            U(2,i) = sign(cyclic_rate(i))*28.8*dt + U(2,i-1);
        end
        cyclic_rate(i) = (U(2,i)-U(2,i-1))/dt;
        
        % Elev
        elev_rate(i) = (delta_e(i)-delta_e(i-1))/dt;
        if abs(elev_rate(i)) > 28.8
            delta_e(i) = sign(elev_rate(i))*28.8*dt + delta_e(i-1);
        end
        elev_rate(i) = (delta_e(i)-delta_e(i-1))/dt;

        % Coll
        coll_rate(i) = (U(1,i)-U(1,i-1))/dt;
        if abs(coll_rate(i)) > 16
            U(1,i) = sign(coll_rate(i))*16*dt + U(1,i-1);
        end
        coll_rate(i) = (U(1,i)-U(1,i-1))/dt;

        % PropColl
        propcoll_rate(i) = (U(3,i)-U(3,i-1))/dt;
        if abs(propcoll_rate(i)) > 16
            U(3,i) = sign(propcoll_rate(i))*16*dt + U(3,i-1);
        end
        propcoll_rate(i) = (U(3,i)-U(3,i-1))/dt;

    end

    % Actuator Constraints
        % LonCyc
        if abs(U(2,i)) > deg2rad(12)
            U(2,i) = sign(U(2,i))*deg2rad(12);
        end
        % Elev
        if abs(delta_e(i)) > deg2rad(20)
            delta_e(i) = sign(delta_e(i))*deg2rad(20);
        end
        % Coll
        if U(1,i) > deg2rad(20)
            U(1,i) = deg2rad(20);
        end
        if U(1,i) < deg2rad(0)
            U(1,i) = deg2rad(0);
        end
        % PropColl
        if U(3,i) > deg2rad(100)
            U(3,i) = deg2rad(100);
        end
        if U(3,i) < deg2rad(-10)
            U(3,i) = deg2rad(-10);
        end


    % Calc Total Required Pitch
    theta_cmd(i+1) = theta_x_cmd(i) + theta_z_cmd(i);

    % State Update
    [xdot(:,i)] = f_xk(vel(:,i), U(:,i), q(i), theta(i), delta_e(i));

    % Euler Integration
    u(i+1) = u(i)+dt*xdot(1,i);
    w(i+1) = w(i)+dt*xdot(2,i);
    q(i+1) = q(i)+dt*xdot(3,i);
    theta(i+1) = theta(i)+dt*q(i);
    U(4,i+1) = U(4,i)+dt*xdot(4,i);
    U(5,i+1) = U(5,i)+dt*xdot(5,i);
    U(6,i+1) = U(6,i)+dt*xdot(6,i);
    t(i+1) = t(i)+dt;

    

end

%%
% EMF Pre-defining
if ManoeuvreSelection == 1
    M_theta_tf = tf([4^2], [1, 2*0.707*4, 4^2]);
    theta_cmd = lsim(M_theta_tf, ref_3211, t(1:end-1));
    theta_cmd = theta_cmd';
    save('ref_3211', 'ref_3211');
    save('theta_cmd', 'theta_cmd');
end

if ManoeuvreSelection == 2
    M_V_z_tf = tf([3^2], [1, 2*0.7*3, 3^2]);
    V_z_cmd = lsim(M_V_z_tf, V_z_ref, t(1:end-1), V_z_ref(1));
    save('V_z_cmd', 'V_z_cmd');

    A = tf2ss([1^2], [1, 2*0.7*1, 1^2]);
    B = [0; -1];
    C = [1 0];
    D = 0;
    sys = ss(A,B,C,D);
    V_x_cmd = lsim(sys, V_x_ref, t(1:end-1), [40,-40*0.7/0.5]);
    V_x_cmd = V_x_cmd';
    save('V_x_cmd', 'V_x_cmd');
end

save('t', 't');


%%
% EMF Controller Performance (Q) Calculation
err = sum(abs(theta - theta_cmd));
V_theta = 1; %deg/in
err_rel = 1/(V_theta*tEnd)*err;
Q_theta = 1-err_rel

err_u = sum(abs(u(1:end-1) - V_x_cmd));
V_u = 1;
err_rel_u = 1/(V_u*tEnd)*err_u;
Q_u = 1-err_rel_u

% Save Q value in the list
Q_values = [Q_values, Q_theta];
r_values = [r_values, r_theta];

% end

disp(Q_values)

%% PLOTTING
% PLOT ACAH Response
% figure(71)
% subplot(3,1,1)
% plot(t(1:end-1), rad2deg(theta_ref), t(1:end-1), rad2deg(theta_cmd), ...
%      t(1:end-1), rad2deg(theta(1:end-1)));
% grid on; legend('3-2-1-1 ref', '\theta_{cmd} model', '\theta_f'); 
% ylabel('\theta [deg]'); xlabel('Time [s]');
% subplot(3,1,2)
% plot(t, rad2deg(U(2,:)), t(1:end-1), rad2deg(delta_e(:))); 
% grid on;
% legend('LonCyc', 'Elev');
% ylabel('Actuator Deflection [deg]')
% xlabel('Time [s]')
% subplot(3,1,3)
% plot(t,rad2deg(q)); grid on; legend('q [deg/s]');

% Model-Following Performance Criterion
% figure(72)
% scatter(r_values, Q_values*100, 100, 'Marker', 's', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'blue', 'LineWidth', 2);
% hold on;
% scatter(r_values, Q_values*100, 100, 'Marker', '.', 'MarkerEdgeColor', 'blue');
% plot(r_values, Q_values*100, '-b'); hold off;
% grid on;
% xlabel('r');
% ylabel('Q [%]');
% title('Optimization of Control System Activity');

% Plot Bob-Up/Down accell decell response
figure(81)
subplot(3,1,1)
plot(t(1:end-1), V_x_ref, t(1:end-1), V_x_cmd, t, u); grid on;
ylabel('V_x [m/s]')
legend('Reference Signal', 'Model Signal', 'Helicopter State')
subplot(3,1,2)
plot(t(1:end-1), V_z_ref, t(1:end-1), V_z_cmd, t, w); grid on;
ylabel('V_z [m/s]')
subplot(3,1,3)
plot(t, rad2deg(theta), t, rad2deg(theta_cmd)); grid on;
ylabel('\theta [deg]')
xlabel('Time [s]')
legend('\theta [deg]', '\theta_{cmd} [deg]')
% 
figure(82)
subplot(4,1,1)
plot(t, rad2deg(U(1,:))); legend('Coll'); grid on;
ylim([7.5,12.5])
subplot(4,1,2)
plot(t, rad2deg(U(2,:))); legend('LonCyc'); grid on;
subplot(4,1,3);
plot(t, rad2deg(U(3,:))); legend('PropColl');grid on;
subplot(4,1,4)
plot(t(1:end-1), rad2deg(rad2deg(delta_e))); legend('Elev'); grid on;

figure(85); 
plot(t, rad2deg(theta), t, rad2deg(theta_cmd), ...
    t(1:end-1), rad2deg(error_theta)); grid on;
ylabel('\theta [deg]')
xlabel('Time [s]')
legend('\theta [deg]', '\theta_{cmd} [deg]')



