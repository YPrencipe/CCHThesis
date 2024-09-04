%% Nonlinear Simulation

clear;clc;close all;

% Load Newest Trim States and Helicopter Parameters
load('trim_saved.mat');load('qpitchhq1.mat');load('qpitchhq2.mat');
load('TrimLinSave.mat'); load('V_transit_save.mat'); load('ref_3211');
load('t'); load('theta_cmd');
coaxial_heli_parameters;

%% Simulation Parameters

% TRIM SPEED
u(1)=0;
w(1)=0;
[~, idx_trimspeed] = min(abs(V_vals - u(1)));

simN = 400;
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
tEnd = 13;
dt = (tEnd-t(1))/simN;

% CA ON or OFF
CA = 1;
V_l_transit = V_transit_save(1);
V_u_transit = V_transit_save(2);

% Control Parameters Pitch
cmd_pitch = tf(4^2, [1 2*0.707*4 4^2]);
inv_pitch = tf(13, [0.26 1]);
if CA == true
    Kp = 50;
    Ki = 50;
    Kd = -33.3; % THE GAINS ARE PROPORTIONATE TO THE B MATRIX FOR CA
end
if CA == false
    Kp = -0.9;
    Ki = -0.9;
    Kd = 0.6;
end




% Reset t
t = 0;

%% Simulation
for i = 1:simN

    
    % Speed definitions
    V(i) = sqrt(u(i)^2+w(i)^2);
    vel(:,i) = [u(i); w(i); 0; 0; 0];

    % Constant Rotor and Prop Collective
    U(1,i) = U(1,1);
    U(3,i) = U(3,1);

    % Select Longitudinal Control Derivatives for Correct Speed
    [~, idx_speed] = min(abs(V_vals - V(i)));
    M_theta_c(i) = TrimLinSave{idx_speed,2}(3,2);
    M_delta_e(i) = TrimLinSave{idx_speed,2}(3,4);


    % LONGITUDINAL PID CONTROLLER
    if t(i) < 1
        theta_cmd(i) = deg2rad(0);

        U(2,i) = U(2,1);
        delta_e(i) = 0;
    else 
        theta_cmd(i) = deg2rad(15);

        gamma_theta(i) = Kp*(theta_cmd(i) - theta(i)) + ...
                      + Ki * trapz(t, (theta_cmd(i) - theta)) + ...
                      Kd*q(i);
        % If not using Control Allocation
        if CA == false
            U(2,i) = gamma_theta(i);
            delta_e(i) = 0;
        end
    end

    % LONGITUDINAL EMF CONTROL
    
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


    % Control Allocation
    epsilon = 0.001;
    if V(i)<V_l_transit
        rho_theta_s(i) = 1;
        rho_delta_e(i) = epsilon;
    elseif V(i)>V_l_transit && V(i)<V_u_transit
        rho_theta_s(i) = 1 - 1/(V_u_transit-V_l_transit)*(V(i)-V_l_transit);
        rho_delta_e(i) = 1/(V_u_transit-V_l_transit)*(V(i)-V_l_transit);
    else
        rho_theta_s(i) = epsilon;
        rho_delta_e(i) = 1;
    end

    if CA == true
        if t(i) > 1
            B_theta = [M_theta_c(i), M_delta_e(i)];
            W_theta = [1/rho_theta_s(i), 0; 0, 1/rho_delta_e(i)];
            U_theta(:,i) = inv(W_theta)*transpose(B_theta)* ...
                inv(B_theta*inv(W_theta)*transpose(B_theta))*gamma_theta(i);
            U(2,i) = U_theta(1,i);
            delta_e(i) = U_theta(2,i);
        end 
    end

    % Actuator Rate Constraints
    if i>1
        cyclic_rate(i) = (U(2,i)-U(2,i-1))/dt;
        if abs(cyclic_rate(i)) > 28.8
            U(2,i) = sign(cyclic_rate(i))*28.8*dt + U(2,i-1);
        end
        cyclic_rate(i) = (U(2,i)-U(2,i-1))/dt;
    end

    % Actuator Constraints
    if abs(U(2,i)) > deg2rad(12)
        U(2,i) = sign(U(2,i))*deg2rad(12);
    end


    % State Update
    [xdot(:,i), M_components(:,i)] = f_xk(vel(:,i), U(:,i), q(i), theta(i), delta_e(i));

    % Euler Integration
    u(i+1) = u(i)+dt*xdot(1,i);
    w(i+1) = w(i)+dt*xdot(2,i);
    q(i+1) = q(i)+dt*xdot(3,i);
    theta(i+1) = theta(i)+dt*q(i);
    U(4,i+1) = U(4,i)+dt*xdot(4,i);
    U(5,i+1) = U(5,i)+dt*xdot(5,i);
    t(i+1) = t(i)+dt;

end

% EMF Pre-defining
M_theta_tf = tf([4^2], [1, 2*0.707*4, 4^2]);
theta_cmd = lsim(M_theta_tf, theta_ref, t(1:end-1));

save('theta_ref', 'theta_ref');
save('theta_cmd', 'theta_cmd');
save('t', 't');

%% PLOTTING
figure(20)
subplot(3,1,1)
plot(t(1:end-1), rad2deg(U(2,1:end-1))); 
grid on; legend('LonCyc [deg]')
subplot(3,1,2)
plot(t, q); grid on; legend('q [rad/s]')
subplot(3,1,3)
plot(t, rad2deg(theta), 'LineWidth', 1.5); grid on;
hold on;
plot(t(1:end-1), rad2deg(theta_cmd), 'm:', 'LineWidth', 1.5); hold off;
legend('\theta_f [deg]')

figure(20)

% First subplot
subplot(3,1,1)
plot(t(1:end-1), rad2deg(U(2,1:end-1)), ...
    t(1:end-1), rad2deg(delta_e)); 
grid on; 
legend('LonCyc [deg]', 'Elev [deg]')

% Second subplot
subplot(3,1,2)
plot(t, rad2deg(q)); 
grid on; 


% Calculate and display maximum value of q
max_q_value = max(q);
q_pk_deg = rad2deg(max_q_value)
disp(['Maximum Value of q: ', num2str(max_q_value)]);

% Plot horizontal red dashed line for max value of q
hold on;
yline(rad2deg(max_q_value), '--r', 'Max Value of q');
hold off;
legend('q [rad/s]')

% Third subplot
subplot(3,1,3)
plot(t, rad2deg(theta), 'LineWidth', 1.5); 
grid on;
hold on;
plot(t(1:end-1), rad2deg(theta_cmd), 'm:', 'LineWidth', 1.5); 
hold off;


% Calculate rise time for theta
percent_10 = 0.1; % 10%
percent_90 = 0.9; % 90%
steady_state_value_theta = mean(rad2deg(theta(end-100:end))); % Assuming the steady state is the mean of last 100 samples
ten_percent_value_theta = percent_10 * steady_state_value_theta;
ninety_percent_value_theta = percent_90 * steady_state_value_theta;
idx_10_theta = find(rad2deg(theta) > ten_percent_value_theta, 1, 'first');
idx_90_theta = find(rad2deg(theta) > ninety_percent_value_theta, 1, 'first');
rise_time_theta = t(idx_90_theta) - t(idx_10_theta);

% Find maximum value of theta
max_value_theta = max(theta);
theta_pk_deg = rad2deg(max_value_theta)

% Find minimum value after the rise time for theta
min_value_after_rise_theta = min(theta(idx_90_theta:end));
rad2deg(min_value_after_rise_theta)

% Plot horizontal red dashed lines for theta
hold on;
yline(rad2deg(max_value_theta), '--r', 'Max Value');
yline(rad2deg(min_value_after_rise_theta), '--r', 'Min Value after Rise Time');
hold off;
legend('\theta_f [deg]')

% Display results for theta
disp(['Rise Time for \theta: ', num2str(rise_time_theta)]);
disp(['Maximum Value for \theta: ', num2str(max_value_theta)]);
disp(['Minimum Value after Rise Time for \theta: ', num2str(min_value_after_rise_theta)]);

% %% PLOT PITCH ATTITUDE QUICKNESS
Qpitch_5deg = max_q_value / (max_value_theta - theta(1))
figure(21)
plot(rad2deg(min_value_after_rise_theta), Qpitch_5deg, 'p', 'MarkerSize', 10)
grid on; hold on;
plot(qpitchhq1(:,1), qpitchhq1(:,2), 'k--', 'LineWidth', 2);
plot(qpitchhq2(:,1), qpitchhq2(:,2), 'k--', 'LineWidth', 2);
% Add text labels for the different levels
text(10, 1.8, 'Level 1', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
text(10, 0.8, 'Level 2', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
text(10, 0.2, 'Level 3', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
otherQvals = [  1.58, 1.55, 1.54, 1.52, 1.50, 1.48, ...
                1.55, 1.53, 1.52, 1.51, 1.50, 1.48, ...
                1.47, 1.46, 1.35, 1.35, 1.37, 1.39];
QvalsX = [  5 10 15 20 25 30 5 10 15 20 25 30 5 10 15 20 25 30];
plot(QvalsX(1:6), otherQvals(1:6), 'p', 'MarkerSize', 10)
plot(QvalsX(7:12), otherQvals(7:12), 'p', 'MarkerSize', 10)
plot(QvalsX(13:18), otherQvals(13:18), 'p', 'MarkerSize', 10)
hold off;
xlim([0, 30]); ylim([0, 3.5]);
xlabel('Minimum attitude change, \Delta\theta_{min} (deg)');
ylabel('$\frac{q_{pk}}{\Delta\theta_{pk}}$ (1/sec)', 'Interpreter', 'latex', 'Rotation', 0); 
title('Pitch Attitude Quickness');

figure(22)
plot(t(1:end-1), rad2deg(gamma_theta)); grid on; legend('\gamma_{\theta} [deg/s]')

%% EMF TRYOUTS
num_theta_tf = 4^2;
den_theta_tf = [1, 2*0.707*4, 4^2];
M_theta_tf = tf(num_theta_tf, den_theta_tf)
[M_theta_ss_A, M_theta_ss_B, M_theta_ss_C, M_theta_ss_D] = tf2ss(num_theta_tf, den_theta_tf);
M_theta_ss = ss(M_theta_ss_A, M_theta_ss_B, M_theta_ss_C, M_theta_ss_D);
step(M_theta_ss); hold on;
M_theta_ss_A(1,2) = -3;
M_theta_ss_C(1,2) = 3;
M_theta_ss = ss(M_theta_ss_A, M_theta_ss_B, M_theta_ss_C, M_theta_ss_D);
step(M_theta_ss); hold off;

% PLOT EMF Control Parameters
figure(71)
plot(t(1:end-1), rad2deg(theta_ref), t(1:end-1), rad2deg(theta_cmd));
grid on; legend('3-2-1-1 ref', '\theta_{cmd} model'); 
ylabel('\theta [deg]'); xlabel('Time [s]');






%%
% %% PLOT LONGITUDINAL PARAMETERS TIME SIMULATION
% 
% figure(21)
%     subplot(4,1,1)
%     plot(t(1:end-1), rad2deg(U(1,1:end-1)), t(1:end-1), rad2deg(U(2,1:end-1)), ...
%         t(1:end-1), rad2deg(delta_e)); 
%     grid on; legend('Coll [deg]', 'LonCyc [deg]', 'Elev [deg]')
% 
%     subplot(4,1,2)
%     plot(t, q); grid on; legend('q [rad/s]')
% 
%     subplot(4,1,3)
%     plot(t, rad2deg(theta)); grid on;
%     hold on;
%     plot(t(1:end-1), rad2deg(theta_cmd), 'g:', 'LineWidth', 1.5);
%     % Initialize variables
%     prev_sign = sign(theta(20));
%     xlocs = [];
%     % Iterate over data points to detect sign changes and plot vertical lines
%     for i = 21:numel(theta)
%         curr_sign = sign(theta(i));
%         if curr_sign ~= prev_sign
%             plot(t(i), rad2deg(theta(i)), 'r.', 'MarkerSize', 10); % Marker to indicate the sign change
%             line([t(i), t(i)], ylim, 'Color', 'r', 'LineStyle', '--'); % Vertical line
%             prev_sign = curr_sign; % Update previous sign
%             xlocs = [xlocs, t(i)]; % Store x-location for vertical line
%         end
%     end
%     legend('theta_f [deg]')
%     hold off;
% 
%     subplot(4, 1, 4);
%     plot(t, u, t, w);
%     grid on;
%     hold on;
%     for i = 1:numel(xlocs)
%         line([xlocs(i), xlocs(i)], ylim, 'Color', 'r', 'LineStyle', '--'); % Vertical line
%     end
%     legend('u [m/s]', 'w [m/s]')
%     hold off;

% %% PLOT TIME HISTORY OF MOMENT COMPONENTS
% 
% figure(22)
% plot(t(1:end-1), M_components(1,:), t(1:end-1), M_components(2,:), ...
%      t(1:end-1), M_components(3,:), t(1:end-1), M_components(4,:), ...
%      t(1:end-1), M_components(5,:), t(1:end-1), M_components(6,:));
% legend('M_{MR_u}', 'M_{MR_l}', 'M_{hinge}', 'M_{flap}', 'M_{fus}', 'M_{ht}')
% 
% 
% close all;














