%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% ADITTIONAL PLOTTING FOR REPORT %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

% Load Heli Parameters
coaxial_heli_parameters;

%% Empennage Stall Modelling
alpha_ht = -360:360;
Cl_ht = zeros(size(alpha_ht));
Cl_vt = zeros(size(alpha_ht));
range_indices = (alpha_ht >= -20 & alpha_ht <= 20);
Cl_ht(range_indices) = interp1([-20, 20], [-1.2, 1.2], alpha_ht(range_indices), 'linear');
% C_l_alpha_h is about 3.4
Cl_vt(range_indices) = interp1([-20, 20], [-1.4, 1.4], alpha_ht(range_indices), 'linear');
% C_l_beta_v is about 4


figure(1)
plot(alpha_ht(330:390), Cl_ht(330:390), alpha_ht(330:390), Cl_vt(330:390), ...
    'LineWidth', 2);
xlim([-30;30])
grid on; legend('Horizontal Stabiliser', 'Vertical Stabiliser', 'Interpreter', 'latex', ...
    'Location', 'northwest')
xlabel('Angle of attack, $\alpha_h$, $\beta_v$ [deg]', 'Interpreter', 'latex')
ylabel('Lift Coefficient, $C_l$ [-]', 'Interpreter', 'latex')


%% Trim Validation
% Parameters with RPM Scheduling OFF
load('trim_saved_6dof_noRPMscheduling.mat');
theta_0 = rad2deg(trim_saved_6dof(1,:)); 
theta_d = rad2deg(trim_saved_6dof(2,:)); 
theta_s = rad2deg(trim_saved_6dof(3,:));  
theta_c = rad2deg(trim_saved_6dof(4,:));  
phi     = rad2deg(trim_saved_6dof(5,:)); 
theta_p = rad2deg(trim_saved_6dof(6,:)); 

% Parameters with RPM Scheduling ON
load('trim_saved_6dof_scheduledRPM.mat');
theta_0_RPM = rad2deg(trim_saved_6dof(1,:)); 
theta_d_RPM = rad2deg(trim_saved_6dof(2,:)); 
theta_s_RPM = rad2deg(trim_saved_6dof(3,:));  
theta_c_RPM = rad2deg(trim_saved_6dof(4,:));  
phi_RPM     = rad2deg(trim_saved_6dof(5,:)); 
theta_p_RPM = rad2deg(trim_saved_6dof(6,:)); 

% Yuqing Qiu Results
load('QiuData.mat');
theta_0_Qiu = QiuData(1,:); 
theta_d_Qiu = QiuData(2,:);
theta_s_Qiu = QiuData(3,:);  
theta_c_Qiu = QiuData(4,:); 
theta_p_Qiu = QiuData(5,:); 

% Berger Results
load('BergData.mat');
theta_0_Berg = BergData(1,:); 
theta_s_Berg = BergData(2,:);  
theta_p_Berg = BergData(3,:); 



% Controls set to 0 for trim
delta_e = zeros(1, length(V_vals));
delta_r = zeros(1, length(V_vals));
theta_f = zeros(1, length(V_vals));



rowsPlot = 3;
colsPlot = 2;


figure(2)
subplot(rowsPlot, colsPlot, 1)
plot(V_vals(1:2:end), theta_0(1:2:end), '-o', ...
    V_vals(1:2:end), theta_0_RPM(1:2:end), '-*', ...
    V_vals(1:2:end), theta_0_Qiu, '--s', ...
    V_vals(1:2:end), theta_0_Berg, '--^', 'LineWidth', 1.2), grid;
ylabel("$\theta_0$ [$\deg$]", 'Interpreter', 'latex', 'FontSize', 14)
ylim([0,20])

subplot(rowsPlot, colsPlot, 2)
plot(V_vals(1:2:end), theta_d(1:2:end)/2, '-o', ...
    V_vals(1:2:end), theta_d_RPM(1:2:end)/2, '-*', ...
    V_vals(1:2:end), theta_d_Qiu, '--s', 'LineWidth', 1.2), grid;
ylim([-1,1])
ylabel("$\theta_d$ [$\deg$]", 'Interpreter', 'latex', 'FontSize', 14)

subplot(rowsPlot, colsPlot, 3)
plot(V_vals(1:2:end), theta_s(1:2:end), '-o', ...
    V_vals(1:2:end), theta_s_RPM(1:2:end), '-*', ...
    V_vals(1:2:end), theta_s_Qiu, '--s', ...
    V_vals(1:2:end), theta_s_Berg, '--^', 'LineWidth', 1.2), grid;
ylabel("$\theta_{1s}$ [$\deg$]", 'Interpreter', 'latex', 'FontSize', 14)

subplot(rowsPlot, colsPlot, 4)
plot(V_vals(1:2:end), theta_c(1:2:end), '-o', ... 
    V_vals(1:2:end), theta_c_RPM(1:2:end), '-*', ...
    V_vals(1:2:end), theta_c_Qiu, '--s', 'LineWidth', 1.2), grid;
ylim([-1,1])
ylabel("$\theta_{1c}$ [$\deg$]", 'Interpreter', 'latex', 'FontSize', 14)

subplot(rowsPlot, colsPlot, 5)
plot(V_vals(1:2:end), phi(1:2:end), '-o', V_vals(1:2:end), phi_RPM(1:2:end), '-*', 'LineWidth', 1.2), grid;
xlabel("Forward Velocity [m/s]", 'Interpreter', 'latex', 'FontSize', 14);
ylabel("$\phi$ [$\deg$]", 'Interpreter', 'latex', 'FontSize', 14)
ylim([-1,1])

subplot(rowsPlot, colsPlot, 6)
plot(V_vals(1:2:end), theta_p(1:2:end), '-o', ...
    V_vals(1:2:end), theta_p_RPM(1:2:end), '-*', ...
    V_vals(1:2:end), theta_p_Qiu, '--s', ...
    V_vals(1:2:end), theta_p_Berg, '--^', 'LineWidth', 1.2), grid;
xlabel("Forward Velocity [m/s]", 'Interpreter', 'latex', 'FontSize', 14);
ylabel("$\theta_p$ [$\deg$]", 'Interpreter', 'latex', 'FontSize', 14)

% subplot(rowsPlot, colsPlot, 7)
% plot(V_vals(1:2:end), rad2deg(delta_e(1:2:end)), '-o'), grid;
% 
% subplot(rowsPlot, colsPlot, 8)
% plot(V_vals(1:2:end), rad2deg(delta_r(1:2:end)), '-o'), grid;
% 
% subplot(rowsPlot, colsPlot, 9)
% plot(V_vals(1:2:end), rad2deg(theta_f(1:2:end)), '-o'), grid;
% xlabel("Forward Velocity [m/s]");
% 
% subplot(rowsPlot, colsPlot, 10)
% plot(V_vals(1:2:end), zeros(1, length(theta_f(1:2:end))), '-o'), grid;
% xlabel("Forward Velocity [m/s]"); 

legend('Model results', 'Model with RPM change', 'Qiu', 'Berger', ...
    'Location', 'northoutside', 'Orientation', 'horizontal', ...
    'Interpreter', 'latex', 'FontSize', 12);

%% Trim Validation For Paper
% Parameters with RPM Scheduling OFF
load('trim_saved_6dof_noRPMscheduling.mat');
theta_0 = rad2deg(trim_saved_6dof(1,:)); 
theta_d = rad2deg(trim_saved_6dof(2,:)); 
theta_s = rad2deg(trim_saved_6dof(3,:));  
theta_c = rad2deg(trim_saved_6dof(4,:));  
phi     = rad2deg(trim_saved_6dof(5,:)); 
theta_p = rad2deg(trim_saved_6dof(6,:)); 

% Parameters with RPM Scheduling ON
load('trim_saved_6dof_scheduledRPM.mat');
theta_0_RPM = rad2deg(trim_saved_6dof(1,:)); 
theta_d_RPM = rad2deg(trim_saved_6dof(2,:)); 
theta_s_RPM = rad2deg(trim_saved_6dof(3,:));  
theta_c_RPM = rad2deg(trim_saved_6dof(4,:));  
phi_RPM     = rad2deg(trim_saved_6dof(5,:)); 
theta_p_RPM = rad2deg(trim_saved_6dof(6,:)); 

% Yuqing Qiu Results
load('QiuData.mat');
theta_0_Qiu = QiuData(1,:); 
theta_d_Qiu = QiuData(2,:);
theta_s_Qiu = QiuData(3,:);  
theta_c_Qiu = QiuData(4,:); 
theta_p_Qiu = QiuData(5,:); 

% Berger Results
load('BergData.mat');
theta_0_Berg = BergData(1,:); 
theta_s_Berg = BergData(2,:);  
theta_p_Berg = BergData(3,:); 



% Controls set to 0 for trim
delta_e = zeros(1, length(V_vals));
delta_r = zeros(1, length(V_vals));
theta_f = zeros(1, length(V_vals));



rowsPlot = 5;
colsPlot = 1;


figure(2)
subplot(rowsPlot, colsPlot, 1)
plot(V_vals(1:2:end), theta_0(1:2:end), '-o', ...
    V_vals(1:2:end), theta_0_Qiu, '--s', ...
    V_vals(1:2:end), theta_0_Berg, '--^'), grid;
ylabel("$\theta_0$ [$\deg$]", 'Interpreter', 'latex', 'FontSize', 14)
ylim([0,20])

subplot(rowsPlot, colsPlot, 2)
plot(V_vals(1:2:end), theta_d(1:2:end)/2, '-o', ...
    V_vals(1:2:end), theta_d_Qiu, '--s'), grid;
ylim([-1,1])
ylabel("$\theta_d$ [$\deg$]", 'Interpreter', 'latex', 'FontSize', 14)

subplot(rowsPlot, colsPlot, 3)
plot(V_vals(1:2:end), theta_s(1:2:end), '-o', ...
    V_vals(1:2:end), theta_s_Qiu, '--s', ...
    V_vals(1:2:end), theta_s_Berg, '--^'), grid;
ylabel("$\theta_{1s}$ [$\deg$]", 'Interpreter', 'latex', 'FontSize', 14)

subplot(rowsPlot, colsPlot, 4)
plot(V_vals(1:2:end), theta_c(1:2:end), '-o', ... 
    V_vals(1:2:end), theta_c_Qiu, '--s'), grid;
ylim([-1,1])
ylabel("$\theta_{1c}$ [$\deg$]", 'Interpreter', 'latex', 'FontSize', 14)

subplot(rowsPlot, colsPlot, 5)
plot(V_vals(1:2:end), theta_p(1:2:end), '-o', ...
    V_vals(1:2:end), theta_p_Qiu, '--s', ...
    V_vals(1:2:end), theta_p_Berg, '--^'), grid;
xlabel("Forward Velocity [m/s]", 'Interpreter', 'latex', 'FontSize', 14);
ylabel("$\theta_p$ [$\deg$]", 'Interpreter', 'latex', 'FontSize', 14)

% subplot(rowsPlot, colsPlot, 7)
% plot(V_vals(1:2:end), rad2deg(delta_e(1:2:end)), '-o'), grid;
% 
% subplot(rowsPlot, colsPlot, 8)
% plot(V_vals(1:2:end), rad2deg(delta_r(1:2:end)), '-o'), grid;
% 
% subplot(rowsPlot, colsPlot, 9)
% plot(V_vals(1:2:end), rad2deg(theta_f(1:2:end)), '-o'), grid;
% xlabel("Forward Velocity [m/s]");
% 
% subplot(rowsPlot, colsPlot, 10)
% plot(V_vals(1:2:end), zeros(1, length(theta_f(1:2:end))), '-o'), grid;
% xlabel("Forward Velocity [m/s]"); 

legend('Model Results', 'Qiu', 'Berger', ...
    'Location', 'northoutside', 'Orientation', 'horizontal', ...
    'Interpreter', 'latex', 'FontSize', 12);


%% Trim Convergence
load('countTrimConvergence_fast.mat')
countFast = count;
load('countTrimConvergence_slow.mat')
countSlow = count;
figure(3)
bar([countFast(1:2:end)', countSlow(1:2:end)'], 'grouped'); grid on;
legend('$ x_{1,{i}} =  x_{final,(i-1)} $',  ...
    '$x_{1,{i}} = x_{1,1}$', 'Interpreter', 'latex', ...
    'FontSize', 14)
xlabel('Trim Point', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Number of iterations until threshold $\epsilon$', 'Interpreter', 'latex', 'FontSize', 14)

%% Flapping Angles for LonCyc Analysis
figure(3)

% Normal tiwst --> theta_tw = -10 deg
load('rotorAngles_normalTwist.mat')
a0_u = angles(1,:); a0_l = angles(2,:); a1_u = angles(3,:); 
a1_l = angles(4,:); a1r_u = angles(5,:); a1r_l = angles(6,:); 
alfa_sp = angles(7,:); alfa_cp = angles(8,:); alfa_dp_u = angles(9,:);
alfa_dp_l = angles(10,:);
plot(V_vals, a0_u, V_vals, a0_l,...
    V_vals, a1_u, V_vals, a1_l, ...
    V_vals, a1r_u, 'k--', V_vals, a1r_l, 'k--', ...
    V_vals, alfa_sp, 'g--o', V_vals, alfa_cp, 'b--o');
grid on;
legend('$a_{0_u}$', '$a_{0_l}$', '$a_{1_u}$', '$a_{1_l}$', '$a_{1r_u}$', '$a_{1r_l}$', ...
    '$\alpha_{sp}$', '$\alpha_{cp}$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Angle [deg]', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Forward velocity [m/s]', 'Interpreter', 'latex', 'FontSize', 14)

figure(4)
% Small tiwst --> theta_tw = -3 deg
load('rotorAngles_smallTwistMin3.mat')
a0_u = angles(1,:); a0_l = angles(2,:); a1_u = angles(3,:); 
a1_l = angles(4,:); a1r_u = angles(5,:); a1r_l = angles(6,:); 
alfa_sp = angles(7,:); alfa_cp = angles(8,:); alfa_dp_u = angles(9,:);
alfa_dp_l = angles(10,:);
plot(V_vals, a0_u, V_vals, a0_l,...
    V_vals, a1_u, V_vals, a1_l, ...
    V_vals, a1r_u, 'k--', V_vals, a1r_l, 'k--', ...
    V_vals, alfa_sp, 'g--o', V_vals, alfa_cp, 'b--o');
grid on;
legend('$a_{0_u}$', '$a_{0_l}$', '$a_{1_u}$', '$a_{1_l}$', '$a_{1r_u}$', '$a_{1r_l}$', ...
    '$\alpha_{sp}$', '$\alpha_{cp}$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Angle [deg]', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Forward velocity [m/s]', 'Interpreter', 'latex', 'FontSize', 14)


%% Stability Derivatives Plotting
% Mu,  Mw, Mq, Lv, Lp, Lq, Np

load('TrimLinSave_6dof.mat')
for i = 1:length(TrimLinSave_6dof)
    Xu(i) = TrimLinSave_6dof{i,1}(1,1);
    Xv(i) = TrimLinSave_6dof{i,1}(1,2);
    Xw(i) = TrimLinSave_6dof{i,1}(1,3);
    Xp(i) = TrimLinSave_6dof{i,1}(1,4);
    Xq(i) = TrimLinSave_6dof{i,1}(1,5);
    Xr(i) = TrimLinSave_6dof{i,1}(1,6);

    Yu(i) = TrimLinSave_6dof{i,1}(2,1);
    Yv(i) = TrimLinSave_6dof{i,1}(2,2);
    Yw(i) = TrimLinSave_6dof{i,1}(2,3);
    Yp(i) = TrimLinSave_6dof{i,1}(2,4);
    Yq(i) = TrimLinSave_6dof{i,1}(2,5);
    Yr(i) = TrimLinSave_6dof{i,1}(2,6);

    Zu(i) = TrimLinSave_6dof{i,1}(3,1);
    Zv(i) = TrimLinSave_6dof{i,1}(3,2);
    Zw(i) = TrimLinSave_6dof{i,1}(3,3);
    Zp(i) = TrimLinSave_6dof{i,1}(3,4);
    Zq(i) = TrimLinSave_6dof{i,1}(3,5);
    Zr(i) = TrimLinSave_6dof{i,1}(3,6);

    Lu(i) = TrimLinSave_6dof{i,1}(4,1);
    Lv(i) = TrimLinSave_6dof{i,1}(4,2);
    Lw(i) = TrimLinSave_6dof{i,1}(4,3);
    Lp(i) = TrimLinSave_6dof{i,1}(4,4);
    Lq(i) = TrimLinSave_6dof{i,1}(4,5);
    Lr(i) = TrimLinSave_6dof{i,1}(4,6);

    Mu(i) = TrimLinSave_6dof{i,1}(5,1);
    Mv(i) = TrimLinSave_6dof{i,1}(5,2);
    Mw(i) = TrimLinSave_6dof{i,1}(5,3);
    Mp(i) = TrimLinSave_6dof{i,1}(5,4);
    Mq(i) = TrimLinSave_6dof{i,1}(5,5);
    Mr(i) = TrimLinSave_6dof{i,1}(5,6);

    Nu(i) = TrimLinSave_6dof{i,1}(6,1);
    Nv(i) = TrimLinSave_6dof{i,1}(6,2);
    Nw(i) = TrimLinSave_6dof{i,1}(6,3);
    Np(i) = TrimLinSave_6dof{i,1}(6,4);
    Nq(i) = TrimLinSave_6dof{i,1}(6,5);
    Nr(i) = TrimLinSave_6dof{i,1}(6,6);
end

figure(5)
subplot(3,3,1)
plot(V_vals(1:2:end), Mu(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$M_u$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,2)
plot(V_vals(1:2:end), Mw(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$M_w$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,3)
plot(V_vals(1:2:end), Mq(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$M_q$', 'Interpreter', 'latex', 'FontSize', 14)

subplot(3,3,4)
plot(V_vals(1:2:end), Lv(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$L_v$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,5)
plot(V_vals(1:2:end), Lp(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$L_p$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,6)
plot(V_vals(1:2:end), Lq(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$L_q$', 'Interpreter', 'latex', 'FontSize', 14)

subplot(3,3,7)
plot(V_vals(1:2:end), Xu(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$X_u$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Airspeed [m/s]', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,8)
plot(V_vals(1:2:end), Zw(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Z_w$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Airspeed [m/s]', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,9)
plot(V_vals(1:2:end), Np(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$N_p$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Airspeed [m/s]', 'Interpreter', 'latex', 'FontSize', 14)


% Stability Derivatives - Longitudinal
figure(6)
subplot(3,3,1)
plot(V_vals(1:2:end), Xu(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$X_u$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,2)
plot(V_vals(1:2:end), Xw(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$X_w$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,3)
plot(V_vals(1:2:end), Xq(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$X_q$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,4)
plot(V_vals(1:2:end), Zu(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Z_u$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,5)
plot(V_vals(1:2:end), Zw(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Z_w$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,6)
plot(V_vals(1:2:end), Zq(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Z_q$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,7)
plot(V_vals(1:2:end), Mu(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$M_u$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Airspeed [m/s]', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,8)
plot(V_vals(1:2:end), Mw(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$M_w$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Airspeed [m/s]', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,9)
plot(V_vals(1:2:end), Mq(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$M_q$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Airspeed [m/s]', 'Interpreter', 'latex', 'FontSize', 14)

% Stability Derivatives - Lateral
figure(7)
subplot(3,3,1)
plot(V_vals(1:2:end), Yv(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Y_v$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,2)
plot(V_vals(1:2:end), Yp(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Y_p$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,3)
plot(V_vals(1:2:end), Yr(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Y_r$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,4)
plot(V_vals(1:2:end), Lv(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$L_v$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,5)
plot(V_vals(1:2:end), Lp(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$L_p$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,6)
plot(V_vals(1:2:end), Lr(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$L_r$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,7)
plot(V_vals(1:2:end), Nv(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$N_v$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Airspeed [m/s]', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,8)
plot(V_vals(1:2:end), Np(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$N_p$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Airspeed [m/s]', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,9)
plot(V_vals(1:2:end), Nr(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$N_r$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Airspeed [m/s]', 'Interpreter', 'latex', 'FontSize', 14)

% Stability Derivatives - Lateral into Longitudinal
figure(8)
subplot(3,3,1)
plot(V_vals(1:2:end), Xv(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$X_v$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,2)
plot(V_vals(1:2:end), Xp(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$X_p$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,3)
plot(V_vals(1:2:end), Xr(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$X_r$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,4)
plot(V_vals(1:2:end), Zv(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Z_v$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,5)
plot(V_vals(1:2:end), Zp(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Z_p$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,6)
plot(V_vals(1:2:end), Zr(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Z_r$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,7)
plot(V_vals(1:2:end), Mv(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$M_v$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Airspeed [m/s]', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,8)
plot(V_vals(1:2:end), Mp(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$M_p$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Airspeed [m/s]', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,9)
plot(V_vals(1:2:end), Mr(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$M_r$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Airspeed [m/s]', 'Interpreter', 'latex', 'FontSize', 14)

% Stability Derivatives - Longitudinal into Lateral
figure(9)
subplot(3,3,1)
plot(V_vals(1:2:end), Yu(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Y_u$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,2)
plot(V_vals(1:2:end), Yw(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Y_w$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,3)
plot(V_vals(1:2:end), Yq(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Y_q$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,4)
plot(V_vals(1:2:end), Lu(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$L_u$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,5)
plot(V_vals(1:2:end), Lw(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$L_w$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,6)
plot(V_vals(1:2:end), Lq(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$L_q$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,7)
plot(V_vals(1:2:end), Nu(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$N_u$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Airspeed [m/s]', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,8)
plot(V_vals(1:2:end), Nw(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$N_w$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Airspeed [m/s]', 'Interpreter', 'latex', 'FontSize', 14)
subplot(3,3,9)
plot(V_vals(1:2:end), Nq(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$N_q$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Airspeed [m/s]', 'Interpreter', 'latex', 'FontSize', 14)


%% Control Derivative Plots
coaxial_heli_parameters  
load('TrimLinSave_6dof.mat')
for i = 1:length(TrimLinSave_6dof)
    % collective
    X_theta_0(i) = TrimLinSave_6dof{i,2}(1,1);
    Y_theta_0(i) = TrimLinSave_6dof{i,2}(2,1);
    Z_theta_0(i) = TrimLinSave_6dof{i,2}(3,1);
    L_theta_0(i) = TrimLinSave_6dof{i,2}(4,1);
    M_theta_0(i) = TrimLinSave_6dof{i,2}(5,1);
    N_theta_0(i) = TrimLinSave_6dof{i,2}(6,1);

    % differential collective
    X_theta_d(i) = TrimLinSave_6dof{i,2}(1,2);
    Y_theta_d(i) = TrimLinSave_6dof{i,2}(2,2);
    Z_theta_d(i) = TrimLinSave_6dof{i,2}(3,2);
    L_theta_d(i) = TrimLinSave_6dof{i,2}(4,2);
    M_theta_d(i) = TrimLinSave_6dof{i,2}(5,2);
    N_theta_d(i) = TrimLinSave_6dof{i,2}(6,2);

    % LonCyc
    X_theta_s(i) = TrimLinSave_6dof{i,2}(1,3);
    Y_theta_s(i) = TrimLinSave_6dof{i,2}(2,3);
    Z_theta_s(i) = TrimLinSave_6dof{i,2}(3,3);
    L_theta_s(i) = TrimLinSave_6dof{i,2}(4,3);
    M_theta_s(i) = TrimLinSave_6dof{i,2}(5,3);
    N_theta_s(i) = TrimLinSave_6dof{i,2}(6,3);

    % LatCyc
    X_theta_c(i) = TrimLinSave_6dof{i,2}(1,4);
    Y_theta_c(i) = TrimLinSave_6dof{i,2}(2,4);
    Z_theta_c(i) = TrimLinSave_6dof{i,2}(3,4);
    L_theta_c(i) = TrimLinSave_6dof{i,2}(4,4);
    M_theta_c(i) = TrimLinSave_6dof{i,2}(5,4);
    N_theta_c(i) = TrimLinSave_6dof{i,2}(6,4);
    
    % differential LatCyc
    X_theta_cdiff(i) = TrimLinSave_6dof{i,2}(1,5);
    Y_theta_cdiff(i) = TrimLinSave_6dof{i,2}(2,5);
    Z_theta_cdiff(i) = TrimLinSave_6dof{i,2}(3,5);
    L_theta_cdiff(i) = TrimLinSave_6dof{i,2}(4,5);
    M_theta_cdiff(i) = TrimLinSave_6dof{i,2}(5,5);
    N_theta_cdiff(i) = TrimLinSave_6dof{i,2}(6,5);

    % PropColl
    X_theta_p(i) = TrimLinSave_6dof{i,2}(1,6);
    Y_theta_p(i) = TrimLinSave_6dof{i,2}(2,6);
    Z_theta_p(i) = TrimLinSave_6dof{i,2}(3,6);
    L_theta_p(i) = TrimLinSave_6dof{i,2}(4,6);
    M_theta_p(i) = TrimLinSave_6dof{i,2}(5,6);
    N_theta_p(i) = TrimLinSave_6dof{i,2}(6,6);

    % elev
    X_delta_e(i) = TrimLinSave_6dof{i,2}(1,7);
    Y_delta_e(i) = TrimLinSave_6dof{i,2}(2,7);
    Z_delta_e(i) = TrimLinSave_6dof{i,2}(3,7);
    L_delta_e(i) = TrimLinSave_6dof{i,2}(4,7);
    M_delta_e(i) = TrimLinSave_6dof{i,2}(5,7);
    N_delta_e(i) = TrimLinSave_6dof{i,2}(6,7);

    % elev
    X_delta_r(i) = TrimLinSave_6dof{i,2}(1,8);
    Y_delta_r(i) = TrimLinSave_6dof{i,2}(2,8);
    Z_delta_r(i) = TrimLinSave_6dof{i,2}(3,8);
    L_delta_r(i) = TrimLinSave_6dof{i,2}(4,8);
    M_delta_r(i) = TrimLinSave_6dof{i,2}(5,8);
    N_delta_r(i) = TrimLinSave_6dof{i,2}(6,8);


end

%collective
figure(10)
subplot(2,3,1)
plot(V_vals(1:2:end), X_theta_0(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$X_{\theta_0}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,2)
plot(V_vals(1:2:end), Y_theta_0(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Y_{\theta_0}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,3)
plot(V_vals(1:2:end), Z_theta_0(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Z_{\theta_0}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,4)
plot(V_vals(1:2:end), L_theta_0(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$L_{\theta_0}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,5)
plot(V_vals(1:2:end), M_theta_0(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$M_{\theta_0}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,6)
plot(V_vals(1:2:end), N_theta_0(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$N_{\theta_0}$', 'Interpreter', 'latex', 'FontSize', 14)

%differential collective
figure(11)
subplot(2,3,1)
plot(V_vals(1:2:end), X_theta_d(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$X_{\theta_d}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,2)
plot(V_vals(1:2:end), Y_theta_d(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Y_{\theta_d}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,3)
plot(V_vals(1:2:end), Z_theta_d(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Z_{\theta_d}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,4)
plot(V_vals(1:2:end), L_theta_d(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$L_{\theta_d}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,5)
plot(V_vals(1:2:end), M_theta_d(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$M_{\theta_d}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,6)
plot(V_vals(1:2:end), rad2deg(N_theta_d(1:2:end)), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$N_{\theta_d}$', 'Interpreter', 'latex', 'FontSize', 14)

%loncyc
figure(12)
subplot(2,3,1)
plot(V_vals(1:2:end), X_theta_s(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$X_{\theta_s}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,2)
plot(V_vals(1:2:end), Y_theta_s(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Y_{\theta_s}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,3)
plot(V_vals(1:2:end), Z_theta_s(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Z_{\theta_s}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,4)
plot(V_vals(1:2:end), L_theta_s(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$L_{\theta_s}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,5)
plot(V_vals(1:2:end), M_theta_s(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$M_{\theta_s}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,6)
plot(V_vals(1:2:end), N_theta_s(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$N_{\theta_s}$', 'Interpreter', 'latex', 'FontSize', 14)

%latcyc
figure(13)
subplot(2,3,1)
plot(V_vals(1:2:end), X_theta_c(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$X_{\theta_c}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,2)
plot(V_vals(1:2:end), Y_theta_c(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Y_{\theta_c}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,3)
plot(V_vals(1:2:end), Z_theta_c(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Z_{\theta_c}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,4)
plot(V_vals(1:2:end), L_theta_c(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$L_{\theta_c}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,5)
plot(V_vals(1:2:end), M_theta_c(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$M_{\theta_c}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,6)
plot(V_vals(1:2:end), N_theta_c(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$N_{\theta_c}$', 'Interpreter', 'latex', 'FontSize', 14)

%differential latcyc
figure(14)
subplot(2,3,1)
plot(V_vals(1:2:end), X_theta_cdiff(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$X_{\Delta \theta_c}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,2)
plot(V_vals(1:2:end), Y_theta_cdiff(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Y_{\Delta \theta_c}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,3)
plot(V_vals(1:2:end), Z_theta_cdiff(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Z_{\Delta \theta_c}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,4)
plot(V_vals(1:2:end), L_theta_cdiff(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$L_{\Delta \theta_c}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,5)
plot(V_vals(1:2:end), M_theta_cdiff(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$M_{\Delta \theta_c}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,6)
plot(V_vals(1:2:end), N_theta_cdiff(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$N_{\Delta \theta_c}$', 'Interpreter', 'latex', 'FontSize', 14)

%pusher propeller
figure(15)
subplot(2,3,1)
plot(V_vals(1:2:end), X_theta_p(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$X_{\theta_p}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,2)
plot(V_vals(1:2:end), Y_theta_p(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Y_{\theta_p}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,3)
plot(V_vals(1:2:end), Z_theta_p(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Z_{\theta_p}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,4)
plot(V_vals(1:2:end), L_theta_p(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$L_{\theta_p}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,5)
plot(V_vals(1:2:end), M_theta_p(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$M_{\theta_p}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,6)
plot(V_vals(1:2:end), N_theta_p(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$N_{\theta_p}$', 'Interpreter', 'latex', 'FontSize', 14)

%pusher propeller
figure(16)
subplot(2,3,1)
plot(V_vals(1:2:end), X_delta_e(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$X_{\delta_e}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,2)
plot(V_vals(1:2:end), Y_delta_e(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Y_{\delta_e}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,3)
plot(V_vals(1:2:end), Z_delta_e(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Z_{\delta_e}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,4)
plot(V_vals(1:2:end), L_delta_e(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$L_{\delta_e}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,5)
plot(V_vals(1:2:end), M_delta_e(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$M_{\delta_e}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,6)
plot(V_vals(1:2:end), N_delta_e(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$N_{\delta_e}$', 'Interpreter', 'latex', 'FontSize', 14)

%pusher propeller
figure(17)
subplot(2,3,1)
plot(V_vals(1:2:end), X_delta_r(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$X_{\delta_r}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,2)
plot(V_vals(1:2:end), Y_delta_r(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Y_{\delta_r}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,3)
plot(V_vals(1:2:end), Z_delta_r(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$Z_{\delta_r}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,4)
plot(V_vals(1:2:end), L_delta_r(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$L_{\delta_r}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,5)
plot(V_vals(1:2:end), M_delta_r(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$M_{\delta_r}$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(2,3,6)
plot(V_vals(1:2:end), N_delta_r(1:2:end), '-o', 'MarkerFaceColor', '#D95319', 'LineWidth', 2); grid on;
ylabel('$N_{\delta_r}$', 'Interpreter', 'latex', 'FontSize', 14)

%% Control Effectiveness
clc;clear; close all;
trim_6dof
close all;
%%

% Longitudinal
transition_min_val = 20;        % [m/s]
[~, idx_transstart] = min(abs(V - transition_min_val));
[~, idx_intersection] = min(abs(abs(M_theta_s) - abs(M_delta_e)));


f = figure(18);
subplot(2,1,1)
plot(V_vals, abs(M_theta_s), V_vals, abs(M_delta_e), '-*', 'LineWidth', 2);
grid on;
% xclabel('$V_f$ [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$CE_\theta$ ($1/s^2$)', 'Interpreter', 'latex', 'FontSize', 14);
hold on;

V_l_transit = V_vals(idx_transstart);
V_u_transit = V_vals(idx_intersection);

line([V_l_transit, V_l_transit], ylim, 'Color', 'r', 'LineStyle', '--');
line([V_u_transit, V_u_transit], ylim, 'Color', 'r', 'LineStyle', '--');
x_shade = [V_vals(idx_transstart), V_vals(idx_intersection)];
y_shade = ylim; 
fill([x_shade, fliplr(x_shade)], [y_shade(1), y_shade(1), y_shade(2), y_shade(2)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
text(mean(x_shade), 190, 'Transition region', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'FontSize', 11); % Adjust the font weight and size of the text
legend('$M_{\theta_s}$', '$M_{\delta_e}$', 'Location', 'southeast', 'FontSize', 14, 'Interpreter', 'latex');
hold off;
xlim([0;125])
title('Control effectiveness of $\theta_{1s}$ and $\delta_e$', 'Interpreter', 'latex', 'FontSize', 14)

% Directional
transition_min_val = 20;        % [m/s]
[~, idx_transstart] = min(abs(V - transition_min_val));
[~, idx_intersection] = min(abs(abs(rad2deg(N_theta_d)) - abs(N_delta_r)));

subplot(2,1,2)
plot(V_vals, rad2deg(N_theta_d), V_vals, abs(N_delta_r), '-*', 'LineWidth', 2);
grid on;
xlabel('$V$ [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$CE_\psi$ ($1/s^2$)', 'Interpreter', 'latex', 'FontSize', 14);
hold on;

V_l_transit = V_vals(idx_transstart);
V_u_transit = V_vals(idx_intersection);

line([V_l_transit, V_l_transit], ylim, 'Color', 'r', 'LineStyle', '--');
line([V_u_transit, V_u_transit], ylim, 'Color', 'r', 'LineStyle', '--');
x_shade = [V_vals(idx_transstart), V_vals(idx_intersection)];
y_shade = ylim; 
fill([x_shade, fliplr(x_shade)], [y_shade(1), y_shade(1), y_shade(2), y_shade(2)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
legend('$N_{\theta_d}$', '$N_{\delta_r}$', 'Location', 'southeast', 'FontSize', 14, 'Interpreter', 'latex');
hold off;
xlim([0;125])
title('Control effectiveness of $\theta_d$ and $\delta_r$', 'Interpreter', 'latex', 'FontSize', 14)
f.Position = [500 200 570 650];


%% Open-Loop Nonlinear Response
% Hover
f = figure(19);
cycData = load('./nonLinSimSaves/0ms_cyc.mat');
cycData = cycData.resultsNonLinOpenLoop;
elevData = load('./nonLinSimSaves/0ms_elev.mat');
elevData = elevData.resultsNonLinOpenLoop;

t = cycData(13,:);

subplot(4,1,1)
plot(t, cycData(1,:),  t, cycData(2,:), t, cycData(3,:), 'LineWidth', 1.5); hold on;
plot(t, elevData(1,:), '--', 'Color', '#0072BD');
plot(t, elevData(2,:), '--', 'Color', '#D95319');
plot(t, elevData(3,:), '--', 'Color', '#EDB120'); hold off;
grid on; legend('Coll [deg]','LonCyc [deg]','Elev [deg]','Interpreter', 'latex', 'FontSize', 10, 'Location', 'SouthEast')
ylim([-1;17])
xlim([0;24])

subplot(4,1,2)
plot(t, cycData(4,:), t, cycData(5,:), t, cycData(6,:),'LineWidth', 1.5); hold on;
plot(t, elevData(4,:), '--', 'Color', '#0072BD');
plot(t, elevData(5,:), '--', 'Color', '#D95319');
plot(t, elevData(6,:), '--', 'Color', '#EDB120'); hold off;
grid on; legend('u [m/s]', 'v [m/s]', 'w [m/s]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'SouthEast')
xlim([0;24])

subplot(4,1,3)
plot(t, rad2deg(cycData(7,:)), t, rad2deg(cycData(8,:)), t, rad2deg(cycData(9,:)),'LineWidth', 1.5); grid on;
hold on;
plot(t, rad2deg(elevData(7,:)), '--', 'Color', '#0072BD');
plot(t, rad2deg(elevData(8,:)), '--', 'Color', '#D95319');
plot(t, rad2deg(elevData(9,:)), '--', 'Color', '#EDB120'); hold off;
legend('p [deg/s]', 'q [deg/s]', 'r [deg/s]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'SouthEast')
xlim([0;24])

subplot(4,1,4)
plot(t, cycData(10,:), t, cycData(11,:), t,cycData(12,:), 'LineWidth', 1.5); grid on;
hold on;
plot(t, elevData(10,:), '--', 'Color', '#0072BD');
plot(t, elevData(11,:), '--', 'Color', '#D95319');
plot(t, elevData(12,:), '--', 'Color', '#EDB120'); hold off;
legend('$\phi_f$ [deg]', '$\theta_f$ [deg]', '$\psi_f$ [deg]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'NorthWest')
xlabel('Simulation runtime [s]', 'Interpreter', 'latex', 'FontSize', 13)
xlim([0;24])

sgtitle('Nonlinear open-loop response at hover', 'Interpreter', 'latex')
f.Position = [500 200 570 650];

% 40 m/s
clear;
f = figure(20);
cycData = load('./nonLinSimSaves/40ms_cyc.mat');
cycData = cycData.resultsNonLinOpenLoop;
elevData = load('./nonLinSimSaves/40ms_elev.mat');
elevData = elevData.resultsNonLinOpenLoop;

t = cycData(13,:);

subplot(4,1,1)
plot(t, cycData(1,:),  t, cycData(2,:), t, cycData(3,:), 'LineWidth', 1.5); hold on;
plot(t, elevData(1,:), '--', 'Color', '#0072BD');
plot(t, elevData(2,:), '--', 'Color', '#D95319');
plot(t, elevData(3,:), '--', 'Color', '#EDB120'); hold off;
grid on; legend('Coll [deg]','LonCyc [deg]','Elev [deg]','Interpreter', 'latex', 'FontSize', 10, 'Location', 'NorthWest')
ylim([-1;17])
xlim([0;24])

subplot(4,1,2)
plot(t, cycData(4,:), t, cycData(5,:), t, cycData(6,:),'LineWidth', 1.5); hold on;
plot(t, elevData(4,:), '--', 'Color', '#0072BD');
plot(t, elevData(5,:), '--', 'Color', '#D95319');
plot(t, elevData(6,:), '--', 'Color', '#EDB120'); hold off;
grid on; legend('u [m/s]', 'v [m/s]', 'w [m/s]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'NorthWest')
xlim([0;24])

subplot(4,1,3)
plot(t, rad2deg(cycData(7,:)), t, rad2deg(cycData(8,:)), t, rad2deg(cycData(9,:)),'LineWidth', 1.5); grid on;
hold on;
plot(t, rad2deg(elevData(7,:)), '--', 'Color', '#0072BD');
plot(t, rad2deg(elevData(8,:)), '--', 'Color', '#D95319');
plot(t, rad2deg(elevData(9,:)), '--', 'Color', '#EDB120'); hold off;
legend('p [deg/s]', 'q [deg/s]', 'r [deg/s]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'NorthWest')
xlim([0;24])

subplot(4,1,4)
plot(t, cycData(10,:), t, cycData(11,:), t,cycData(12,:), 'LineWidth', 1.5); grid on;
hold on;
plot(t, elevData(10,:), '--', 'Color', '#0072BD');
plot(t, elevData(11,:), '--', 'Color', '#D95319');
plot(t, elevData(12,:), '--', 'Color', '#EDB120'); hold off;
legend('$\phi_f$ [deg]', '$\theta_f$ [deg]', '$\psi_f$ [deg]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'NorthWest')
xlabel('Simulation runtime [s]', 'Interpreter', 'latex', 'FontSize', 13)
xlim([0;24])

sgtitle('Nonlinear open-loop response at 40 m/s', 'Interpreter', 'latex')
f.Position = [500 200 570 650];

% 80 m/s
clear;
f = figure(21);
cycData = load('./nonLinSimSaves/65ms_cyc.mat');
cycData = cycData.resultsNonLinOpenLoop;
elevData = load('./nonLinSimSaves/65ms_elev.mat');
elevData = elevData.resultsNonLinOpenLoop;

t = cycData(13,:);

subplot(4,1,1)
plot(t, cycData(1,:),  t, cycData(2,:), t, cycData(3,:), 'LineWidth', 1.5); hold on;
plot(t, elevData(1,:), '--', 'Color', '#0072BD');
plot(t, elevData(2,:), '--', 'Color', '#D95319');
plot(t, elevData(3,:), '--', 'Color', '#EDB120'); hold off;
grid on; legend('Coll [deg]','LonCyc [deg]','Elev [deg]','Interpreter', 'latex', 'FontSize', 10, 'Location', 'West')
ylim([-1;17])
xlim([0;12])

subplot(4,1,2)
plot(t, cycData(4,:), t, cycData(5,:), t, cycData(6,:),'LineWidth', 1.5); hold on;
plot(t, elevData(4,:), '--', 'Color', '#0072BD');
plot(t, elevData(5,:), '--', 'Color', '#D95319');
plot(t, elevData(6,:), '--', 'Color', '#EDB120'); hold off;
grid on; legend('u [m/s]', 'v [m/s]', 'w [m/s]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'West')
xlim([0;12])

subplot(4,1,3)
plot(t, rad2deg(cycData(7,:)), t, rad2deg(cycData(8,:)), t, rad2deg(cycData(9,:)),'LineWidth', 1.5); grid on;
hold on;
plot(t, rad2deg(elevData(7,:)), '--', 'Color', '#0072BD');
plot(t, rad2deg(elevData(8,:)), '--', 'Color', '#D95319');
plot(t, rad2deg(elevData(9,:)), '--', 'Color', '#EDB120'); hold off;
legend('p [deg/s]', 'q [deg/s]', 'r [deg/s]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'West')
xlim([0;12])

subplot(4,1,4)
plot(t, cycData(10,:), t, cycData(11,:), t,cycData(12,:), 'LineWidth', 1.5); grid on;
hold on;
plot(t, elevData(10,:), '--', 'Color', '#0072BD');
plot(t, elevData(11,:), '--', 'Color', '#D95319');
plot(t, elevData(12,:), '--', 'Color', '#EDB120'); hold off;
legend('$\phi_f$ [deg]', '$\theta_f$ [deg]', '$\psi_f$ [deg]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'West')
xlabel('Simulation runtime [s]', 'Interpreter', 'latex', 'FontSize', 13)
xlim([0;12])


sgtitle('Nonlinear open-loop response at 65 m/s', 'Interpreter', 'latex')
f.Position = [500 200 570 650];

% 100 m/s
clear;
f = figure(22);
cycData = load('./nonLinSimSaves/100ms_cyc.mat');
cycData = cycData.resultsNonLinOpenLoop;
elevData = load('./nonLinSimSaves/100ms_elev.mat');
elevData = elevData.resultsNonLinOpenLoop;

t = cycData(13,:);

subplot(4,1,1)
plot(t, cycData(1,:),  t, cycData(2,:), t, cycData(3,:), 'LineWidth', 1.5); hold on;
plot(t, elevData(1,:), '--', 'Color', '#0072BD');
plot(t, elevData(2,:), '--', 'Color', '#D95319');
plot(t, elevData(3,:), '--', 'Color', '#EDB120'); hold off;
grid on; legend('Coll [deg]','LonCyc [deg]','Elev [deg]','Interpreter', 'latex', 'FontSize', 10, 'Location', 'West')
ylim([-1;17])
xlim([0;4.5])

subplot(4,1,2)
plot(t, cycData(4,:), t, cycData(5,:), t, cycData(6,:),'LineWidth', 1.5); hold on;
plot(t, elevData(4,:), '--', 'Color', '#0072BD');
plot(t, elevData(5,:), '--', 'Color', '#D95319');
plot(t, elevData(6,:), '--', 'Color', '#EDB120'); hold off;
grid on; legend('u [m/s]', 'v [m/s]', 'w [m/s]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'West')
xlim([0;4.5])

subplot(4,1,3)
plot(t, rad2deg(cycData(7,:)), t, rad2deg(cycData(8,:)), t, rad2deg(cycData(9,:)),'LineWidth', 1.5); grid on;
hold on;
plot(t, rad2deg(elevData(7,:)), '--', 'Color', '#0072BD');
plot(t, rad2deg(elevData(8,:)), '--', 'Color', '#D95319');
plot(t, rad2deg(elevData(9,:)), '--', 'Color', '#EDB120'); hold off;
legend('p [deg/s]', 'q [deg/s]', 'r [deg/s]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'West')
xlim([0;4.5])

subplot(4,1,4)
plot(t, cycData(10,:), t, cycData(11,:), t,cycData(12,:), 'LineWidth', 1.5); grid on;
hold on;
plot(t, elevData(10,:), '--', 'Color', '#0072BD');
plot(t, elevData(11,:), '--', 'Color', '#D95319');
plot(t, elevData(12,:), '--', 'Color', '#EDB120'); hold off;
legend('$\phi_f$ [deg]', '$\theta_f$ [deg]', '$\psi_f$ [deg]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'West')
xlabel('Simulation runtime [s]', 'Interpreter', 'latex', 'FontSize', 13)
xlim([0;4.5])

sgtitle('Nonlinear open-loop response at 100 m/s', 'Interpreter', 'latex')
f.Position = [500 200 570 650];

%% Linear Open-Loop Response
clear;
f = figure(23);
cycData = load('./LinSimSaves/hov_cyc.mat');
cycData = cycData.out;
elevData = load('./LinSimSaves/hov_elev.mat');
elevData = elevData.out;

tCyc = cycData.tout;
tElev = elevData.tout;

subplot(4,1,1)
plot(tCyc, rad2deg(cycData.theta_0(:,2)),  tCyc, rad2deg(cycData.theta_1s(:,2)), tCyc, rad2deg(cycData.delta_e(:,2)), 'LineWidth', 1.5); hold on;
plot(tElev, rad2deg(elevData.theta_0(:,2)), '--', 'Color', '#0072BD');
plot(tElev, rad2deg(elevData.theta_1s(:,2)), '--', 'Color', '#D95319');
plot(tElev, rad2deg(elevData.delta_e(:,2)), '--', 'Color', '#EDB120'); hold off;
grid on; legend('Coll [deg]','LonCyc [deg]','Elev [deg]','Interpreter', 'latex', 'FontSize', 10, 'Location', 'NorthWest')
ylim([-1;17])
xlim([0;tCyc(end)+0.02])


subplot(4,1,2)
plot(tCyc, cycData.u(:,2), tCyc, cycData.v(:,2), tCyc, cycData.w(:,2),'LineWidth', 1.5); hold on;
plot(tElev, elevData.u(:,2), '--', 'Color', '#0072BD');
plot(tElev, elevData.v(:,2), '--', 'Color', '#D95319');
plot(tElev, elevData.w(:,2), '--', 'Color', '#EDB120'); hold off;
grid on; legend('u [m/s]', 'v [m/s]', 'w [m/s]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'SouthWest')
xlim([0;tCyc(end)])


subplot(4,1,3)
plot(tCyc, cycData.p(:,2), tCyc, cycData.q(:,2), tCyc, cycData.r(:,2),'LineWidth', 1.5); grid on;
hold on;
plot(tElev, elevData.p(:,2), '--', 'Color', '#0072BD');
plot(tElev, elevData.q(:,2), '--', 'Color', '#D95319');
plot(tElev, elevData.r(:,2), '--', 'Color', '#EDB120'); hold off;
legend('p [deg/s]', 'q [deg/s]', 'r [deg/s]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'SouthWest')
xlim([0;tCyc(end)])


subplot(4,1,4)
plot(tCyc, cycData.phi(:,2), tCyc, cycData.theta(:,2), tCyc,cycData.psi(:,2), 'LineWidth', 1.5); grid on;
hold on;
plot(tElev, elevData.phi(:,2), '--', 'Color', '#0072BD');
plot(tElev, elevData.theta(:,2), '--', 'Color', '#D95319');
plot(tElev, elevData.psi(:,2), '--', 'Color', '#EDB120'); hold off;
legend('$\phi_f$ [deg]', '$\theta_f$ [deg]', '$\psi_f$ [deg]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'SouthWest')
xlabel('Simulation runtime [s]', 'Interpreter', 'latex', 'FontSize', 13)
xlim([0;tCyc(end)])

sgtitle('Linear open-loop response at hover', 'Interpreter', 'latex')
f.Position = [500 200 570 650];


% 40 m/s
clear;
f = figure(24);
cycData = load('./LinSimSaves/40ms_cyc.mat');
cycData = cycData.out;
elevData = load('./LinSimSaves/40ms_elev.mat');
elevData = elevData.out;

tCyc = cycData.tout;
tElev = elevData.tout;

subplot(4,1,1)
plot(tCyc, rad2deg(cycData.theta_0(:,2)),  tCyc, rad2deg(cycData.theta_1s(:,2)), tCyc, rad2deg(cycData.delta_e(:,2)), 'LineWidth', 1.5); hold on;
plot(tElev, rad2deg(elevData.theta_0(:,2)), '--', 'Color', '#0072BD');
plot(tElev, rad2deg(elevData.theta_1s(:,2)), '--', 'Color', '#D95319');
plot(tElev, rad2deg(elevData.delta_e(:,2)), '--', 'Color', '#EDB120'); hold off;
grid on; legend('Coll [deg]','LonCyc [deg]','Elev [deg]','Interpreter', 'latex', 'FontSize', 10, 'Location', 'NorthWest')
ylim([-1;17])
xlim([0;tCyc(end)+0.02])

subplot(4,1,2)
plot(tCyc, cycData.u(:,2), tCyc, cycData.v(:,2), tCyc, cycData.w(:,2),'LineWidth', 1.5); hold on;
plot(tElev, elevData.u(:,2), '--', 'Color', '#0072BD');
plot(tElev, elevData.v(:,2), '--', 'Color', '#D95319');
plot(tElev, elevData.w(:,2), '--', 'Color', '#EDB120'); hold off;
grid on; legend('u [m/s]', 'v [m/s]', 'w [m/s]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'SouthWest')
xlim([0;tCyc(end)])

subplot(4,1,3)
plot(tCyc, cycData.p(:,2), tCyc, cycData.q(:,2), tCyc, cycData.r(:,2),'LineWidth', 1.5); grid on;
hold on;
plot(tElev, elevData.p(:,2), '--', 'Color', '#0072BD');
plot(tElev, elevData.q(:,2), '--', 'Color', '#D95319');
plot(tElev, elevData.r(:,2), '--', 'Color', '#EDB120'); hold off;
legend('p [deg/s]', 'q [deg/s]', 'r [deg/s]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'SouthWest')
xlim([0;tCyc(end)])

subplot(4,1,4)
plot(tCyc, cycData.phi(:,2), tCyc, cycData.theta(:,2), tCyc,cycData.psi(:,2), 'LineWidth', 1.5); grid on;
hold on;
plot(tElev, elevData.phi(:,2), '--', 'Color', '#0072BD');
plot(tElev, elevData.theta(:,2), '--', 'Color', '#D95319');
plot(tElev, elevData.psi(:,2), '--', 'Color', '#EDB120'); hold off;
legend('$\phi_f$ [deg]', '$\theta_f$ [deg]', '$\psi_f$ [deg]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'SouthWest')
xlabel('Simulation runtime [s]', 'Interpreter', 'latex', 'FontSize', 13)
xlim([0;tCyc(end)])

sgtitle('Linear open-loop response at 40 m/s', 'Interpreter', 'latex')
f.Position = [500 200 570 650];


% 80 m/s
clear;
f = figure(25);
cycData = load('./LinSimSaves/65ms_cyc.mat');
cycData = cycData.out;
elevData = load('./LinSimSaves/65ms_elev.mat');
elevData = elevData.out;

tCyc = cycData.tout;
tElev = elevData.tout;

subplot(4,1,1)
plot(tCyc, rad2deg(cycData.theta_0(:,2)),  tCyc, rad2deg(cycData.theta_1s(:,2)), tCyc, rad2deg(cycData.delta_e(:,2)), 'LineWidth', 1.5); hold on;
plot(tElev, rad2deg(elevData.theta_0(:,2)), '--', 'Color', '#0072BD');
plot(tElev, rad2deg(elevData.theta_1s(:,2)), '--', 'Color', '#D95319');
plot(tElev, rad2deg(elevData.delta_e(:,2)), '--', 'Color', '#EDB120'); hold off;
grid on; legend('Coll [deg]','LonCyc [deg]','Elev [deg]','Interpreter', 'latex', 'FontSize', 10, 'Location', 'NorthWest')
ylim([-1;17])
xlim([0;11])

subplot(4,1,2)
plot(tCyc, cycData.u(:,2), tCyc, cycData.v(:,2), tCyc, cycData.w(:,2),'LineWidth', 1.5); hold on;
plot(tElev, elevData.u(:,2), '--', 'Color', '#0072BD');
plot(tElev, elevData.v(:,2), '--', 'Color', '#D95319');
plot(tElev, elevData.w(:,2), '--', 'Color', '#EDB120'); hold off;
grid on; legend('u [m/s]', 'v [m/s]', 'w [m/s]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'SouthWest')
xlim([0;11])

subplot(4,1,3)
plot(tCyc, cycData.p(:,2), tCyc, cycData.q(:,2), tCyc, cycData.r(:,2),'LineWidth', 1.5); grid on;
hold on;
plot(tElev, elevData.p(:,2), '--', 'Color', '#0072BD');
plot(tElev, elevData.q(:,2), '--', 'Color', '#D95319');
plot(tElev, elevData.r(:,2), '--', 'Color', '#EDB120'); hold off;
legend('p [deg/s]', 'q [deg/s]', 'r [deg/s]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'SouthWest')
xlim([0;11])

subplot(4,1,4)
plot(tCyc, cycData.phi(:,2), tCyc, cycData.theta(:,2), tCyc,cycData.psi(:,2), 'LineWidth', 1.5); grid on;
hold on;
plot(tElev, elevData.phi(:,2), '--', 'Color', '#0072BD');
plot(tElev, elevData.theta(:,2), '--', 'Color', '#D95319');
plot(tElev, elevData.psi(:,2), '--', 'Color', '#EDB120'); hold off;
legend('$\phi_f$ [deg]', '$\theta_f$ [deg]', '$\psi_f$ [deg]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'SouthWest')
xlabel('Simulation runtime [s]', 'Interpreter', 'latex', 'FontSize', 13)
xlim([0;11])

sgtitle('Linear open-loop response at 65 m/s', 'Interpreter', 'latex')
f.Position = [500 200 570 650];


% 100 m/s
clear;
f = figure(26);
cycData = load('./LinSimSaves/100ms_cyc.mat');
cycData = cycData.out;
elevData = load('./LinSimSaves/100ms_elev.mat');
elevData = elevData.out;

tCyc = cycData.tout;
tElev = elevData.tout;

subplot(4,1,1)
plot(tCyc, rad2deg(cycData.theta_0(:,2)),  tCyc, rad2deg(cycData.theta_1s(:,2)), tCyc, rad2deg(cycData.delta_e(:,2)), 'LineWidth', 1.5); hold on;
plot(tElev, rad2deg(elevData.theta_0(:,2)), '--', 'Color', '#0072BD');
plot(tElev, rad2deg(elevData.theta_1s(:,2)), '--', 'Color', '#D95319');
plot(tElev, rad2deg(elevData.delta_e(:,2)), '--', 'Color', '#EDB120'); hold off;
grid on; legend('Coll [deg]','LonCyc [deg]','Elev [deg]','Interpreter', 'latex', 'FontSize', 10, 'Location', 'NorthWest')
ylim([-1;17])
xlim([0;5])

subplot(4,1,2)
plot(tCyc, cycData.u(:,2), tCyc, cycData.v(:,2), tCyc, cycData.w(:,2),'LineWidth', 1.5); hold on;
plot(tElev, elevData.u(:,2), '--', 'Color', '#0072BD');
plot(tElev, elevData.v(:,2), '--', 'Color', '#D95319');
plot(tElev, elevData.w(:,2), '--', 'Color', '#EDB120'); hold off;
grid on; legend('u [m/s]', 'v [m/s]', 'w [m/s]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'SouthWest')
xlim([0;5])

subplot(4,1,3)
plot(tCyc, cycData.p(:,2), tCyc, cycData.q(:,2), tCyc, cycData.r(:,2),'LineWidth', 1.5); grid on;
hold on;
plot(tElev, elevData.p(:,2), '--', 'Color', '#0072BD');
plot(tElev, elevData.q(:,2), '--', 'Color', '#D95319');
plot(tElev, elevData.r(:,2), '--', 'Color', '#EDB120'); hold off;
legend('p [deg/s]', 'q [deg/s]', 'r [deg/s]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'SouthWest')
xlim([0;5])

subplot(4,1,4)
plot(tCyc, cycData.phi(:,2), tCyc, cycData.theta(:,2), tCyc,cycData.psi(:,2), 'LineWidth', 1.5); grid on;
hold on;
plot(tElev, elevData.phi(:,2), '--', 'Color', '#0072BD');
plot(tElev, elevData.theta(:,2), '--', 'Color', '#D95319');
plot(tElev, elevData.psi(:,2), '--', 'Color', '#EDB120'); hold off;
legend('$\phi_f$ [deg]', '$\theta_f$ [deg]', '$\psi_f$ [deg]', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'SouthWest')
xlabel('Simulation runtime [s]', 'Interpreter', 'latex', 'FontSize', 13)
xlim([0;5])

sgtitle('Linear open-loop response at 100 m/s', 'Interpreter', 'latex')
f.Position = [500 200 570 650];

%% LOS
clc; close all; clear;
coaxial_heli_parameters;
LOS = 0.00002*V_vals.^2;

figure(27)
plot(V_vals, LOS, 'LineWidth', 2); grid on;
legend('LOS Schedule [Ferguson]', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'NorthWest')
xlabel('Forward speed [m/s]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('LOS [-]', 'Interpreter', 'latex', 'FontSize', 14)


%% Example 5 deg step pitch input for command model

% Parameters
t_end = 4;            % End time for simulation
t_step = 0.01;        % Time step for simulation
t = 0:t_step:t_end;   % Time vector
step_time = 1;        % Time when step input becomes active
step_value = 5;       % Value of the step input

% Create step input reference command
ref_command = zeros(size(t));
ref_command(t >= step_time) = step_value;

% Transfer function parameters
wn = 4;               % Natural frequency
zeta = 0.707;         % Damping ratio
num = wn^2;           % Numerator of transfer function
den = [1 2*zeta*wn wn^2];  % Denominator of transfer function

% Create transfer function system
sys = tf(num, den);

% Simulate system response to step input
[y, t_out] = lsim(sys, ref_command, t);

% Plot the reference command and system response
figure(28);
plot(t, ref_command, '--', 'LineWidth', 1.5, 'Color', '#0072BD'); hold on;
plot(t_out, y, '-', 'LineWidth', 2, 'Color', '#D95319');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\theta$ [deg]', 'Interpreter', 'latex', 'FontSize', 14);
ylim([0,6.5])
% title('Step Input Reference Command and System Response');
legend('Reference Command', 'EMF Model Response', 'Interpreter', 'latex', 'FontSize', 10);
grid on;



%% coupled step-input tuned plots
clear; clc; close all;

f = figure(29);
tunedData = load('stepTuningSimulink.mat');
tunedData = tunedData.out;

t = tunedData.tout;

plot(t, rad2deg(tunedData.phi_cmd(:,2)), '--', 'Color', '#0072BD', 'LineWidth', 1.5);  hold on;
plot(t, rad2deg(tunedData.theta_cmd(:,2)), '--', 'Color', '#D95319', 'LineWidth', 1.5);
plot(t, rad2deg(tunedData.psi_cmd(:,2)), '--', 'Color', '#EDB120', 'LineWidth', 1.5);

plot(t, rad2deg(tunedData.phi_c(:,2)), '-.', 'Color', '#0072BD', 'LineWidth', 2);
plot(t, rad2deg(tunedData.theta_c(:,2)), '-.',  'Color', '#D95319', 'LineWidth', 2);
plot(t, rad2deg(tunedData.psi_c(:,2)), '-.', 'Color', '#EDB120', 'LineWidth', 2);

g7 = plot(t, rad2deg(tunedData.phi(:,2)), 'Color', '#0072BD', 'LineWidth', 2);
g8 = plot(t, rad2deg(tunedData.theta(:,2)),  'Color', '#D95319', 'LineWidth', 2);
g9 = plot(t, rad2deg(tunedData.psi(:,2)),  'Color', '#EDB120', 'LineWidth', 2);

grid on; hold off;
legend([g7, g8, g9], {'$\phi$', '$\theta$', '$\psi$'}, 'Interpreter', 'latex', 'FontSize', 13);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Attitude [deg]', 'Interpreter', 'latex', 'FontSize', 14);

% performance value calculation
err_theta = sum(abs(tunedData.theta(:,2) - tunedData.theta_c(:,2)));
err_phi = sum(abs(tunedData.phi(:,2) - tunedData.phi_c(:,2)));
err_psi = sum(abs(tunedData.psi(:,2) - tunedData.psi_c(:,2)));

V = 1; %deg/in
err_rel_theta = 1/(V*t(end))*err_theta;
Q_theta = 1-err_rel_theta
err_rel_phi = 1/(V*t(end))*err_phi;
Q_phi = 1-err_rel_phi
err_rel_psi = 1/(V*t(end))*err_psi;
Q_psi = 1-err_rel_psi


%% Attitude Quickness
% Pitch

load('./attitudeQuicknessSaves/5degpitch.mat')
t = savedStates(1,:);
lonCyc5 = savedStates(2,:);
elev5 = savedStates(3,:);
q5 = savedStates(4,:);
theta5 = savedStates(5,:);

load('./attitudeQuicknessSaves/20degpitch.mat')
t = savedStates(1,:);
lonCyc20 = savedStates(2,:);
elev20 = savedStates(3,:);
q20 = savedStates(4,:);
theta20 = savedStates(5,:);


figure(30)
subplot(3,1,1)
plot(t(1:end-1), rad2deg(lonCyc5(1:end-1)), 'LineWidth', 1.5); hold on;
plot(t, rad2deg(elev5), 'LineWidth', 1.5); 
plot(t(1:end-1), rad2deg(lonCyc20(1:end-1)), '-.', 'Color', '#0072BD', 'LineWidth', 2); 
plot(t, rad2deg(elev20), '-.', 'Color', '#D95319', 'LineWidth', 2);
hold off;
grid on; 
legend('LonCyc [deg]', 'Elev [deg]', 'Interpreter', ...
    'latex', 'FontSize', 10, 'Location', 'SouthEast')
ylabel('Input [deg]', 'Interpreter', 'latex', 'FontSize', 14)

subplot(3,1,2)
plot(t, rad2deg(q5), 'LineWidth', 2); 
grid on; 
max_q_value = max(q5);
q_pk_deg = rad2deg(max_q_value)
disp(['Maximum Value of q: ', num2str(max_q_value), 'Interpreter', 'latex']);
hold on;
yline(rad2deg(max_q_value), '--r');
hold off;
ylabel('q [deg/s]', 'Interpreter', 'latex', 'FontSize', 14)
ylim([rad2deg(min(q5))*1.2, rad2deg(max_q_value)*1.4]);

subplot(3,1,3)
plot(t, rad2deg(theta5), 'LineWidth', 2); 
grid on;
ylim([rad2deg(min(theta5))*1.2, rad2deg(max(theta5)*1.4)])
ylabel('Attitude [deg]', 'Interpreter', 'latex', 'FontSize', 14)


% Calculate rise time for theta
percent_10 = 0.1; % 10%
percent_90 = 0.9; % 90%
steady_state_value_theta = mean(rad2deg(theta5(end-100:end))); % Assuming the steady state is the mean of last 100 samples
ten_percent_value_theta = percent_10 * steady_state_value_theta;
ninety_percent_value_theta = percent_90 * steady_state_value_theta;
idx_10_theta = find(rad2deg(theta5) > ten_percent_value_theta, 1, 'first');
idx_90_theta = find(rad2deg(theta5) > ninety_percent_value_theta, 1, 'first');
rise_time_theta = t(idx_90_theta) - t(idx_10_theta);

% Find maximum value of theta
max_value_theta = max(theta5);
theta_pk_deg = rad2deg(max_value_theta)

% Find minimum value after the rise time for theta
min_value_after_rise_theta = min(theta5(idx_90_theta:end));
rad2deg(min_value_after_rise_theta)

% Plot horizontal red dashed lines for theta
hold on;
yline(rad2deg(max_value_theta), '--r');
yline(rad2deg(min_value_after_rise_theta), '--r');
hold off;
legend('$\theta$ [deg]', 'Interpreter', 'latex', 'FontSize', 10, ...
    'Location', 'SouthEast')




























