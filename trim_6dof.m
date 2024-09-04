%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename:     coaxial_6dof_trim
% Project:      MSc Thesis
% Supervisor:   M.D. Pavel
% Author:       Ynias Prencipe 
% Student Nr.:  4777158
% 
% Description:  Trim using the Newton-Rhapson algorithm for the coaxial
% helicopter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Set-up
% Load coaxial helicopter parameters
coaxial_heli_parameters;

% Initial Values
theta_0(1) = deg2rad(10);
theta_d(1) = deg2rad(0);
theta_s(1) = deg2rad(0);
theta_c(1) = deg2rad(0);
theta_cdiff(1) = deg2rad(0);
theta_f(1) = deg2rad(0);
theta_p(1) = deg2rad(0);
delta_e(1) = deg2rad(0);
delta_r(1) = deg2rad(0);
t(1)=0;
u(1)=0;
w(1)=0;
p(1)=0;
q(1)=0;
r(1)=0;
phi(1)=0;
Omega_init = 40;
lambda_0_u(1)=sqrt(mass/2*abs(g)/(area*2*rho))/(Omega_init*R);
lambda_0_l(1)=sqrt(mass/2*abs(g)/(area*2*rho))/(Omega_init*R);
lambda_p(1) = 0.002;
vi_hov = sqrt(1/2*mass*g/(2*rho*pi*R^2));

% Define trim variable vector
x_k(:,1) = [theta_0(1); theta_d(1); theta_s(1); theta_c(1); phi(1); theta_p(1); ...
    lambda_0_u(1); lambda_0_l(1); lambda_p(1)];


%% Trim Routine

for i = 1:length(V_vals)
    i
    V(i) = V_vals(i);
    u(i) = V(i)*cos(atan2(w(i),u(i)));
    v(i) = 0;
    w(i) = V(i)*sin(atan2(w(i),u(i)));
    theta_f(i) = 0;
    delta_e(i) = 0;
    delta_r(i) = 0;
    % theta_cdiff(i) = 0;

    %% Rotor RPM Scheduling
    if V(i) < 70
        Omega = 40;
    else
        % Omega = 40 - 5/30*(V(i)-70);
        Omega = 40;
    end

    %%
    if i>1
        % x_k(:,i) = x_k(:,i-1);      % initial guess of each trim iteration is the previous trim state - for quick convergence
        x_k(:,i) = x_k(:,1);
    end

    %%
    if i>1
        theta_cdiff(i) = theta_cdiff(i-1);
    else
        theta_cdiff(i) = 0;
    end

    %% Define State Variable Names
    theta_0_u(i) = x_k(1,i)/2; theta_d(i)=x_k(2,i); theta_s(i) = x_k(3,i); 
    theta_c(i) = x_k(4,i); phi(i) = x_k(5,i); theta_p(i) = x_k(6,i); 
    lambda_0_u(i) = x_k(7,i); lambda_0_l(i) = x_k(8,i); lambda_0_p = x_k(9,i);
   
    %%
    % Calculate Angles
    alfa_sp(i) = atan2(w(i),u(i));
    alfa_cp(i) = alfa_sp(i) + theta_s(i);

    % Speed variables
    mu(i) = V(i)/(Omega*R);
    mu_x(i) = V(i)/(Omega*R) * cos(alfa_cp(i));
    mu_z(i) = V(i)/(Omega*R) * sin(alfa_cp(i));
    vel(:,i) = [u(i); v(i); w(i); mu(i); mu_x(i); mu_z(i)];

    %%
    lambda_c(i) = mu(i) * sin(alfa_cp(i));

    %%
    count(i) = 1;
    f_sum(1) = 10;
    % Newton-Rhapson Trim Algorithm 
    while f_sum(count(i)) > 0.0005 % convergence said to occur at 0.0005 according to leishman p 97

        count(i) = count(i)+1;
        tistep = deg2rad(0.01);          % small increment for Newton-Rhapson

        % calculate f(:,i) using x_k
        x_ki(:,i) = x_k(:,i);
        f(:,i) = f_xk6(vel(:,i), x_ki(:,i), p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));

        % calculate x_k+1 for perturbation in x(1,i)
        ti1 = [tistep; 0; 0; 0; 0; 0; 0; 0; 0];
        x_k1(:,i) = x_ki(:,i) + ti1;
        % calculate f_xk1_ti(:,i) using x_k+1 for FIRST column of jacobian
        f_xk1_ti(:,i) = f_xk6(vel(:,i), x_k1(:,i), p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));

        % calculate x_k+1 for perturbation in x(2,i)
        ti2 = [0; tistep; 0; 0; 0; 0; 0; 0; 0];
        x_k2(:,i) = x_ki(:,i) + ti2;
        % calculate f_xk2_ti(:,i) using x_k+1 for SECOND column of jacobian
        f_xk2_ti(:,i) = f_xk6(vel(:,i), x_k2(:,i), p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
        
        % calculate x_k+1 for perturbation in x(3,i)
        ti3 = [0; 0; tistep; 0; 0; 0; 0; 0; 0];
        x_k3(:,i) = x_ki(:,i) + ti3;
        % calculate f_xk_ti(:,i) using x_k+1 for THIRD column of jacobian
        f_xk3_ti(:,i) = f_xk6(vel(:,i), x_k3(:,i), p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));

        % calculate x_k+1 for perturbation in x(4,i)
        ti4 = [0; 0; 0; tistep; 0; 0; 0; 0; 0];
        x_k4(:,i) = x_ki(:,i) + ti4;
        f_xk4_ti(:,i) = f_xk6(vel(:,i), x_k4(:,i), p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));

        % calculate x_k+1 for perturbation in x(5,i)
        ti5 = [0; 0; 0; 0; tistep; 0; 0; 0; 0];
        x_k5(:,i) = x_ki(:,i) + ti5;
        f_xk5_ti(:,i) = f_xk6(vel(:,i), x_k5(:,i), p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));

        % calculate x_k+1 for perturbation in x(6,i)
        ti6 = [0; 0; 0; 0; 0; tistep; 0; 0; 0];
        x_k6(:,i) = x_ki(:,i) + ti6;
        f_xk6_ti(:,i) = f_xk6(vel(:,i), x_k6(:,i), p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));

        % calculate x_k+1 for perturbation in x(7,i)
        ti7 = [0; 0; 0; 0; 0; 0; tistep; 0; 0];
        x_k7(:,i) = x_ki(:,i) + ti7;
        f_xk7_ti(:,i) = f_xk6(vel(:,i), x_k7(:,i), p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));

        % calculate x_k+1 for perturbation in x(8,i)
        ti8 = [0; 0; 0; 0; 0; 0; 0; tistep; 0];
        x_k8(:,i) = x_ki(:,i) + ti8;
        f_xk8_ti(:,i) = f_xk6(vel(:,i), x_k8(:,i), p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));

        % calculate x_k+1 for perturbation in x(9,i)
        ti9 = [0; 0; 0; 0; 0; 0; 0; 0; tistep];
        x_k9(:,i) = x_ki(:,i) + ti9;
        f_xk9_ti(:,i) = f_xk6(vel(:,i), x_k9(:,i), p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
        
        % assemble all columns for Jacobian
        ti = tistep * ones(size(f(:,i)));
        Jac = [(f_xk1_ti(:,i)-f(:,i))./ti, (f_xk2_ti(:,i)-f(:,i))./ti, (f_xk3_ti(:,i)-f(:,i))./ti, ...
            (f_xk4_ti(:,i)-f(:,i))./ti, (f_xk5_ti(:,i)-f(:,i))./ti, (f_xk6_ti(:,i)-f(:,i))./ti, ...
            (f_xk7_ti(:,i)-f(:,i))./ti, (f_xk8_ti(:,i)-f(:,i))./ti, (f_xk9_ti(:,i)-f(:,i))./ti];
        
        % calculate new input vector according to Jacobian
        x_k(:,i) = x_ki(:,i) - inv(Jac) * f(:,i);

        % f(:,i)
        [f(:,i), output_other(:,i), Mw_components(:,i), angles(:,i), ...
            inflow(:,i), latAngles(:,i), LOSparams(:,i)] = ...
            f_xk6(vel(:,i), x_k(:,i), p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));

        f_sum(count(i)) = sum(abs(f(:,i)),"all");
        
    end

    % T_tot(i) = abs(LOSparams(1,i));
    % L_MR_u(i) = LOSparams(2,i);
    % L_MR_l(i) = LOSparams(3,i);
    % LOS(i) = L_MR_u(i)/(T_tot(i)*R);
    % LOSdes(i) = 0.00002*V(i)^2;
    % LOSdiff(i) = LOSdes(i) - LOS(i);
    % 
    % while abs(LOSdiff(i)) > 0.001
    %     theta_cdiff(i) = theta_cdiff(i) + deg2rad(0.01);
    %     [f(:,i), output_other(:,i), Mw_components(:,i), angles(:,i), ...
    %         inflow(:,i), latAngles(:,i), LOSparams(:,i)] = ...
    %         f_xk6(vel(:,i), x_k(:,i), p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    % 
    %     T_tot(i) = abs(LOSparams(1,i));
    %     L_MR_u(i) = LOSparams(2,i);
    %     L_MR_l(i) = LOSparams(3,i);
    %     LOS(i) = L_MR_u(i)/(T_tot(i)*R);
    %     LOSdes(i) = 0.00002*V(i)^2;
    %     LOSdiff(i) = LOSdes(i) - LOS(i);
    % end

    udot(i) = f(1,i);
    vdot(i) = f(2,i);
    wdot(i) = f(3,i);
    pdot(i) = f(4,i);
    qdot(i) = f(5,i);
    rdot(i) = f(6,i);
    lambda_0_u_dot(i) = f(7,i);
    lambda_0_l_dot(i) = f(8,i);
    lambda_p_dot(i) = f(9,i);

    u(i+1) = V(i)*cos(atan2(w(i),u(i)));
    w(i+1) = V(i)*sin(atan2(w(i),u(i)));
    p(i+1) = 0;
    q(i+1) = 0;
    r(i+1) = 0;

    ytrim(:,i) = [u(i); v(i); w(i); p(i); q(i); r(i)];
    

end

%%
% SAVE x_k(:,i)
if exist('trim_saved_6dof.mat', 'file') == 2
    load('trim_saved_6dof.mat', 'trim_saved_6dof');
    % Append the new x_k to the existing vector
    trim_saved_6dof = x_k;
else
    % Create a new vector to store x_k
    trim_saved_6dof = x_k;
end

% Save the updated x_k to a file
save('trim_saved_6dof.mat', 'trim_saved_6dof');

save('ytrim.mat', 'ytrim')

% save('rotorAngles_smallTwistMin3.mat', 'angles')

%% TRIM PLOT

figure(10)
theta_0 = rad2deg(x_k(1,:)); 
theta_d = rad2deg(x_k(2,:)); 
theta_s = rad2deg(x_k(3,:)); 
theta_c = rad2deg(x_k(4,:)); 
phi     = rad2deg(x_k(5,:)); 
theta_p = rad2deg(x_k(6,:));
subplot(5,2,1)
plot(mu, theta_0, '-*'), grid;
legend('\theta_0'); title('Mean Collective')
subplot(5,2,2)
plot(mu, theta_d, '-*'), grid;
legend('\theta_d'); title('Differential Collective')
subplot(5,2,3)
plot(mu, theta_s, '-*'), grid;
legend('\theta_s'); title('Longitudinal Cyclic')
subplot(5,2,4)
plot(mu, theta_c, '-*'), grid;
legend('\theta_c'); title('Lateral Cyclic')
subplot(5,2,5)
plot(mu, phi, '-*'), grid;
legend('\phi'); title('Roll Angle')
subplot(5,2,6)
plot(mu, theta_p, '-*'), grid;
legend('\theta_p'); title('Pusher Prop Collective')
subplot(5,2,7)
plot(mu, rad2deg(delta_e), '-*'), grid;
legend('\delta_e'); title('Elevator Deflection (Scheduled)')
subplot(5,2,8)
plot(mu, rad2deg(delta_r), '-*'), grid;
legend('\delta_r'); title('Rudder Deflection (Scheduled)')
subplot(5,2,9)
plot(mu, rad2deg(theta_f), '-*'), grid;
legend('\theta_f'); xlabel("Advance Ratio \mu [-]"); title('Fuselage Pitch (Scheduled)')
subplot(5,2,10)
plot(mu, zeros(length(theta_f)), '-*'), grid;
legend('-'); xlabel("Advance Ratio \mu [-]"); title('EMPTY')

%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ROTOR ANGLES PLOT
% 
% figure(11)
% a0_u = angles(1,:); a0_l = angles(2,:); a1_u = angles(3,:); 
% a1_l = angles(4,:); a1r_u = angles(5,:); a1r_l = angles(6,:); 
% alfa_sp = angles(7,:); alfa_cp = angles(8,:); alfa_dp_u = angles(9,:);
% alfa_dp_l = angles(10,:);
% plot(V_vals, a0_u, V_vals, a0_l, ...
%     V_vals, a1_u, V_vals, a1_l, ...
%     V_vals, a1r_u, 'k--', V_vals, a1r_l, 'k--', ...
%     V_vals, alfa_sp, 'g--o', V_vals, alfa_cp, 'b--o');
% grid on;
% legend('$a_{0_u}$', '$a_{0_l}$', '$a_{1_u}$', '$a_{1_l}$', '$a_{1r_u}$', '$a_{1r_l}$', ...
%     '$\alpha_{sp}$', '$\alpha_{cp}$', 'Interpreter', 'latex');
% ylabel('Angle [deg]', 'Interpreter', 'latex')
% xlabel('Forward velocity [m/s]', 'Interpreter', 'latex')
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Linearisation 
lin_step = 0.1;
lin_systems = cell(length(V_vals), 2);
% at V = 0.1m/s
for i = 1:length(V_vals)
    i
    trim_vals_0 = x_k(:,i);

    %%%%%%%%%%%%%%%%%%% A-MATRIX %%%%%%%%%%%%%%%%%%%
    % Row 1 (X-force --> udot state)
    vel_lin0 = [lin_step; 0; 0; 0; 0; 0];
    f_lin_dist_0_pos = f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_u(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; lin_step; 0; 0; 0; 0];
    f_lin_dist_0_pos = f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_v(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; 0; lin_step; 0; 0; 0];
    f_lin_dist_0_pos = f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_w(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i)+0.01, q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i)-0.01, q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_p(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i)+0.01, r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i)-0.01, r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_q(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i)+0.01, theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i)-0.01, theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_r(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    % Row 2 (Z-force --> vdot state)
    vel_lin0 = [lin_step; 0; 0; 0; 0; 0];
    f_lin_dist_0_pos = f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Y_u(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; lin_step; 0; 0; 0; 0];
    f_lin_dist_0_pos = f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Y_v(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; 0; lin_step; 0; 0; 0];
    f_lin_dist_0_pos = f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Y_w(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i)+0.01, q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i)-0.01, q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Y_p(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i)+0.01, r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i)-0.01, r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Y_q(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i)+0.01, theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i)-0.01, theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Y_r(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    % Row 3 (Z-force --> wdot state)
    vel_lin0 = [lin_step; 0; 0; 0; 0; 0];
    f_lin_dist_0_pos = f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    Z_u(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; lin_step; 0; 0; 0; 0];
    f_lin_dist_0_pos = f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    Z_v(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; 0; lin_step; 0; 0; 0];
    f_lin_dist_0_pos = f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    Z_w(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i)+0.01, q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i)-0.01, q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    Z_p(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i)+0.01, r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i)-0.01, r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    Z_q(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i)+0.01, theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i)-0.01, theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    Z_r(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    % Row 4 (L-moment --> pdot state)
    vel_lin0 = [lin_step; 0; 0; 0; 0; 0];
    f_lin_dist_0_pos = f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(4,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(4,1);
    L_u(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; lin_step; 0; 0; 0; 0];
    f_lin_dist_0_pos = f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(4,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(4,1);
    L_v(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; 0; lin_step; 0; 0; 0];
    f_lin_dist_0_pos = f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(4,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(4,1);
    L_w(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i)+0.01, q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(4,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i)-0.01, q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(4,1);
    L_p(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i)+0.01, r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(4,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i)-0.01, r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(4,1);
    L_q(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i)+0.01, theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(4,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i)-0.01, theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(4,1);
    L_r(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    % Row 5 (M-moment --> qdot state)
    vel_lin0 = [lin_step; 0; 0; 0; 0; 0];
    f_lin_dist_0_pos = f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(5,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(5,1);
    M_u(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; lin_step; 0; 0; 0; 0];
    f_lin_dist_0_pos = f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(5,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(5,1);
    M_v(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; 0; lin_step; 0; 0; 0];
    [f_lin_dist_0_pos, M_moments_Mw_pos, Mw_components_pos] = ...
        f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos_Mw = f_lin_dist_0_pos(5,1);
    [f_lin_dist_0_neg, M_moments_Mw_neg, Mw_components_neg] = ...
        f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg_Mw = f_lin_dist_0_neg(5,1);
    M_w(i) = (f_lin_dist_0_pos_Mw - f_lin_dist_0_neg_Mw)/(2*lin_step);
    M_w_components(:,i) = (Mw_components_pos - Mw_components_neg)./(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i)+0.01, q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(5,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i)-0.01, q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(5,1);
    M_p(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i)+0.01, r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(5,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i)-0.01, r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(5,1);
    M_q(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i)+0.01, theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(5,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i)-0.01, theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(5,1);
    M_r(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    % Row 6 (N-moment --> rdot state)
    vel_lin0 = [lin_step; 0; 0; 0; 0; 0];
    f_lin_dist_0_pos = f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(6,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(6,1);
    N_u(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; lin_step; 0; 0; 0; 0];
    f_lin_dist_0_pos = f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(6,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(6,1);
    N_v(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; 0; lin_step; 0; 0; 0];
    f_lin_dist_0_pos = f_xk6(vel(:,i)+vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(6,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i)-vel_lin0, trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(6,1);
    N_w(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i)+0.01, q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(6,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i)-0.01, q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(6,1);
    N_p(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i)+0.01, r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(6,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i)-0.01, r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(6,1);
    N_q(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i)+0.01, theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(6,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i)-0.01, theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(6,1);
    N_r(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    % Assemble A-matrix
    A_lin = [X_u(i), X_v(i), X_w(i), X_p(i), X_q(i), X_r(i);
            Y_u(i), Y_v(i), Y_w(i), Y_p(i), Y_q(i), Y_r(i);
            Z_u(i), Z_v(i), Z_w(i), Z_p(i), Z_q(i), Z_r(i);
            L_u(i), L_v(i), L_w(i), L_p(i), L_q(i), L_r(i);
            M_u(i), M_v(i), M_w(i), M_p(i), M_q(i), M_r(i);
            N_u(i), N_v(i), N_w(i), N_p(i), N_q(i), N_r(i)]

    % %%%%%%%%%%%%%%%%%%% B-MATRIX %%%%%%%%%%%%%%%%%%%
    lin_step = deg2rad(0.1);

    % Row 1 (X-force)
    b_step = [lin_step; 0; 0; 0; 0; 0; 0; 0; 0];   % input in collective
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_theta_0(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; lin_step; 0; 0; 0; 0; 0; 0; 0];   % input in differential collective
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_theta_d(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; lin_step; 0; 0; 0; 0; 0; 0]; % input in longitudinal cyclic
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_theta_s(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; 0; lin_step; 0; 0; 0; 0; 0]; % input in lateral cyclic
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_theta_c(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i)+lin_step);
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i)-lin_step);
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_theta_cdiff(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; 0; 0; 0; lin_step; 0; 0; 0]; % input in lateral cyclic
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_theta_p(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i)+lin_step, delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i)-lin_step, delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_delta_e(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i)+lin_step, theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i)-lin_step, theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_delta_r(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    % Row 2 (Y-force)
    b_step = [lin_step; 0; 0; 0; 0; 0; 0; 0; 0];   % input in collective
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Y_theta_0(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; lin_step; 0; 0; 0; 0; 0; 0; 0];   % input in differential collective
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Y_theta_d(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; lin_step; 0; 0; 0; 0; 0; 0]; % input in longitudinal cyclic
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Y_theta_s(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; 0; lin_step; 0; 0; 0; 0; 0]; % input in lateral cyclic
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Y_theta_c(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i)+lin_step);
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i)-lin_step);
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Y_theta_cdiff(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; 0; 0; 0; lin_step; 0; 0; 0]; % input in prop coll
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Y_theta_p(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i)+lin_step, delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i)-lin_step, delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Y_delta_e(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i)+lin_step, theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i)-lin_step, theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Y_delta_r(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);
    
    % Row 3 (Z-Force)
    b_step = [lin_step; 0; 0; 0; 0; 0; 0; 0; 0];   % input in collective
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    Z_theta_0(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; lin_step; 0; 0; 0; 0; 0; 0; 0];   % input in differential collective
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    Z_theta_d(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; lin_step; 0; 0; 0; 0; 0; 0]; % input in longitudinal cyclic
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    Z_theta_s(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; 0; lin_step; 0; 0; 0; 0; 0]; % input in lateral cyclic
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    Z_theta_c(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i)+lin_step);
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i)-lin_step);
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    Z_theta_cdiff(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; 0; 0; 0; lin_step; 0; 0; 0]; % input in prop coll
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    Z_theta_p(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i)+lin_step, delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i)-lin_step, delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    Z_delta_e(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i)+lin_step, theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i)-lin_step, theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    Z_delta_r(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    % Row 4 (L-Moment)
    b_step = [lin_step; 0; 0; 0; 0; 0; 0; 0; 0];   % input in collective
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(4,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(4,1);
    L_theta_0(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; lin_step; 0; 0; 0; 0; 0; 0; 0];   % input in differential collective
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(4,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(4,1);
    L_theta_d(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; lin_step; 0; 0; 0; 0; 0; 0]; % input in longitudinal cyclic
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(4,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(4,1);
    L_theta_s(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; 0; lin_step; 0; 0; 0; 0; 0]; % input in lateral cyclic
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(4,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(4,1);
    L_theta_c(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i)+lin_step);
    f_lin_dist_0_pos = f_lin_dist_0_pos(4,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i)-lin_step);
    f_lin_dist_0_neg = f_lin_dist_0_neg(4,1);
    L_theta_cdiff(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; 0; 0; 0; lin_step; 0; 0; 0]; % input in prop coll
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(4,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(4,1);
    L_theta_p(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i)+lin_step, delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(4,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i)-lin_step, delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(4,1);
    L_delta_e(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i)+lin_step, theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(4,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i)-lin_step, theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(4,1);
    L_delta_r(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    % Row 5 (M-Moment)
    b_step = [lin_step; 0; 0; 0; 0; 0; 0; 0; 0];   % input in collective
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(5,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(5,1);
    M_theta_0(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; lin_step; 0; 0; 0; 0; 0; 0; 0];   % input in differential collective
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(5,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(5,1);
    M_theta_d(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; lin_step; 0; 0; 0; 0; 0; 0]; % input in longitudinal cyclic
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(5,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(5,1);
    M_theta_s(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; 0; lin_step; 0; 0; 0; 0; 0]; % input in lateral cyclic
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(5,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(5,1);
    M_theta_c(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i)+lin_step);
    f_lin_dist_0_pos = f_lin_dist_0_pos(5,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i)-lin_step);
    f_lin_dist_0_neg = f_lin_dist_0_neg(5,1);
    M_theta_cdiff(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; 0; 0; 0; lin_step; 0; 0; 0]; % input in prop coll
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(5,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(5,1);
    M_theta_p(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i)+lin_step, delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(5,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i)-lin_step, delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(5,1);
    M_delta_e(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i)+lin_step, theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(5,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i)-lin_step, theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(5,1);
    M_delta_r(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    % Row 6 (N-Moment)
    b_step = [lin_step; 0; 0; 0; 0; 0; 0; 0; 0];   % input in collective
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(6,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(6,1);
    N_theta_0(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; lin_step; 0; 0; 0; 0; 0; 0; 0];   % input in differential collective
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(6,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(6,1);
    N_theta_d(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; lin_step; 0; 0; 0; 0; 0; 0]; % input in longitudinal cyclic
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(6,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(6,1);
    N_theta_s(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; 0; lin_step; 0; 0; 0; 0; 0]; % input in lateral cyclic
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(6,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(6,1);
    N_theta_c(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i)+lin_step);
    f_lin_dist_0_pos = f_lin_dist_0_pos(6,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i)-lin_step);
    f_lin_dist_0_neg = f_lin_dist_0_neg(6,1);
    N_theta_cdiff(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; 0; 0; 0; lin_step; 0; 0; 0]; % input in prop coll
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0+b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(6,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0-b_step, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(6,1);
    N_theta_p(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i)+lin_step, delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(6,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i)-lin_step, delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(6,1);
    N_delta_e(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i)+lin_step, theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(6,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i), delta_e(i), delta_r(i)-lin_step, theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(6,1);
    N_delta_r(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);
    
    % Theta X and Z force influence
    f_lin_dist_0_pos = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i)+lin_step, delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk6(vel(:,i), trim_vals_0, p(i), q(i), r(i), theta_f(i)-lin_step, delta_e(i), delta_r(i), theta_cdiff(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    X_theta_f(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    % Assemble B-MATRIX
    B_lin = [X_theta_0(i), X_theta_d(i), X_theta_s(i), X_theta_c(i), X_theta_cdiff(i), X_theta_p(i), X_delta_e(i), X_delta_r(i);
            Y_theta_0(i), Y_theta_d(i), Y_theta_s(i), Y_theta_c(i), Y_theta_cdiff(i), Y_theta_p(i), Y_delta_e(i), Y_delta_r(i);
            Z_theta_0(i), Z_theta_d(i), Z_theta_s(i), Z_theta_c(i), Z_theta_cdiff(i), Z_theta_p(i), Z_delta_e(i), Z_delta_r(i);
            L_theta_0(i), L_theta_d(i), L_theta_s(i), L_theta_c(i), L_theta_cdiff(i), L_theta_p(i), L_delta_e(i), L_delta_r(i);
            M_theta_0(i), M_theta_d(i), M_theta_s(i), M_theta_c(i), M_theta_cdiff(i), M_theta_p(i), M_delta_e(i), M_delta_r(i);
            N_theta_0(i), N_theta_d(i), N_theta_s(i), N_theta_c(i), N_theta_cdiff(i), N_theta_p(i), N_delta_e(i), N_delta_r(i)]

    % Assemble C-Matrix
    C_lin = eye(size(A_lin));

    % Assemble D-Matrix
    D_lin = zeros(size(B_lin));       % #columns = #inputs


    % Save Current Trim SS System in lin_systems
    TrimLin_6dof{i, 1} = A_lin;
    TrimLin_6dof{i, 2} = B_lin;
    TrimLin_6dof{i, 3} = C_lin;
    TrimLin_6dof{i, 4} = D_lin;
end

sys6DoF = ss(TrimLin_6dof{9, 1}, TrimLin_6dof{9, 2}, TrimLin_6dof{9, 3}, ...
    TrimLin_6dof{9, 4});

% Save Linear Systems around Trim states
if exist('TrimLinSave_6dof.mat', 'file') == 2
    load('TrimLinSave_6dof.mat', 'TrimLinSave_6dof');
    % Append the new x_k to the existing vector
    TrimLinSave_6dof = TrimLin_6dof;
else
    % Create a new vector to store x_k
    TrimLinSave_6dof = TrimLin_6dof;
end

% Save the updated x_k to a file
save('TrimLinSave_6dof.mat', 'TrimLinSave_6dof');
save('A_matrix_6dof.mat', 'A_lin');
save('B_matrix_6dof.mat', 'B_lin');

%%
stabDerivs = [X_u; X_v; X_w; X_p; X_q; X_r;...
    Y_u; Y_v; Y_w; Y_p; Y_q; Y_r;...
    Z_u; Z_v; Z_w; Z_p; Z_q; Z_r;...
    L_u; L_v; L_w; L_p; L_q; L_r;...
    M_u; M_v; M_w; M_p; M_q; M_r;...
    N_u; N_v; N_w; N_p; N_q; N_r];
save('stabDerivs.mat', 'stabDerivs');
contrDerivs = [X_theta_0; X_theta_d; X_theta_s; X_theta_c; X_theta_cdiff; X_theta_p; X_delta_e; X_delta_r; ...
    Y_theta_0; Y_theta_d; Y_theta_s; Y_theta_c; Y_theta_cdiff; Y_theta_p; Y_delta_e; Y_delta_r; ...
    Z_theta_0; Z_theta_d; Z_theta_s; Z_theta_c; Z_theta_cdiff; Z_theta_p; Z_delta_e; Z_delta_r; ...
    L_theta_0; L_theta_d; L_theta_s; L_theta_c; L_theta_cdiff; L_theta_p; L_delta_e; L_delta_r; ...
    M_theta_0; M_theta_d; M_theta_s; M_theta_c; M_theta_cdiff; M_theta_p; M_delta_e; M_delta_r; ...
    N_theta_0; N_theta_d; N_theta_s; N_theta_c; N_theta_cdiff; N_theta_p; N_delta_e; N_delta_r];
save('contrDerivs.mat', 'contrDerivs');

V_vals_kts = V_vals.*1.944;

%% STABILITY DERIVATIVES PLOTS
% figure(10)
% subplot(3,3,1)
% plot(V_vals_kts, X_u, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_u'); legend('X_u')
% subplot(3,3,2)
% plot(V_vals_kts, X_w, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_w'); legend('X_w')
% subplot(3,3,3)
% plot(V_vals_kts, X_q, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_q'); legend('X_q')
% subplot(3,3,4)
% plot(V_vals_kts, Z_u, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_u'); legend('Z_u')
% subplot(3,3,5)
% plot(V_vals_kts, Z_w, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_w'); legend('Z_w')
% subplot(3,3,6)
% plot(V_vals_kts, Z_q, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_q'); legend('Z_q')
% subplot(3,3,7)
% plot(V_vals_kts, M_u, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_u'); legend('M_u')
% subplot(3,3,8)
% plot(V_vals_kts, M_w, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_w'); legend('M_w')
% subplot(3,3,9)
% plot(V_vals_kts, M_q, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_q'); legend('M_q') 
% 
% figure(20)
% subplot(6,6,1)
% plot(V_vals_kts, X_u, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_u'); legend('X_u')
% subplot(6,6,2)
% plot(V_vals_kts, X_v, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_v'); legend('X_v')
% subplot(6,6,3)
% plot(V_vals_kts, X_w, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_w'); legend('X_w')
% subplot(6,6,4)
% plot(V_vals_kts, X_p, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_p'); legend('X_p')
% subplot(6,6,5)
% plot(V_vals_kts, X_q, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_q'); legend('X_q')
% subplot(6,6,6)
% plot(V_vals_kts, X_r, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_r'); legend('X_r')
% 
% subplot(6,6,7)
% plot(V_vals_kts, Y_u, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Y_u'); legend('Y_u')
% subplot(6,6,8)
% plot(V_vals_kts, Y_v, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Y_v'); legend('Y_v')
% subplot(6,6,9)
% plot(V_vals_kts, Y_w, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Y_w'); legend('Y_w')
% subplot(6,6,10)
% plot(V_vals_kts, Y_p, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Y_p'); legend('Y_p')
% subplot(6,6,11)
% plot(V_vals_kts, Y_q, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Y_q'); legend('Y_q')
% subplot(6,6,12)
% plot(V_vals_kts, Y_r, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Y_r'); legend('Y_r')
% 
% subplot(6,6,13)
% plot(V_vals_kts, Z_u, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_u'); legend('Z_u')
% subplot(6,6,14)
% plot(V_vals_kts, Z_v, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_v'); legend('Z_v')
% subplot(6,6,15)
% plot(V_vals_kts, Z_w, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_w'); legend('Z_w')
% subplot(6,6,16)
% plot(V_vals_kts, Z_p, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_p'); legend('Z_p')
% subplot(6,6,17)
% plot(V_vals_kts, Z_q, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_q'); legend('Z_q')
% subplot(6,6,18)
% plot(V_vals_kts, Z_r, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_r'); legend('Z_r')
% 
% subplot(6,6,19)
% plot(V_vals_kts, L_u, '--*'); grid on; xlabel('V_f [kts]'); ylabel('L_u'); legend('L_u')
% subplot(6,6,20)
% plot(V_vals_kts, L_v, '--*'); grid on; xlabel('V_f [kts]'); ylabel('L_v'); legend('L_v')
% subplot(6,6,21)
% plot(V_vals_kts, L_w, '--*'); grid on; xlabel('V_f [kts]'); ylabel('L_w'); legend('L_w')
% subplot(6,6,22)
% plot(V_vals_kts, L_p, '--*'); grid on; xlabel('V_f [kts]'); ylabel('L_p'); legend('L_p')
% subplot(6,6,23)
% plot(V_vals_kts, L_q, '--*'); grid on; xlabel('V_f [kts]'); ylabel('L_q'); legend('L_q')
% subplot(6,6,24)
% plot(V_vals_kts, L_r, '--*'); grid on; xlabel('V_f [kts]'); ylabel('L_r'); legend('L_r')
% 
% subplot(6,6,25)
% plot(V_vals_kts, M_u, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_u'); legend('M_u')
% subplot(6,6,26)
% plot(V_vals_kts, M_v, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_v'); legend('M_v')
% subplot(6,6,27)
% plot(V_vals_kts, M_w, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_w'); legend('M_w')
% subplot(6,6,28)
% plot(V_vals_kts, M_p, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_p'); legend('M_p')
% subplot(6,6,29)
% plot(V_vals_kts, M_q, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_q'); legend('M_q')
% subplot(6,6,30)
% plot(V_vals_kts, M_r, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_r'); legend('M_r')
% 
% subplot(6,6,31)
% plot(V_vals_kts, N_u, '--*'); grid on; xlabel('V_f [kts]'); ylabel('N_u'); legend('N_u')
% subplot(6,6,32)
% plot(V_vals_kts, N_v, '--*'); grid on; xlabel('V_f [kts]'); ylabel('N_v'); legend('N_v')
% subplot(6,6,33)
% plot(V_vals_kts, N_w, '--*'); grid on; xlabel('V_f [kts]'); ylabel('N_w'); legend('N_w')
% subplot(6,6,34)
% plot(V_vals_kts, N_p, '--*'); grid on; xlabel('V_f [kts]'); ylabel('N_p'); legend('N_p')
% subplot(6,6,35)
% plot(V_vals_kts, N_q, '--*'); grid on; xlabel('V_f [kts]'); ylabel('N_q'); legend('N_q')
% subplot(6,6,36)
% plot(V_vals_kts, N_r, '--*'); grid on; xlabel('V_f [kts]'); ylabel('N_r'); legend('N_r')

%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CONTROL DERIVATIVES PLOTS
% 
% figure(21)
% subplot(6,8,1)
% plot(V_vals_kts, X_theta_0, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_theta_0'); legend('X_theta_0')
% subplot(6,8,2)
% plot(V_vals_kts, X_theta_d, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_theta_d'); legend('X_theta_d')
% subplot(6,8,3)
% plot(V_vals_kts, X_theta_s, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_theta_s'); legend('X_theta_s')
% subplot(6,8,4)
% plot(V_vals_kts, X_theta_c, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_delta_c'); legend('X_delta_c')
% subplot(6,8,5)
% plot(V_vals_kts, X_theta_cdiff, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_delta_cdiff'); legend('X_delta_cdiff')
% subplot(6,8,6)
% plot(V_vals_kts, X_theta_p, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_theta_p'); legend('X_theta_p')
% subplot(6,8,7)
% plot(V_vals_kts, X_delta_e, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_delta_e'); legend('X_delta_e')
% subplot(6,8,8)
% plot(V_vals_kts, X_delta_r, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_delta_r'); legend('X_delta_r')
% 
% subplot(6,8,9)
% plot(V_vals_kts, Y_theta_0, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Y_theta_0'); legend('Y_theta_0')
% subplot(6,8,10)
% plot(V_vals_kts, Y_theta_d, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Y_theta_d'); legend('Y_theta_d')
% subplot(6,8,11)
% plot(V_vals_kts, Y_theta_s, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Y_theta_s'); legend('Y_theta_s')
% subplot(6,8,12)
% plot(V_vals_kts, Y_theta_c, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Y_delta_c'); legend('Y_delta_c')
% subplot(6,8,13)
% plot(V_vals_kts, Y_theta_cdiff, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Y_delta_cdiff'); legend('Y_delta_cdiff')
% subplot(6,8,14)
% plot(V_vals_kts, Y_theta_p, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Y_theta_p'); legend('Y_theta_p')
% subplot(6,8,15)
% plot(V_vals_kts, Y_delta_e, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Y_delta_e'); legend('Y_delta_e')
% subplot(6,8,16)
% plot(V_vals_kts, Y_delta_r, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Y_delta_r'); legend('Y_delta_r')
% 
% subplot(6,8,17)
% plot(V_vals_kts, Z_theta_0, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_theta_0'); legend('Z_theta_0')
% subplot(6,8,18)
% plot(V_vals_kts, Z_theta_d, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_theta_d'); legend('Z_theta_d')
% subplot(6,8,19)
% plot(V_vals_kts, Z_theta_s, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_theta_s'); legend('Z_theta_s')
% subplot(6,8,20)
% plot(V_vals_kts, Z_theta_c, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_delta_c'); legend('Z_delta_c')
% subplot(6,8,21)
% plot(V_vals_kts, Z_theta_cdiff, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_delta_cdiff'); legend('Z_delta_cdiff')
% subplot(6,8,22)
% plot(V_vals_kts, Z_theta_p, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_theta_p'); legend('Z_theta_p')
% subplot(6,8,23)
% plot(V_vals_kts, Z_delta_e, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_delta_e'); legend('Z_delta_e')
% subplot(6,8,24)
% plot(V_vals_kts, Z_delta_r, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_delta_r'); legend('Z_delta_r')
% 
% subplot(6,8,25)
% plot(V_vals_kts, L_theta_0, '--*'); grid on; xlabel('V_f [kts]'); ylabel('L_theta_0'); legend('L_theta_0')
% subplot(6,8,26)
% plot(V_vals_kts, L_theta_d, '--*'); grid on; xlabel('V_f [kts]'); ylabel('L_theta_d'); legend('L_theta_d')
% subplot(6,8,27)
% plot(V_vals_kts, L_theta_s, '--*'); grid on; xlabel('V_f [kts]'); ylabel('L_theta_s'); legend('L_theta_s')
% subplot(6,8,28)
% plot(V_vals_kts, L_theta_c, '--*'); grid on; xlabel('V_f [kts]'); ylabel('L_delta_c'); legend('L_delta_c')
% subplot(6,8,29)
% plot(V_vals_kts, L_theta_cdiff, '--*'); grid on; xlabel('V_f [kts]'); ylabel('L_delta_cdiff'); legend('L_delta_cdiff')
% subplot(6,8,30)
% plot(V_vals_kts, L_theta_p, '--*'); grid on; xlabel('V_f [kts]'); ylabel('L_theta_p'); legend('L_theta_p')
% subplot(6,8,31)
% plot(V_vals_kts, L_delta_e, '--*'); grid on; xlabel('V_f [kts]'); ylabel('L_delta_e'); legend('L_delta_e')
% subplot(6,8,32)
% plot(V_vals_kts, L_delta_r, '--*'); grid on; xlabel('V_f [kts]'); ylabel('L_delta_r'); legend('L_delta_r')
% 
% subplot(6,8,33)
% plot(V_vals_kts, M_theta_0, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_theta_0'); legend('M_theta_0')
% subplot(6,8,34)
% plot(V_vals_kts, M_theta_d, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_theta_d'); legend('M_theta_d')
% subplot(6,8,35)
% plot(V_vals_kts, M_theta_s, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_theta_s'); legend('M_theta_s')
% subplot(6,8,36)
% plot(V_vals_kts, M_theta_c, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_delta_c'); legend('M_delta_c')
% subplot(6,8,37)
% plot(V_vals_kts, M_theta_cdiff, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_delta_cdiff'); legend('M_delta_cdiff')
% subplot(6,8,38)
% plot(V_vals_kts, M_theta_p, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_theta_p'); legend('M_theta_p')
% subplot(6,8,39)
% plot(V_vals_kts, M_delta_e, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_delta_e'); legend('M_delta_e')
% subplot(6,8,40)
% plot(V_vals_kts, M_delta_r, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_delta_r'); legend('M_delta_r')
% 
% subplot(6,8,41)
% plot(V_vals_kts, N_theta_0, '--*'); grid on; xlabel('V_f [kts]'); ylabel('N_theta_0'); legend('N_theta_0')
% subplot(6,8,42)
% plot(V_vals_kts, N_theta_d, '--*'); grid on; xlabel('V_f [kts]'); ylabel('N_theta_d'); legend('N_theta_d')
% subplot(6,8,43)
% plot(V_vals_kts, N_theta_s, '--*'); grid on; xlabel('V_f [kts]'); ylabel('N_theta_s'); legend('N_theta_s')
% subplot(6,8,44)
% plot(V_vals_kts, N_theta_c, '--*'); grid on; xlabel('V_f [kts]'); ylabel('N_delta_c'); legend('N_delta_c')
% subplot(6,8,45)
% plot(V_vals_kts, N_theta_cdiff, '--*'); grid on; xlabel('V_f [kts]'); ylabel('N_delta_cdiff'); legend('N_delta_cdiff')
% subplot(6,8,46)
% plot(V_vals_kts, N_theta_p, '--*'); grid on; xlabel('V_f [kts]'); ylabel('N_theta_p'); legend('N_theta_p')
% subplot(6,8,47)
% plot(V_vals_kts, N_delta_e, '--*'); grid on; xlabel('V_f [kts]'); ylabel('N_delta_e'); legend('N_delta_e')
% subplot(6,8,48)
% plot(V_vals_kts, N_delta_r, '--*'); grid on; xlabel('V_f [kts]'); ylabel('N_delta_r'); legend('N_delta_r')

%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CONTROL EFFECTIVENESS PLOTS
% % Find the index where V_vals_kts is closest to 10
% transition_min_val = 20;        % [m/s]
% [~, idx_transstart] = min(abs(V - transition_min_val));
% 
% % Find the intersection point between M_theta_c and M_delta_e
% [~, idx_intersection] = min(abs(abs(M_theta_s) - abs(M_delta_e)));
% 
% % Create the plot
% figure(91);
% plot(V_vals_kts, abs(M_theta_s), V_vals_kts, abs(M_delta_e), '--*');
% grid on;
% xlabel('V_f [kts]');
% ylabel('CE_\theta (1/s^2)');
% hold on;
% % Plot vertical lines at x = 40 and x = 121.6
% 
% V_l_transit = V_vals_kts(idx_transstart);
% V_u_transit = V_vals_kts(idx_intersection);
% V_transit = [V_l_transit, V_u_transit];
% 
% % Save Transition Speeds
% if exist('V_transit_save.mat', 'file') == 2
%     load('V_transit_save.mat', 'V_transit_save');
%     % Append the new x_k to the existing vector
%     V_transit_save = V_transit;
% else
%     % Create a new vector to store x_k
%     V_transit_save = V_transit;
% end
% 
% % Save the updated x_k to a file
% save('V_transit_save.mat', 'V_transit_save');
% 
% line([V_l_transit, V_l_transit], ylim, 'Color', 'r', 'LineStyle', '--');
% line([V_u_transit, V_u_transit], ylim, 'Color', 'r', 'LineStyle', '--');
% % Shade the region between the two lines
% x_shade = [V_vals_kts(idx_transstart), V_vals_kts(idx_intersection)];
% y_shade = ylim; % Use ylim to define y-values for the shaded region
% fill([x_shade, fliplr(x_shade)], [y_shade(1), y_shade(1), y_shade(2), y_shade(2)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% % Add text 'Transition region' at the middle of the shaded region
% text(mean(x_shade), 190, 'Transition region', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'FontSize', 11); % Adjust the font weight and size of the text
% legend('M_{\theta_s}', 'M_{\delta_e}', 'Location', 'southeast', 'FontSize', 12);
% hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STACKED BAR PLOT ON M_W_COMPONENTS THROUGHOUT ENTIRE FORWARD VELOCITY SPEED RANGE
% 
% % figure(21)
% % hp = bar(V_vals_kts, M_w_components', 'stacked');
% % cm = colororder;
% % xlabel('Forward speed u [m/s]'); ylabel('M_w [rad/s m]');
% % legend('M_{MR_u}', 'M_{MR_l}', 'M_{hinge}', 'M_{flap}', 'M_{fus}', 'M_{ht}', 'M_e');
% % title('Contributions to the static stability derivative M_w');
% % grid on;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%

%% M_W COMPONENTS BAR PLOT FOR SELECTED TRIM SPEEDS -- FOR COMPARISON WITH NLR PLOTS
% 
M_w_MR_u = M_w_components(1,:); M_w_MR_l = M_w_components(2,:); 
M_w_hinge = M_w_components(3,:); M_w_flap = M_w_components(4,:); 
M_w_fus = M_w_components(5,:); M_w_ht = M_w_components(6,:);

M_w_rot = M_w_MR_u + M_w_MR_l + M_w_hinge + M_w_flap;

figure(22)
% bary = [M_w(1), M_w(6), M_w(11), M_w(16), M_w(21), M_w(26); 
%     M_w_MR_u(1), M_w_MR_u(6), M_w_MR_u(11), M_w_MR_u(16), M_w_MR_u(21), M_w_MR_u(26);
%     M_w_MR_l(1), M_w_MR_l(6), M_w_MR_l(11), M_w_MR_l(16), M_w_MR_l(21), M_w_MR_l(26);
%     M_w_hinge(1), M_w_hinge(6), M_w_hinge(11), M_w_hinge(16), M_w_hinge(21), M_w_hinge(26);
%     M_w_flap(1), M_w_flap(6), M_w_flap(11), M_w_flap(16), M_w_flap(21), M_w_flap(26);
%     M_w_fus(1), M_w_fus(6), M_w_fus(11), M_w_fus(16), M_w_fus(21), M_w_fus(26);
%     M_w_ht(1), M_w_ht(6), M_w_ht(11), M_w_ht(16), M_w_ht(21), M_w_ht(26)];
bary = [M_w(1), M_w(6), M_w(11), M_w(16), M_w(21), M_w(26); 
    M_w_rot(1), M_w_rot(6), M_w_rot(11), M_w_rot(16), M_w_rot(21), M_w_rot(26);
    M_w_fus(1), M_w_fus(6), M_w_fus(11), M_w_fus(16), M_w_fus(21), M_w_fus(26);
    M_w_ht(1), M_w_ht(6), M_w_ht(11), M_w_ht(16), M_w_ht(21), M_w_ht(26)];
b = bar(bary);
s=6;
colours = [];
for j = 1:numel(b)
    colours = [colours; [1-j/s, j/s, 1]]; % Adjust the color definition as needed
    set(b(j), 'FaceColor', colours(j,:))
end
legend('0 m/s', '25 m/s', '50 m/s', '75 m/s','100 m/s', '125 m/s');
% name={'Total';'URotor';'LRotor';'Hinge';'Spring';'Fuselage';'HStab + Elev'};
name={'Total';'Rotor';'Fuselage';'HStab'};
set(gca,'xticklabel', name);
% xlabel('Component Contribution to $M_w$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$M_w$', 'Interpreter', 'latex', 'FontSize', 14);
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% figure(4)
% hp = bar(bary);
% 
% hatchfill2(hp(1),'single','HatchAngle',0,'hatchcolor','k');
% hatchfill2(hp(2),'cross','HatchAngle',45,'hatchcolor','k');
% hatchfill2(hp(3),'single','HatchAngle',45,'hatchcolor','k');
% hatchfill2(hp(4),'single','HatchAngle',-45,'hatchcolor','k');
% hatchfill2(hp(5),'cross','HatchAngle',30,'hatchcolor','k');
% hatchfill2(hp(6),'cross','HatchAngle',-30,'hatchcolor','k');
% 
% for k = 1:numel(hp)
%     hp(k).FaceColor = 'none';
% end
% 
% legend('0 kts', '50 kts', '100 kts', '150 kts','200 kts', '250 kts');
% name={'Total';'URotor';'LRotor';'Hinge';'Flap';'Fuselage';'HStab'};
% set(gca,'xticklabel',name);
% xlabel('Component Contribution to M_w'); ylabel('M_w');
% grid on;

%% NATURAL MODES OF MOTION
% Longitudinal Modes
%%% Phugoid
% Lynx Values for Testing
% X_u = -0.02;
% M_u = 0.047;
% M_q = -1.9;

a_phug = 1;
b_phug = (X_u + g * M_u ./ (M_q.^2));
c_phug = -g.*M_u./M_q;
discriminant_phug = b_phug.^2 - 4*a_phug*c_phug;
lambda_phug_1 = (-b_phug + sqrt(discriminant_phug)) / (2*a_phug);
lambda_phug_2 = (-b_phug - sqrt(discriminant_phug)) / (2*a_phug);


%%% Short Period
U_e = V_vals;
a_sp = 1;
b_sp = -(Z_w+M_q);
c_sp = Z_w.*M_q - (U_e + Z_q).*M_w;
discriminant_sp = b_sp.^2 - 4*a_sp*c_sp;
lambda_sp_1 = (-b_sp + sqrt(discriminant_sp)) / (2*a_sp);
lambda_sp_2 = (-b_sp - sqrt(discriminant_sp)) / (2*a_sp);


% PLOTTING NATURAL MODES
% Generate colors based on the entry index
num_points = numel(lambda_phug_1);
colors = jet(num_points);  % Use jet colormap for varying colors

figure(23);
    subplot(2,2,1)
    scatter(real(lambda_phug_1), imag(lambda_phug_1), [], colors, 's', 'filled', 'DisplayName', 'phug');
    hold on;
    scatter(real(lambda_phug_2), imag(lambda_phug_2), [], colors, 's', 'filled');
    % load('phugHQ1.mat'); plot(phugHQ1(:,1), phugHQ1(:,2), '-.k');
    % load('phugHQ2.mat'); plot(phugHQ2(:,1), phugHQ2(:,2), '-.k');
    % xlabel('$\mu$ [1/sec]', 'Interpreter', 'latex'); 
    ylabel('$\omega$ [rad/sec]', 'Interpreter', 'latex', 'FontSize', 12); 
    title('Phugoid', 'Interpreter', 'latex', 'FontSize', 12); grid on;
    xlim([-0.2, 0.2]); %ylim([0, max(phugHQ2(:,2))]);
    hold off;

    subplot(2,2,2)
    scatter(real(lambda_sp_1), imag(lambda_sp_1), [], colors, 's', 'filled', 'DisplayName', 'sp');
    hold on;
    scatter(real(lambda_sp_2), imag(lambda_sp_2), [], colors, 's', 'filled');
    % colorbar(colors)
    % xlabel('$\mu$ [1/sec]', 'Interpreter', 'latex'); 
    % ylabel('$\omega$ [rad/sec]', 'Interpreter', 'latex'); 
    title('Short Period', 'Interpreter', 'latex', 'FontSize', 12); grid on;hold off;

    subplot(2,2,3)
    scatter(Z_w, 0, [], colors, 's', 'filled', 'DisplayName', 'Mq');
    xlabel('$\mu$ [1/sec]', 'Interpreter', 'latex', 'FontSize', 12); 
    % ylabel('$\omega$ [rad/sec]', 'Interpreter', 'latex'); 
    title('Heave Subsidence', 'Interpreter', 'latex', 'FontSize', 12); grid on;
    
    subplot(2,2,4)
    scatter(M_q, 0, [], colors, 's', 'filled', 'DisplayName', 'Mq');
    xlabel('$\mu$ [1/sec]', 'Interpreter', 'latex', 'FontSize', 12); 
    ylabel('$\omega$ [rad/sec]', 'Interpreter', 'latex', 'FontSize', 12); 
    title('Pitch Subsidence', 'Interpreter', 'latex', 'FontSize', 12); grid on;

    hcb = colorbar; colormap(jet); clim([0 120]); 
    set(hcb, 'Position', [0.92, 0.1, 0.02, 0.8]);


    

% figure()
% scatter(real(lambda_phug_1), imag(lambda_phug_1), [], colors, 's', 'filled', 'DisplayName', 'sp'); hold on; 
% scatter(real(lambda_phug_2), imag(lambda_phug_2), [], colors, 's', 'filled', 'DisplayName', 'sp'); hold off; grid on;

% figure()
% scatter(real(lambda_sp_1), imag(lambda_sp_1), [], colors, 's', 'filled', 'DisplayName', 'sp'); hold on; 
% scatter(real(lambda_sp_2), imag(lambda_sp_2), [], colors, 's', 'filled', 'DisplayName', 'sp'); hold off; grid on;

springterm = Z_w.*M_q - M_w.*(Z_q+V_vals);
figure(1)
plot(V_vals, springterm, '-*', 'LineWidth', 2, 'Color', '#D95319'); grid on;
xlabel('Trim speed $U_e$ [m/s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Spring Term Equation Result', 'Interpreter', 'latex', 'FontSize', 12);


indices = [1, 6, 11, 16, 21, 26];
%heave subsidence
eig1 = Z_w(indices)
% pitch subsidence
eig2 = M_q(indices)
% phugoid
eig3 = lambda_phug_1(indices)
% dutch roll
eig4 = lambda_sp_1(indices)



%% CONTROL ALLOCATION

% epsilon = 0.001;
% 
% for i=1:length(V)
%     V_kts(i) = V(i)*1.944;
% 
%     if V_kts(i)<V_l_transit
%         rho_theta_s(i) = 1;
%         rho_delta_e(i) = epsilon;
% 
%     elseif V_kts(i)>=V_l_transit && V_kts(i)<V_u_transit
%         rho_theta_s(i) = 1 - 1/(V_u_transit-V_l_transit)*(V_kts(i)-V_l_transit);
%         rho_delta_e(i) = 1/(V_u_transit-V_l_transit)*(V_kts(i)-V_l_transit);
% 
%     else
%         rho_theta_s(i) = epsilon;
%         rho_delta_e(i) = 1;
%     end
% end
% 
% rho = [rho_theta_s; rho_delta_e];
% 
% % PLOT CONTROL ALLOCATION DISTRIBUTION FOR PITCH
% figure(41)
% bar(V_vals_kts(1:2:12), rho(:,1:2:12), 'stacked');
% xlabel('Forward velocity $u$ [kts]', 'Interpreter', 'latex', 'FontSize', 14); 
% ylabel('Control Weight', 'Interpreter', 'latex', 'FontSize', 14);
% legend('$\theta_s$', '$\delta_e$', 'Interpreter', 'latex', 'FontSize', 14);
% ylim([0,1]);
% % title('Control allocation weights for pitch');
% grid on;
% figure(42)
% plot(V_vals_kts(1:12), rho_theta_s(1:12), V_vals_kts(1:12), rho_delta_e(1:12), '-.', 'LineWidth', 2); grid on; 
% ylabel('Control Weight'); xlabel('Forward velocity u [kts]');
% title('Control allocation weights for pitch');
% legend('\theta_s', '\delta_e')

%% Handling Quality Analysis
% figure(61)
% pitch_tf_fit = tf(4^2, [1 2*0.707*4 4^2]);
% 
% TrimPoint = 8;
% sys_lin = ss(TrimLinSave_6dof{TrimPoint,1},TrimLinSave_6dof{TrimPoint,2}, ...
%     TrimLinSave_6dof{TrimPoint,3},TrimLinSave_6dof{TrimPoint,4});
% 
% bode(sys_lin(5, 3)); hold on;
% bode(sys_lin(5, 7));
% % bode(sys_lin{1, 2}(3,2)); bode(sys_lin{1, 3}(3,2));
% % [GM, PM, w_gm, w_pm] = margin(sys_lin(5, 3));
% 
% 
% bode(pitch_tf_fit , 'r--'); bode(tf(-15, [0.30 1]), 'k--');
% hold off; grid on;
% legend('linSys', 'command model', 'inverse model')
% 
% [mag, phase, wout] = bode(sys_lin(5, 3));
% 
% % Convert magnitude to dB
% magdB = 20*log10(mag);
% 
% % Convert phase to degrees
% phase_deg = squeeze(phase); % Extract phase data
% phase_deg = phase_deg(:); % Reshape to a column vector
% phase_deg = phase_deg(1:length(wout)); % Extract valid part of phase data
% 
% % Find the frequency where phase is closest to 45 degrees
% desiredPhase = 180-45; % desired phase in degrees
% [~, idx] = min(abs(phase_deg - desiredPhase));
% 
% % Extract corresponding frequency
% w_bw_phase = wout(idx);
% 
% % Display the result
% disp(['Frequency for phase of 135 degrees: ', num2str(w_bw_phase), ' rad/s']);

%% Cross-Coupling
% Select Trim Speed
% TrimPoint = 8;
% sys_lin = ss(TrimLinSave_6dof{TrimPoint,1},TrimLinSave_6dof{TrimPoint,2}, ...
%     TrimLinSave_6dof{TrimPoint,3},TrimLinSave_6dof{TrimPoint,4});
% 
% % Pitch due to roll q/delta_lat --> sys_lin(5,4)
% bode(sys_lin(5, 4)); hold on; grid on;
% 
% % Roll due to pitch p/delta_lon --> sys_lin(5,3), sys_lin(5,7)
% bode(sys_lin(5, 3));
% bode(sys_lin(5, 7));
% hold off;
% legend('q/p', 'p/q LonCyc', 'p/q Elev')


%%
%% Cross-Coupling
% Select Trim Speed
% TrimPoint = 8;
% sys_lin = ss(TrimLinSave_6dof{TrimPoint,1}, TrimLinSave_6dof{TrimPoint,2}, ...
%     TrimLinSave_6dof{TrimPoint,3}, TrimLinSave_6dof{TrimPoint,4});
% 
% % Define frequency range for the Bode plot
% w = logspace(0, 2, 1000); % Frequency range from 1 to 100 rad/s
% 
% % Plot Bode plots and hold on to add more plots
% figure;
% [mag_qp, phase_qp, wout_qp] = bode(sys_lin(5, 4), w); 
% [mag_pq_lon, phase_pq_lon, wout_pq_lon] = bode(sys_lin(5, 3), w); 
% [mag_pq_elev, phase_pq_elev, wout_pq_elev] = bode(sys_lin(5, 7), w);
% 
% % Convert magnitude to dB
% mag_qp_db = 20*log10(squeeze(mag_qp));
% mag_pq_lon_db = 20*log10(squeeze(mag_pq_lon));
% mag_pq_elev_db = 20*log10(squeeze(mag_pq_elev));
% 
% % Plot the Bode magnitude plots
% subplot(2,1,1);
% semilogx(wout_qp, mag_qp_db, 'b'); hold on;
% semilogx(wout_pq_lon, mag_pq_lon_db, 'r');
% semilogx(wout_pq_elev, mag_pq_elev_db, 'g');
% grid on;
% ylabel('Magnitude (dB)');
% legend('q/\delta_{lat}', 'p/\delta_{lon}', 'p/\delta_{elev}');
% title('Bode Plot');
% 
% % Add vertical lines at 1 Hz and 10 Hz
% xline(1, '--k', '1 Hz');
% xline(10, '--k', '10 Hz');
% 
% % Plot the Bode phase plots
% subplot(2,1,2);
% semilogx(wout_qp, squeeze(phase_qp), 'b'); hold on;
% semilogx(wout_pq_lon, squeeze(phase_pq_lon), 'r');
% semilogx(wout_pq_elev, squeeze(phase_pq_elev), 'g');
% grid on;
% ylabel('Phase (deg)');
% xlabel('Frequency (rad/s)');
% legend('q/\delta_{lat}', 'p/\delta_{lon}', 'p/\delta_{elev}');
% 
% % Add vertical lines at 1 Hz and 10 Hz
% xline(1, '--k', '1 Hz');
% xline(10, '--k', '10 Hz');
% 
% % Calculate the average gain value between 1 Hz and 10 Hz
% f1 = 1; % 1 Hz
% f2 = 10; % 10 Hz
% 
% % Find the indices corresponding to 1 Hz and 10 Hz
% index_f1_qp = find(wout_qp >= f1, 1, 'first');
% index_f2_qp = find(wout_qp <= f2, 1, 'last');
% 
% index_f1_pq_lon = find(wout_pq_lon >= f1, 1, 'first');
% index_f2_pq_lon = find(wout_pq_lon <= f2, 1, 'last');
% 
% index_f1_pq_elev = find(wout_pq_elev >= f1, 1, 'first');
% index_f2_pq_elev = find(wout_pq_elev <= f2, 1, 'last');
% 
% % Calculate the average magnitude in dB
% avg_gain_qp = mean(mag_qp_db(index_f1_qp:index_f2_qp));
% avg_gain_pq_lon = mean(mag_pq_lon_db(index_f1_pq_lon:index_f2_pq_lon));
% avg_gain_pq_elev = mean(mag_pq_elev_db(index_f1_pq_elev:index_f2_pq_elev));
% 
% % Display the average gain values
% disp(['Average gain of q/\delta_{lat} between 1 Hz and 10 Hz: ', num2str(avg_gain_qp), ' dB']);
% disp(['Average gain of p/\delta_{lon} between 1 Hz and 10 Hz: ', num2str(avg_gain_pq_lon), ' dB']);
% disp(['Average gain of p/\delta_{elev} between 1 Hz and 10 Hz: ', num2str(avg_gain_pq_elev), ' dB']);


%% Fan Plot
% RPM = 0:1:400;
% Omega = RPM.*(2*pi/60);
% Omega_rotor_nominal = 35;
% normalized_rotorspeed = Omega./Omega_rotor_nominal;
% 
% % hinge, no spring
% e = 0.45515;       % TUNE THIS PARAMETER FOR CORRECT flapping frequency 
%                 % BERGER HAS 1.49/rev, NLR 1.5/rev
% Southwell = 1+3*e/(2*(1-e));
% omega_nr2 = 0;
% flapfreq = omega_nr2 + Southwell*Omega.^2;    %in rad/s
% flapfreq_norm = sqrt(flapfreq)./Omega_rotor_nominal;
% 
% % hinge and spring
% e_spring_ratio = 0.04;
% e_spring = e_spring_ratio*R;
% K_b = 900000;       % Tunes the starting location of the plot
% I_b = 1/3*m_bl*R^2;
% M_b = -500;
% Southwell = 1+3*e_spring/(2*(1-e_spring));
% omega_nr2_spring = K_b/I_b;
% % flapfreq_spring = omega_nr2_spring + Southwell*Omega.^2;    %in rad/s
% flapfreq_spring = K_b/I_b + (1+e_spring*M_b/I_b)*Omega.^2;
% flapfreq_norm_spring = sqrt(flapfreq_spring)./Omega_rotor_nominal;
% 
% figure(81);
% plot(normalized_rotorspeed, flapfreq_norm, 'LineWidth', 2);
% hold on;
% plot(normalized_rotorspeed, flapfreq_norm_spring, 'LineWidth', 2)
% 
% load('flapFreqBerger')
% flapFreqBerger = sgolayfilt(flapFreqBerger, 2, 25); % Polynomial order 2, frame length 11
% plot(flapFreqBerger(:,1), flapFreqBerger(:,2), 'r--', 'LineWidth', 1.5)
% 
% maxperrev = 5;
% revlineyval = zeros(maxperrev, size(normalized_rotorspeed, 2));
% 
% for i = 1:maxperrev
%     revlineyval(i,:) = i * normalized_rotorspeed;
% end

plot(normalized_rotorspeed, revlineyval, 'k-.');
grid on;

% Add text labels
for i = 1:maxperrev
    text(normalized_rotorspeed(end), revlineyval(i, end), ['\bf', num2str(i), '/rev   '], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

[~, idx_closest_to_1] = min(abs(normalized_rotorspeed - 1));
flapmodefreq = flapfreq_norm(idx_closest_to_1);
disp(['Flap mode Frequency: ', num2str(flapmodefreq), '/rev']);

xlabel('Normalized Rotor Speed (\Omega/\Omega_0)');
ylabel('Normalized Frequency (\omega/\Omega_0)'); ylim([0,3]);
yline(flapmodefreq, '--');

legend('Hinge Only', 'Hinge and Spring', '1st Flapping Berger', 'Location', 'northwest')

hold off;




