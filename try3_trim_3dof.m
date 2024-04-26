%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename:     coaxial_3dof_trim
% Project:      MSc Thesis
% Supervisor:   M.D. Pavel
% Author:       Ynias Prencipe 
% Student Nr.:  4777158
% 
% Description:  Trim using the Newton-Rhapson algorithm for the coaxial 
% helicopter, no pusher prop or elevator added yet
% 
% Assumptions:  - Uniform inflows on both upper and lower rotors
%               - Upper rotor inflow is fully added to lower rotor inflow
%               - No rotor interference from lower to upper
%               - No wake contraction
% 
% Issues to be solved:  - CT lower drops with speed , and becomes negative
%                       - Strange behaviour of inflow at low speeds
%                       - BEM thrust coefficient (CT) in wrong plane?
% 
% To Add:       - Pusher propeller
%               - Elevator + Horizontal Stabilizer / Stabilator
%               ? Wake contraction effects
%               ? Rotor interference factor in function of advance ratio
%               - take into account the flight envelope in which this works
%               , look at momentum theory definitions for the regions
%               (especially descent etc.] for which its not valid anymore
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Set-up
% Load coaxial helicopter parameters
coaxial_heli_parameters;

% Initial Values
theta_0(1) = deg2rad(0);
theta_c(1) = deg2rad(0);
theta_f(1)=deg2rad(0);
theta_p(1) = deg2rad(0);
delta_e(1) = deg2rad(0);
t(1)=0;
u(1)=0;
w(1)=0;
q(1)=0;
Omega_init = 40;
lambda_0_u(1)=sqrt(mass/2*abs(g)/(area*2*rho))/(Omega_init*R);
lambda_0_l(1)=sqrt(mass/2*abs(g)/(area*2*rho))/(Omega_init*R);
lambda_p(1) = 0.002;
vi_hov = sqrt(1/2*mass*g/(2*rho*pi*R^2));

% Define trim variable vector
x_k(:,1) = [theta_0(1); theta_c(1); theta_p(1); lambda_0_u(1); lambda_0_l(1); lambda_p(1)];


%% Trim Routine

for i = 1:length(V_vals)
    i
    V(i) = V_vals(i);
    u(i) = V(i)*cos(atan2(w(i),u(i)));
    v(i) = 0;
    w(i) = V(i)*sin(atan2(w(i),u(i)));
    p(i) = 0;
    r(i) = 0;
    % theta_f(i) = theta_f_trim(i);
    theta_f(i) = 0;
    delta_e(i) = 0;

    %% Rotor RPM Scheduling
    if V(i) < 70
        Omega = 40;
    else
        % Omega = 40 - 5/30*(V(i)-70);
        Omega = 40;
    end

    %%
    if i>1
        x_k(:,i) = x_k(:,i-1);      % initial guess of each trim iteration is the previous trim state - for quick convergence
    end

    %% Define State Variable Names
    theta_0_u(i) = x_k(1,i)/2; theta_c(i) = x_k(2,i); 
    theta_p(i) = x_k(3,i); lambda_0_u(i) = x_k(4,i); lambda_0_l(i) = x_k(5,i);
    lambda_0_p = x_k(6,i);
   
    %%
    % Calculate Angles
    alfa_sp(i) = atan2(w(i),u(i));
    alfa_cp(i) = alfa_sp(i) + x_k(2,i);

    % Speed variables
    mu(i) = V(i)/(Omega*R);
    mu_x(i) = V(i)/(Omega*R) * cos(alfa_cp(i));
    mu_z(i) = V(i)/(Omega*R) * sin(alfa_cp(i));
    vel(:,i) = [u(i); w(i); mu(i); mu_x(i); mu_z(i)];

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
        f(:,i) = f_xk(vel(:,i), x_ki(:,i),  q(i), theta_f(i), delta_e(i));

        % calculate x_k+1 for perturbation in x(1,i)
        ti1 = [tistep; 0; 0; 0; 0; 0];
        x_k1(:,i) = x_ki(:,i) + ti1;
        % calculate f_xk1_ti(:,i) using x_k+1 for FIRST column of jacobian
        f_xk1_ti(:,i) = f_xk(vel(:,i), x_k1(:,i), q(i), theta_f(i), delta_e(i));

        % calculate x_k+1 for perturbation in x(2,i)
        ti2 = [0; tistep; 0; 0; 0; 0];
        x_k2(:,i) = x_ki(:,i) + ti2;
        % calculate f_xk2_ti(:,i) using x_k+1 for SECOND column of jacobian
        f_xk2_ti(:,i) = f_xk(vel(:,i), x_k2(:,i), q(i), theta_f(i), delta_e(i));
        
        % calculate x_k+1 for perturbation in x(3,i)
        ti3 = [0; 0; tistep; 0; 0; 0];
        x_k3(:,i) = x_ki(:,i) + ti3;
        % calculate f_xk_ti(:,i) using x_k+1 for THIRD column of jacobian
        f_xk3_ti(:,i) = f_xk(vel(:,i), x_k3(:,i), q(i), theta_f(i), delta_e(i));

        % calculate x_k+1 for perturbation in x(4,i)
        ti4 = [0; 0; 0; tistep; 0; 0];
        x_k4(:,i) = x_ki(:,i) + ti4;
        % calculate f_xk_ti(:,i) using x_k+1 for THIRD column of jacobian
        f_xk4_ti(:,i) = f_xk(vel(:,i), x_k4(:,i), q(i), theta_f(i), delta_e(i));

        % calculate x_k+1 for perturbation in x(5,i)
        ti5 = [0; 0; 0; 0; tistep; 0];
        x_k5(:,i) = x_ki(:,i) + ti5;
        % calculate f_xk_ti(:,i) using x_k+1 for THIRD column of jacobian
        f_xk5_ti(:,i) = f_xk(vel(:,i), x_k5(:,i), q(i), theta_f(i), delta_e(i));

        % calculate x_k+1 for perturbation in x(5,i)
        ti6 = [0; 0; 0; 0; 0; tistep];
        x_k6(:,i) = x_ki(:,i) + ti6;
        % calculate f_xk_ti(:,i) using x_k+1 for THIRD column of jacobian
        f_xk6_ti(:,i) = f_xk(vel(:,i), x_k6(:,i), q(i), theta_f(i), delta_e(i));
        
        % assemble all 3 columns for Jacobian
        ti = tistep * ones(size(f(:,i)));
        Jac = [(f_xk1_ti(:,i)-f(:,i))./ti, (f_xk2_ti(:,i)-f(:,i))./ti, (f_xk3_ti(:,i)-f(:,i))./ti, ...
            (f_xk4_ti(:,i)-f(:,i))./ti, (f_xk5_ti(:,i)-f(:,i))./ti, (f_xk6_ti(:,i)-f(:,i))./ti];
        
        % calculate new input vector according to Jacobian
        x_k(:,i) = x_ki(:,i) - inv(Jac) * f(:,i);

        % f(:,i)
        [f(:,i), output_other(:,i), Mw_components(:,i), angles(:,i), inflow(:,i)] = f_xk(vel(:,i), x_k(:,i), q(i), theta_f(i), delta_e(i));

        f_sum(count(i)) = sum(abs(f(:,i)),"all");
        
    end

    udot(i) = f(1,i);
    wdot(i) = f(2,i);
    qdot(i) = f(3,i);
    lambda_0_u_dot(i) = f(4,i);

    u(i+1) = V(i)*cos(atan2(w(i),u(i)));
    w(i+1) = V(i)*sin(atan2(w(i),u(i)));
    q(i+1) = 0;
    

end

%% SAVE x_k(:,i
if exist('trim_saved.mat', 'file') == 2
    load('trim_saved.mat', 'trim_saved');
    % Append the new x_k to the existing vector
    trim_saved = x_k;
else
    % Create a new vector to store x_k
    trim_saved = x_k;
end

% Save the updated x_k to a file
save('trim_saved.mat', 'trim_saved');

%% TRIM PLOT

% figure(10)
% theta_0 = rad2deg(x_k(1,:)); 
% theta_c = rad2deg(x_k(2,:)); 
% theta_p = rad2deg(x_k(3,:));
% subplot(5,1,1)
% plot(mu, theta_0, '-*'), grid;
% legend('\theta_0'); title('Mean Collective')
% subplot(5,1,2)
% plot(mu, theta_c, '-*'), grid;
% legend('\theta_c'); title('Longitudinal Cyclic')
% subplot(5,1,3)
% plot(mu, theta_p, '-*'), grid;
% legend('\theta p'); title('Pusher Prop Collective')
% subplot(5,1,4)
% plot(mu, rad2deg(delta_e), '-*'), grid;
% legend('\delta e'); title('Elevator Deflection (Scheduled)')
% subplot(5,1,5)
% plot(mu, rad2deg(theta_f), '-*'), grid;
% legend('\theta f'); xlabel("Advance Ratio \mu [-]"); title('Fuselage Pitch (Scheduled)')

%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ROTOR ANGLES PLOT

% figure(11)
% V_vals_kts = V_vals.*1.944;
% a0_u = angles(1,:); a0_l = angles(2,:); a1_u = angles(3,:); 
% a1_l = angles(4,:); a1r_u = angles(5,:); a1r_l = angles(6,:); 
% alfa_sp = angles(7,:); alfa_cp = angles(8,:); alfa_dp_u = angles(9,:);
% alfa_dp_l = angles(10,:);
% plot(V_vals_kts, a0_u, V_vals_kts, a0_l, ...
%     V_vals_kts, a1_u, V_vals_kts, a1_l, ...
%     V_vals_kts, a1r_u, 'k--', V_vals_kts, a1r_l, 'k--', ...
%     V_vals_kts, alfa_sp, 'g--o', V_vals_kts, alfa_cp, 'b--o');
% grid on;
% legend('a0_u', 'a0_l', 'a1_u', 'a1_l', 'a1r_u', 'a1r_l', 'alfa_{sp}', 'alfa_{cp}');

%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Some inflow trials
lambda_0_u = inflow(1,:); lambda_0_l = inflow(2,:); 

v_u = lambda_0_u*Omega*R;
v_l = lambda_0_l*Omega*R;

delta_l2u = zeros(1,length(mu));
delta_u2l = zeros(1,length(mu));

for i = 1:length(mu)
    if (0<=mu(i)) && (mu(i)<= 0.316)
        delta_l2u(i) = -2.15*mu(i) + 0.68;
    else
        % delta_l2u = -2.15*mu(i) + 0.68;     % this is not correct but 0 gives issues
        delta_l2u(i) = 0;
    end
    if (0<=mu(i)) && (mu(i)<=0.381)
        delta_u2l(i) = -3.81*mu(i) + 1.45;    
    else
        % delta_u2l = -3.81*mu(i) + 1.45;      % not correct but 0 gives issues
        delta_u2l(i) = 0;                         % 0 gives jumps in C_T and a_1
    end
end

v_i_u = v_u + delta_l2u.*v_l; % + (K_u*v_u+K_l*delta_u*v_l)*r/R*cos(psi_u)
v_i_l = v_l + delta_u2l.*v_u; 

% figure(12)
% plot(V_vals_kts, v_u, V_vals_kts, v_l, ...
%     V_vals_kts, v_i_u, '--', V_vals_kts, v_i_l, '--');
% grid on;
% legend('lambda_0_u', 'lambda_0_l');

% %%
% C_T_Glau_u = 2*v_i_u/(Omega*R) .* sqrt( (mu.*cos(alfa_dp_u)).^2 + (mu.*sin(alfa_dp_u) + v_i_u/(Omega*R)).^2 );
% C_T_Glau_l = 2*v_i_l/(Omega*R) .* sqrt( (mu.*cos(alfa_dp_l)).^2 + (mu.*sin(alfa_dp_l) + v_i_l/(Omega*R)).^2 );
% C_T_BEM_u = sigma*c_l_a/2*( (1/3+mu_x.^2/2).*theta_0 + (1+mu_x.^2)/8*twist_r + lambda_u/2 );
% C_T_BEM_l = sigma*c_l_a/2*( (1/3+mu_x.^2/2).*theta_0 + (1+mu_x.^2)/8*twist_r + lambda_l/2 );
% 
% figure(13)
% plot(V_vals_kts, C_T_Glau_u, V_vals_kts, C_T_Glau_l);
% grid on;
% legend('C_T Glau U', 'C_T Glau L');












%% Linearization 
lin_step = 0.1;
lin_systems = cell(length(V_vals), 2);
% at V = 0.1m/s
for i = 1:length(V_vals)
    i
    trim_vals_0 = x_k(:,i);

    %%%%%%%%%%%%%%%%%%% A-MATRIX %%%%%%%%%%%%%%%%%%%
    % Row 1 (X-force --> udot state)
    vel_lin0 = [lin_step; 0; 0; 0; 0];
    f_lin_dist_0_pos = f_xk(vel(:,i)+vel_lin0, trim_vals_0, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk(vel(:,i)-vel_lin0, trim_vals_0, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_u(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; lin_step; 0; 0; 0];
    f_lin_dist_0_pos = f_xk(vel(:,i)+vel_lin0, trim_vals_0, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk(vel(:,i)-vel_lin0, trim_vals_0, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_w(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0, q(i)+0.01, theta_f(i), delta_e(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0, q(i)-0.01, theta_f(i), delta_e(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_q(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    % Row 2 (Z-force --> wdot state)
    vel_lin0 = [lin_step; 0; 0; 0; 0];
    f_lin_dist_0_pos = f_xk(vel(:,i)+vel_lin0, trim_vals_0, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk(vel(:,i)-vel_lin0, trim_vals_0, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Z_u(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; lin_step; 0; 0; 0];
    f_lin_dist_0_pos = f_xk(vel(:,i)+vel_lin0, trim_vals_0, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk(vel(:,i)-vel_lin0, trim_vals_0, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Z_w(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    qdiff = 0.01;
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0, q(i) + qdiff, theta_f(i), delta_e(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0, q(i) - qdiff, theta_f(i), delta_e(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Z_q(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*qdiff);

    % Row 3 (M-moment --> qdot state)
     
    vel_lin0 = [lin_step; 0; 0; 0; 0];
    f_lin_dist_0_pos = f_xk(vel(:,i)+vel_lin0, trim_vals_0, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk(vel(:,i)-vel_lin0, trim_vals_0, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    M_u(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; lin_step; 0; 0; 0];
    [f_lin_dist_0_pos, M_moments_Mw_pos, Mw_components_pos] = ...
        f_xk(vel(:,i)+vel_lin0, trim_vals_0, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_pos_Mw = f_lin_dist_0_pos(3,1);
    [f_lin_dist_0_neg, M_moments_Mw_neg, Mw_components_neg] = ...
        f_xk(vel(:,i)-vel_lin0, trim_vals_0, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_neg_Mw = f_lin_dist_0_neg(3,1);
    M_w(i) = (f_lin_dist_0_pos_Mw - f_lin_dist_0_neg_Mw)/(2*lin_step);
    M_w_components(:,i) = (Mw_components_pos - Mw_components_neg)./(2*lin_step);

    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0, q(i)+0.01, theta_f(i), delta_e(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0, q(i)-0.01, theta_f(i), delta_e(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    M_q(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    % Assemble A-matrix
    A_lin = [X_u(i), X_w(i), X_q(i);
            Z_u(i), Z_w(i), Z_q(i);
            M_u(i), M_w(i), M_q(i)]

    %%%%%%%%%%%%%%%%%%% B-MATRIX %%%%%%%%%%%%%%%%%%%
    lin_step = deg2rad(0.01);

    % Row 1 (X-force --> qdot state)
    b_step = [lin_step; 0; 0; 0; 0; 0];   % input in collective
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_theta_0(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; lin_step; 0; 0; 0; 0];   % input in forward cyclic
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_theta_c(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; lin_step; 0; 0; 0];
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_theta_p(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    delta_e_diff = lin_step;
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0, q(i), theta_f(i), delta_e(i)+delta_e_diff);
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0, q(i), theta_f(i), delta_e(i)-delta_e_diff);
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_delta_e(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    % Row 2 (Z-force --> qdot state)
    b_step = [lin_step; 0; 0; 0; 0; 0];   % input in collective
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Z_theta_0(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; lin_step; 0; 0; 0; 0];   % input in forward cyclic
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Z_theta_c(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; lin_step; 0; 0; 0];
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Z_theta_p(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    delta_e_diff = lin_step;
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0, q(i), theta_f(i), delta_e(i)+delta_e_diff);
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0, q(i), theta_f(i), delta_e(i)-delta_e_diff);
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Z_delta_e(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    % Row 3 (M-moment --> qdot state)
    b_step = [lin_step; 0; 0; 0; 0; 0];   % input in collective
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    M_theta_0(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; lin_step; 0; 0; 0; 0];   % input in forward cyclic
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    M_theta_c(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; lin_step; 0; 0; 0];
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i), theta_f(i), delta_e(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    M_theta_p(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    delta_e_diff = lin_step;
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0, q(i), theta_f(i), delta_e(i)+delta_e_diff);
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0, q(i), theta_f(i), delta_e(i)-delta_e_diff);
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    M_delta_e(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    % Assemble B-MATRIX
    B_lin = [X_theta_0(i), X_theta_c(i), X_theta_p(i), X_delta_e(i);
            Z_theta_0(i), Z_theta_c(i), Z_theta_p(i), Z_delta_e(i);
            M_theta_0(i), M_theta_c(i), M_theta_p(i), M_delta_e(i)]

    % Assemble C-Matrix
    C_lin = eye(size(A_lin));

    % Assemble D-Matrix
    D_lin = zeros(size(B_lin));       % #columns = #inputs

    % Analysis with state-space
    sys_lin{i} = ss(A_lin, B_lin, C_lin, D_lin);

    % pzplot(sys_lin)

    % Save Current Trim SS System in lin_systems
    TrimLin{i, 1} = A_lin;
    TrimLin{i, 2} = B_lin;
    TrimLin{i, 3} = C_lin;
    TrimLin{i, 4} = D_lin;
end

% Save Linear Systems around Trim states
if exist('TrimLinSave.mat', 'file') == 2
    load('TrimLinSave.mat', 'TrimLinSave');
    % Append the new x_k to the existing vector
    TrimLinSave = TrimLin;
else
    % Create a new vector to store x_k
    TrimLinSave = TrimLin;
end

% Save the updated x_k to a file
save('TrimLinSave.mat', 'TrimLinSave');


V_vals_kts = V_vals.*1.944;

%% STABILITY DERIVATIVES PLOTS

figure(10)
subplot(3,3,1)
plot(V_vals_kts, X_u, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_u'); legend('X_u')
subplot(3,3,2)
plot(V_vals_kts, X_w, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_w'); legend('X_w')
subplot(3,3,3)
plot(V_vals_kts, X_q, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_q'); legend('X_q')
subplot(3,3,4)
plot(V_vals_kts, Z_u, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_u'); legend('Z_u')
subplot(3,3,5)
plot(V_vals_kts, Z_w, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_w'); legend('Z_w')
subplot(3,3,6)
plot(V_vals_kts, Z_q, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_q'); legend('Z_q')
subplot(3,3,7)
plot(V_vals_kts, M_u, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_u'); legend('M_u')
subplot(3,3,8)
plot(V_vals_kts, M_w, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_w'); legend('M_w')
subplot(3,3,9)
plot(V_vals_kts, M_q, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_q'); legend('M_q')

%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CONTROL DERIVATIVES PLOTS

figure(11)
subplot(3,4,1)
plot(V_vals_kts, X_theta_0, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_theta_0'); legend('X_theta_0')
subplot(3,4,2)
plot(V_vals_kts, X_theta_c, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_theta_c'); legend('X_theta_c')
subplot(3,4,3)
plot(V_vals_kts, X_theta_p, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_theta_p'); legend('X_theta_p')
subplot(3,4,4)
plot(V_vals_kts, X_delta_e, '--*'); grid on; xlabel('V_f [kts]'); ylabel('X_delta_e'); legend('X_delta_e')
subplot(3,4,5)
plot(V_vals_kts, Z_theta_0, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_theta_0'); legend('Z_theta_0')
subplot(3,4,6)
plot(V_vals_kts, Z_theta_c, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_theta_c'); legend('Z_theta_c')
subplot(3,4,7)
plot(V_vals_kts, Z_theta_p, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_theta_p'); legend('Z_theta_p')
subplot(3,4,8)
plot(V_vals_kts, Z_delta_e, '--*'); grid on; xlabel('V_f [kts]'); ylabel('Z_delta_e'); legend('Z_delta_e')
subplot(3,4,9)
plot(V_vals_kts, M_theta_0, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_theta_0'); legend('M_theta_0')
subplot(3,4,10)
plot(V_vals_kts, M_theta_c, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_theta_c'); legend('M_theta_c')
subplot(3,4,11)
plot(V_vals_kts, M_theta_p, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_theta_p'); legend('M_theta_p')
subplot(3,4,12)
plot(V_vals_kts, M_delta_e, '--*'); grid on; xlabel('V_f [kts]'); ylabel('M_delta_e'); legend('M_delta_e')

%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CONTROL EFFECTIVENESS PLOTS
% Find the index where V_vals_kts is closest to 10
transition_min_val = 20;        % [m/s]
[~, idx_transstart] = min(abs(V - transition_min_val));

% Find the intersection point between M_theta_c and M_delta_e
[~, idx_intersection] = min(abs(abs(M_theta_c) - abs(M_delta_e)));

% Create the plot
figure(91);
plot(V_vals_kts, abs(M_theta_c), V_vals_kts, abs(M_delta_e), '--*');
grid on;
xlabel('V_f [kts]');
ylabel('CE_\theta (1/s^2)');
hold on;
% Plot vertical lines at x = 40 and x = 121.6

V_l_transit = V_vals_kts(idx_transstart);
V_u_transit = V_vals_kts(idx_intersection);
V_transit = [V_l_transit, V_u_transit];

% Save Transition Speeds
if exist('V_transit_save.mat', 'file') == 2
    load('V_transit_save.mat', 'V_transit_save');
    % Append the new x_k to the existing vector
    V_transit_save = V_transit;
else
    % Create a new vector to store x_k
    V_transit_save = V_transit;
end

% Save the updated x_k to a file
save('V_transit_save.mat', 'V_transit_save');

line([V_l_transit, V_l_transit], ylim, 'Color', 'r', 'LineStyle', '--');
line([V_u_transit, V_u_transit], ylim, 'Color', 'r', 'LineStyle', '--');
% Shade the region between the two lines
x_shade = [V_vals_kts(idx_transstart), V_vals_kts(idx_intersection)];
y_shade = ylim; % Use ylim to define y-values for the shaded region
fill([x_shade, fliplr(x_shade)], [y_shade(1), y_shade(1), y_shade(2), y_shade(2)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% Add text 'Transition region' at the middle of the shaded region
text(mean(x_shade), 190, 'Transition region', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'FontSize', 11); % Adjust the font weight and size of the text
legend('M_{\theta_c}', 'M_{\delta_e}', 'Location', 'southeast', 'FontSize', 12);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STACKED BAR PLOT ON M_W_COMPONENTS THROUGHOUT ENTIRE FORWARD VELOCITY SPEED RANGE

% figure(21)
% hp = bar(V_vals_kts, M_w_components', 'stacked');
% cm = colororder;
% xlabel('Forward speed u [m/s]'); ylabel('M_w [rad/s m]');
% legend('M_{MR_u}', 'M_{MR_l}', 'M_{hinge}', 'M_{flap}', 'M_{fus}', 'M_{ht}', 'M_e');
% title('Contributions to the static stability derivative M_w');
% grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%

%% M_W COMPONENTS BAR PLOT FOR SELECTED TRIM SPEEDS -- FOR COMPARISON WITH NLR PLOTS

M_w_MR_u = M_w_components(1,:); M_w_MR_l = M_w_components(2,:); 
M_w_hinge = M_w_components(3,:); M_w_flap = M_w_components(4,:); 
M_w_fus = M_w_components(5,:); M_w_ht = M_w_components(6,:); 
M_w_e = M_w_components(7,:);

figure(22)
bary = [M_w(1), M_w(6), M_w(11), M_w(16), M_w(21), M_w(26); 
    M_w_MR_u(1), M_w_MR_u(6), M_w_MR_u(11), M_w_MR_u(16), M_w_MR_u(21), M_w_MR_u(26);
    M_w_MR_l(1), M_w_MR_l(6), M_w_MR_l(11), M_w_MR_l(16), M_w_MR_l(21), M_w_MR_l(26);
    M_w_hinge(1), M_w_hinge(6), M_w_hinge(11), M_w_hinge(16), M_w_hinge(21), M_w_hinge(26);
    M_w_flap(1), M_w_flap(6), M_w_flap(11), M_w_flap(16), M_w_flap(21), M_w_flap(26);
    M_w_fus(1), M_w_fus(6), M_w_fus(11), M_w_fus(16), M_w_fus(21), M_w_fus(26);
    M_w_ht(1), M_w_ht(6), M_w_ht(11), M_w_ht(16), M_w_ht(21), M_w_ht(26)];
b = bar(bary);
s=6;
colours = [];
for j = 1:numel(b)
    colours = [colours; [1-j/s, j/s, 1]]; % Adjust the color definition as needed
    set(b(j), 'FaceColor', colours(j,:))
end
legend('0 kts', '50 kts', '100 kts', '150 kts','200 kts', '250 kts');
name={'Total';'URotor';'LRotor';'Hinge';'Flap';'Fuselage';'HStab + Elev'};
set(gca,'xticklabel',name);
xlabel('Component Contribution to M_w'); ylabel('M_w');
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

w_phug2 = -g.*M_u./M_q;
a_phug = 1;
b_phug = - (X_u + g * M_u ./ (M_q.^2));
c_phug = w_phug2;
discriminant_phug = b_phug.^2 - 4*a_phug*c_phug;
lambda_phug_1 = (-b_phug + sqrt(discriminant_phug)) / (2*a_phug);
lambda_phug_2 = (-b_phug - sqrt(discriminant_phug)) / (2*a_phug);

%%% Short Period
U_e = V_vals;
w_sp2 = abs(Z_w.*M_q - (Z_q+U_e).*M_w);
eta_sp = -(Z_w+M_q)/(2.*sqrt(w_sp2));
a_sp = 1;
b_sp = -(Z_w+M_q);
c_sp = w_sp2;
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
    load('phugHQ1.mat'); plot(phugHQ1(:,1), phugHQ1(:,2), '-.k');
    load('phugHQ2.mat'); plot(phugHQ2(:,1), phugHQ2(:,2), '-.k');
    xlabel('\mu [1/sec]'); ylabel('\omega [rad/sec]'); title('Phugoid'); grid on;
    xlim([min(phugHQ1(:,1)), 0.2]); ylim([0, max(phugHQ2(:,2))]);
    hold off;

    subplot(2,2,2)
    scatter(real(lambda_sp_1), imag(lambda_sp_1), [], colors, 's', 'filled', 'DisplayName', 'sp');
    hold on;
    scatter(real(lambda_sp_2), imag(lambda_sp_2), [], colors, 's', 'filled');
    xlabel('\mu [1/sec]'); ylabel('\omega [rad/sec]'); title('Short Period'); grid on;hold off;

    subplot(2,2,3)
    scatter(M_q, 0, [], colors, 's', 'filled', 'DisplayName', 'Mq');
    xlabel('\mu [1/sec]'); ylabel('\omega [rad/sec]'); title('Pitch Subsidence'); grid on;

    subplot(2,2,4)
    scatter(Z_w, 0, [], colors, 's', 'filled', 'DisplayName', 'Mq');
    xlabel('\mu [1/sec]'); ylabel('\omega [rad/sec]'); title('Heave Subsidence'); grid on;

% figure()
% scatter(real(lambda_phug_1), imag(lambda_phug_1), [], colors, 's', 'filled', 'DisplayName', 'sp'); hold on; 
% scatter(real(lambda_phug_2), imag(lambda_phug_2), [], colors, 's', 'filled', 'DisplayName', 'sp'); hold off; grid on;

% figure()
% scatter(real(lambda_sp_1), imag(lambda_sp_1), [], colors, 's', 'filled', 'DisplayName', 'sp'); hold on; 
% scatter(real(lambda_sp_2), imag(lambda_sp_2), [], colors, 's', 'filled', 'DisplayName', 'sp'); hold off; grid on;

% springterm = Z_w.*M_q - (Z_q+V_vals).*M_w;
% figure()
% plot(springterm, '--*'); grid on;




%% CONTROL ALLOCATION

epsilon = 0.001;

for i=1:length(V)
    V_kts(i) = V(i)*1.944;

    if V_kts(i)<V_l_transit
        rho_theta_s(i) = 1;
        rho_delta_e(i) = epsilon;

    elseif V_kts(i)>=V_l_transit && V_kts(i)<V_u_transit
        rho_theta_s(i) = 1 - 1/(V_u_transit-V_l_transit)*(V_kts(i)-V_l_transit);
        rho_delta_e(i) = 1/(V_u_transit-V_l_transit)*(V_kts(i)-V_l_transit);

    else
        rho_theta_s(i) = epsilon;
        rho_delta_e(i) = 1;
    end
end

rho = [rho_theta_s; rho_delta_e];

% PLOT CONTROL ALLOCATION DISTRIBUTION FOR PITCH
% figure(41)
% bar(V_vals_kts(3:2:20), rho(:,3:2:20), 'stacked');
% xlabel('Forward velocity u [kts]'); ylabel('Control Weight');
% legend('\theta_s', '\delta_e');
% ylim([0,1]);
% title('Control allocation weights for pitch');
% grid on;
% figure(42)
% plot(V_vals_kts, rho_theta_s, V_vals_kts, rho_delta_e, '-.', 'LineWidth', 2); grid on; 
% ylabel('Control Weight'); xlabel('Forward velocity u [kts]');
% title('Control allocation weights for pitch');
% legend('\theta_s', '\delta_e')

%% Handling Quality Analysis
% figure(61)
% pitch_tf_fit = tf(4^2, [1 2*0.707*4 4^2]);
% bode(sys_lin{1, 1}(3,2)); hold on;
% % bode(sys_lin{1, 2}(3,2)); bode(sys_lin{1, 3}(3,2));
% [GM, PM, w_gm, w_pm] = margin(sys_lin{1, 1}(3,2));
% 
% bode(pitch_tf_fit , 'r--'); bode(tf(-13, [0.26 1]), 'k--');
% hold off; grid on;
% legend('linSys 0kts', 'command model', 'inverse model')

[mag, phase, wout] = bode(sys_lin{1, 1}(3,2));

% Convert magnitude to dB
magdB = 20*log10(mag);

% Convert phase to degrees
phase_deg = squeeze(phase); % Extract phase data
phase_deg = phase_deg(:); % Reshape to a column vector
phase_deg = phase_deg(1:length(wout)); % Extract valid part of phase data

% Find the frequency where phase is closest to 45 degrees
desiredPhase = 180-45; % desired phase in degrees
[~, idx] = min(abs(phase_deg - desiredPhase));

% Extract corresponding frequency
w_bw_phase = wout(idx);

% Display the result
disp(['Frequency for phase of 135 degrees: ', num2str(w_bw_phase), ' rad/s']);


%% Fan Plot
RPM = 0:1:400;
Omega = RPM.*(2*pi/60);
Omega_rotor_nominal = 35;
normalized_rotorspeed = Omega./Omega_rotor_nominal;

% hinge, no spring
e = 0.45515;       % TUNE THIS PARAMETER FOR CORRECT flapping frequency 
                % BERGER HAS 1.39/rev, NLR 1.5/rev
Southwell = 1+3*e/(2*(1-e));
omega_nr2 = 0;
flapfreq = omega_nr2 + Southwell*Omega.^2;    %in rad/s
flapfreq_norm = sqrt(flapfreq)./Omega_rotor_nominal;

% hinge and spring
e_spring_ratio = 0.04;
e_spring = e_spring_ratio*R;
K_b = 900000;       % Tunes the starting location of the plot
I_b = 1/3*m_bl*R^2;
M_b = -500;
frac = 1800;
Southwell = 1+3*e_spring/(2*(1-e_spring));
omega_nr2_spring = K_b/I_b;
% flapfreq_spring = omega_nr2_spring + Southwell*Omega.^2;    %in rad/s
flapfreq_spring = K_b/I_b + (1+e_spring*M_b/I_b)*Omega.^2;
flapfreq_norm_spring = sqrt(flapfreq_spring)./Omega_rotor_nominal;

figure(81);
plot(normalized_rotorspeed, flapfreq_norm, normalized_rotorspeed, flapfreq_norm_spring, 'LineWidth', 2);
hold on;

maxperrev = 5;
revlineyval = zeros(maxperrev, size(normalized_rotorspeed, 2));

for i = 1:maxperrev
    revlineyval(i,:) = i * normalized_rotorspeed;
end

plot(normalized_rotorspeed, revlineyval, 'k-.');
grid on;

% Add text labels
for i = 1:maxperrev
    text(normalized_rotorspeed(end), revlineyval(i, end), ['\bf', num2str(i), '/rev   '], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

[~, idx_closest_to_1] = min(abs(normalized_rotorspeed - 1));
flapmodefreq = flapfreq_norm(idx_closest_to_1);
disp(['Flap mode Frequency: ', num2str(flapmodefreq), '/rev']);

xlabel('Normalized Rotor Speed \omega/\omega_0 (-)');
ylabel('Frequency (Hz)'); ylim([0,3]);
title('1/rev, 2/rev, ... Lines');

yline(flapmodefreq, '--');

hold off;





















    