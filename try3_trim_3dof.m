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
% d_theta_0(1) = deg2rad(0.1); not needed in 3dof , only for 6dof for yaw
% stability
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
% x_k = zeros(4,length(V_vals));
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

    %% Rotor RPM Scheduling
    if V(i) < 70
        Omega = 40;
    else
        Omega = 40 - 5/30*(V(i)-70);
    end

    %%
    if i>1
        x_k(:,i) = x_k(:,i-1);      % initial guess of each trim iteration is the previous trim state - for quick convergence
    end

    %% Define State Variable Names
    theta_0_u(i) = x_k(1,i)/2; theta_0_l(i) = x_k(1,i)/2;
    theta_c(i) = x_k(3,i); lambda_0_u(i) = x_k(4,i); lambda_0_l(i) = x_k(5,i);
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
        f(:,i) = f_xk(vel(:,i), x_ki(:,i),  q(i),  theta_f(i));

        % calculate x_k+1 for perturbation in x(1,i)
        ti1 = [tistep; 0; 0; 0; 0; 0];
        x_k1(:,i) = x_ki(:,i) + ti1;
        % calculate f_xk1_ti(:,i) using x_k+1 for FIRST column of jacobian
        f_xk1_ti(:,i) = f_xk(vel(:,i), x_k1(:,i), q(i),  theta_f(i));

        % calculate x_k+1 for perturbation in x(2,i)
        ti2 = [0; tistep; 0; 0; 0; 0];
        x_k2(:,i) = x_ki(:,i) + ti2;
        % calculate f_xk2_ti(:,i) using x_k+1 for SECOND column of jacobian
        f_xk2_ti(:,i) = f_xk(vel(:,i), x_k2(:,i), q(i),  theta_f(i));
        
        % calculate x_k+1 for perturbation in x(3,i)
        ti3 = [0; 0; tistep; 0; 0; 0];
        x_k3(:,i) = x_ki(:,i) + ti3;
        % calculate f_xk_ti(:,i) using x_k+1 for THIRD column of jacobian
        f_xk3_ti(:,i) = f_xk(vel(:,i), x_k3(:,i), q(i),  theta_f(i));

        % calculate x_k+1 for perturbation in x(4,i)
        ti4 = [0; 0; 0; tistep; 0; 0];
        x_k4(:,i) = x_ki(:,i) + ti4;
        % calculate f_xk_ti(:,i) using x_k+1 for THIRD column of jacobian
        f_xk4_ti(:,i) = f_xk(vel(:,i), x_k4(:,i), q(i),  theta_f(i));

        % calculate x_k+1 for perturbation in x(5,i)
        ti5 = [0; 0; 0; 0; tistep; 0];
        x_k5(:,i) = x_ki(:,i) + ti5;
        % calculate f_xk_ti(:,i) using x_k+1 for THIRD column of jacobian
        f_xk5_ti(:,i) = f_xk(vel(:,i), x_k5(:,i), q(i),  theta_f(i));

        % calculate x_k+1 for perturbation in x(5,i)
        ti6 = [0; 0; 0; 0; 0; tistep];
        x_k6(:,i) = x_ki(:,i) + ti6;
        % calculate f_xk_ti(:,i) using x_k+1 for THIRD column of jacobian
        f_xk6_ti(:,i) = f_xk(vel(:,i), x_k6(:,i), q(i),  theta_f(i));
        
        % assemble all 3 columns for Jacobian
        ti = tistep * ones(size(f(:,i)));
        Jac = [(f_xk1_ti(:,i)-f(:,i))./ti, (f_xk2_ti(:,i)-f(:,i))./ti, (f_xk3_ti(:,i)-f(:,i))./ti, ...
            (f_xk4_ti(:,i)-f(:,i))./ti, (f_xk5_ti(:,i)-f(:,i))./ti, (f_xk6_ti(:,i)-f(:,i))./ti];
        
        % calculate new input vector according to Jacobian
        x_k(:,i) = x_ki(:,i) - inv(Jac) * f(:,i);

        % f(:,i)
        [f(:,i), output_other(:,i)] = f_xk(vel(:,i), x_k(:,i), q(i), theta_f(i));

        f_sum(count(i)) = sum(abs(f(:,i)),"all");
        
       

    end
    % Jac

    

    %%
    % Defining the final calculated values, not necessary for trim but nice
    % for analysis and debugging purposes
    lambda_u(i) = (mu_z(i) - x_k(4,i))/2;
    % lambda_0_l(i) = (1+(-3.81*mu_x(i)+1.45))*x_k(4,i);
    lambda_0_l(i) = 2*x_k(4,i);
    lambda_l(i) = (mu_z(i) - lambda_0_l(i))/2;
    
    C_T_Glau_u(i) = x_k(4,i) * sqrt(mu(i)^2+(mu_z(i)+x_k(4,i))^2);
    C_T_BEM_u(i) = sigma*c_l_a/2*( (1/3+mu_x(i)^2/2)*x_k(1,i)/2 + (1+mu_x(i)^2)/8*twist_r + lambda_u(i)/2 );
    C_T_BEM_l(i) = sigma*c_l_a/2*( (1/3+mu_x(i)^2/2)*x_k(1,i)/2 + (1+mu_x(i)^2)/8*twist_r + lambda_l(i)/2 );
    
    T_u(i) = C_T_BEM_u(i) * rho * (Omega*R)^2 * pi*R^2;
    % T_l(i) = C_T_BEM_l(i) * rho * (Omega*R)^2 * pi*R^2;
    % if (0<=mu(i)) && (mu(i)<=0.381)
    %     delta_u2l(i) = -3.81*mu(i) + 1.45;    
    % else
    %     % delta_u2l = -3.81*mu(i) + 1.45;      % not correct but 0 gives issues
    %     delta_u2l(i) = 0.2;                         % 0 gives jumps in C_T and a_1
    % end
    T_l(i) = T_u(i);
    
    % Flapping Angles
    a0_u(i) = lock/(8*nu_b2) * (x_k(1,i)/2*(1+mu_x(i)^2) + 4/3*lambda_u(i) + twist_r*(4/5+2/3*mu_x(i)) - 4/3*mu_x(i)*x_k(2,i) );
    a0_l(i) = lock/(8*nu_b2) * (x_k(1,i)/2*(1+mu_x(i)^2) + 4/3*lambda_l(i) + twist_r*(4/5+2/3*mu_x(i)) - 4/3*mu_x(i)*x_k(2,i) );
    a1_u(i) = (8/3*mu_x(i)*x_k(1,i)/2 + 2*mu_x(i)*lambda_u(i) ...
        - 16/lock*q(i)/Omega + 2*twist_r*mu_x(i) - (1+3/2*mu_x(i)^2)*x_k(2,i)) / (1-1/2*mu_x(i)^2); 
    a1_l(i) = (8/3*mu_x(i)*x_k(1,i)/2 + 2*mu_x(i)*lambda_l(i) ...
        - 16/lock*q(i)/Omega + 2*twist_r*mu_x(i) - (1+3/2*mu_x(i)^2)*x_k(2,i)) / (1-1/2*mu_x(i)^2); 

    % Rotor plane angles
    alfa_sp(i) = atan2(w(i),u(i));
    alfa_cp(i) = alfa_sp(i) + x_k(2,i);
    alfa_dp_u(i) = alfa_cp(i) - a1_u(i);
    alfa_dp_l(i) = alfa_cp(i) - a1_l(i);
    a1r_u = x_k(2,i) - a1_u(i);
    a1r_l = x_k(2,i) - a1_l(i);

    udot(i) = f(1,i);
    wdot(i) = f(2,i);
    qdot(i) = f(3,i);
    lambda_0_u_dot(i) = f(4,i);

    u(i+1) = V(i)*cos(atan2(w(i),u(i)));
    w(i+1) = V(i)*sin(atan2(w(i),u(i)));
    q(i+1) = 0;
    
    % H-Forces
    CHdp_u(i) = sigma*0.015*mu_x(i)/4 + sigma*c_l_a/4*((a1_u(i)*mu_x(i)^2/2 + mu_x(i)*lambda_u(i))*x_k(1,i)/2 ...
        + q(i)/Omega*(-a0_u(i)/3) + (a0_u(i)^2+a1_u(i)^2)*mu_x(i)/2);
    CHdp_l(i) = sigma*0.015*mu_x(i)/4 + sigma*c_l_a/4*((a1_l(i)*mu_x(i)^2/2 + mu_x(i)*lambda_l(i))*x_k(1,i)/2 ...
        + q(i)/Omega*(-a0_l(i)/3) + (a0_l(i)^2+a1_l(i)^2)*mu_x(i)/2);

    Hdp_u(i) = CHdp_u(i) * rho * (Omega*R)^2 * pi*R^2;
    Hdp_l(i) = CHdp_l(i) * rho * (Omega*R)^2 * pi*R^2;

    % Fuselage Drag Force
    D_fus(i) = 1/2*rho*V(i)^2*CDS;

    % Pusher Prop Testing
    alfa_sp_p(i) = atan2(w(i),u(i));

    mu_xp(i) = sqrt(v(i)^2 + (w(i) + Kp*Omega*R*(x_k(4,i)+lambda_0_l(i)) + q(i)*l_p)^2)/(omega_p*R_p);
    lambda_p(i) = -(u(i)-q(i)*h_p-r(i)*d_p)/(omega_p*R_p) - x_k(5,i);

    C_T_Glau_p(i) = 2*x_k(5,i)*sqrt( (mu(i)*sin(alfa_sp_p(i)))^2 + ...
        (mu(i)*cos(alfa_sp_p(i))+x_k(5,i))^2 );       % this is like climbing flight for a rotor
    C_T_BEM_p(i) = sigma_p*c_l_a_p/2*( (1/3+mu_xp(i)^2/2)*x_k(3,i) + (1+mu_xp(i)^2)/8*twist_p +  lambda_p(i)/2 );
    Tp(i) = C_T_BEM_p(i) * rho * (omega_p*R_p)^2 * pi*R_p^2;

    L_h = 0;

end

% PLOTTING
%%
% figure(10)
% subplot(4,1,1)
% plot(mu, rad2deg(x_k(1,:)), '-*'), grid;
% legend('\theta_0'); title('Mean Collective')
% subplot(4,1,2)
% plot(mu, rad2deg(x_k(2,:)), '-*'), grid;
% legend('\theta_c'); title('Longitudinal Cyclic')
% subplot(4,1,3)
% plot(mu, rad2deg(x_k(3,:)), '-*'), grid;
% legend('\theta p'); title('Pusher Prop Collective')
% subplot(4,1,4)
% plot(mu, rad2deg(theta_f), '-*'), grid;
% legend('\theta f'); xlabel("Advance Ratio \mu [-]"); title('Fuselage Pitch (Scheduled)')
% 
% figure(11)
% subplot(2,1,1)
% plot(mu, x_k(4,:), mu, x_k(5,:)), grid;
% legend('lambda_0_u', 'lambda_0_l');
% subplot(2,1,2)
% plot(mu, T_u, mu, T_l), grid;
% legend('T_u', 'T_l');

figure(10)
plot(u(1:end-1)*1.94384, output_other(1,:)); 
hold on;
yline(0.9, '--');
grid on;
ylabel('M_{tip} , Omega/Omega_{hover}'); xlabel('Forward Speed u [kts]');
hold off;


%%
% figure(2)
% subplot(3,1,1)
% plot(mu, Hdp_u, mu, Hdp_l), grid, title('H-Forces');
% ylabel('H-Force [N]'), xlabel('Advance Ratio \mu [-]');
% legend('H-force U', 'H-force L');
% subplot(3,1,2)
% plot(mu, T_u, '-*', mu, T_l, '--o', mu, W, ...
%     mu, (T_u+T_l), '--o'), grid, title('Thrust Forces and Weight');
% ylabel('Rotor Thrust [N]'), xlabel('Advance Ratio \mu [-]');
% legend('T u', 'T l', 'Weight', 'T tot')
% subplot(3,1,3)
% plot(mu, Tp, '-*', mu, D_fus+Hdp_l+Hdp_u, mu, T_u.*sin(alfa_dp_u+gamma_s), '-o', mu, T_l.*sin(alfa_dp_l+gamma_s), '--o');
% grid, title('Propeller Thrust Compared to Drag');
% ylabel('Force[N]'), xlabel('Advance Ratio \mu [-]');
% legend('T p', 'Drag (Dfus+Hl+Hu) [N]');

%%
% figure(2)
% plot(mu, rad2deg(theta_p), '-*'), grid, title('Pusher Propeller Collective')
% legend('\theta p', 'Location', 'southeast'); xlabel("Advance Ratio \mu [-]")

%%
% figure(2)
% plot(mu, x_k(4,:), mu, 2*x_k(4,:),mu,lambda_c, ...
%     mu, sqrt(-V.^2/2 + sqrt((V.^2/2).^2 + vi_hov^4))/(Omega*R), ...
%     mu, 2*sqrt(-V.^2/2 + sqrt((V.^2/2).^2 + vi_hov^4))/(Omega*R) ),grid
% legend('lambda i u mean', 'lambda i l mean', 'lambda c', '^4 upper', '^4 lower'); xlabel("Advance Ratio \mu [-]")

%%
% figure(3)
% plot(mu, C_T_Glau_u, '-*', mu, C_T_Glau_l, mu, C_T_BEM_u, mu, C_T_BEM_l),grid
% legend('CT Glau lower', 'CT Glau lower', 'CT BEM upper', 'CT BEM lower'); xlabel("Advance Ratio \mu [-]")

% figure(4)
% plot(V, C_T_Glau_u, '-*', V, C_T_Glau_l, V, C_T_BEM_u, V, C_T_BEM_l, V, CTp),grid
% legend('CT Glau upper', 'CT Glau lower', 'CT BEM upper', 'CT BEM lower', 'C T Prop'); xlabel("Advance Ratio \mu [-]")

%%
% figure(5)
% plot(mu, T_u/g, mu, T_l/g),grid
% legend('T upper [kg]', 'T lower [kg]'); xlabel("Advance Ratio \mu [-]")

% figure(6)
% plot(V, T_u/g, V, T_l/g),grid
% legend('T upper [kg]', 'T lower [kg]'); xlabel("Speed V [-]")

%%
% figure(7)
% plot(mu, rad2deg(a1_u), '--', mu, rad2deg(a1_l), '--', ...
%     mu, rad2deg(alfa_cp), '-*', mu, rad2deg(alfa_sp), '-*', ...
%     mu, rad2deg(alfa_dp_u), '-*', mu, rad2deg(alfa_dp_l), '-*',...
%     mu, rad2deg(x_k(2,:)), '--o')
% % ylim([-10, 5])
% legend("a 1 u", "a 1 l ", ...
%     "\alpha_{cp}", "\alpha_{sp}", ...
%     "\alpha_{dp u}" , "\alpha_{dp l}", ...
%     "\theta_c", 'Location', 'northwest'), grid
% xlabel("Advance Ratio \mu [-]")


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
    f_lin_dist_0_pos = f_xk(vel(:,i)+vel_lin0, trim_vals_0, q(i),  theta_f(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk(vel(:,i)-vel_lin0, trim_vals_0, q(i),  theta_f(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_u(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; lin_step; 0; 0; 0];
    f_lin_dist_0_pos = f_xk(vel(:,i)+vel_lin0, trim_vals_0, q(i),  theta_f(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk(vel(:,i)-vel_lin0, trim_vals_0, q(i),  theta_f(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_w(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0, q(i)+0.01,  theta_f(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0, q(i)-0.01,  theta_f(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_q(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    % Row 2 (Z-force --> wdot state)
    vel_lin0 = [lin_step; 0; 0; 0; 0];
    f_lin_dist_0_pos = f_xk(vel(:,i)+vel_lin0, trim_vals_0, q(i),  theta_f(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk(vel(:,i)-vel_lin0, trim_vals_0, q(i),  theta_f(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Z_u(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; lin_step; 0; 0; 0];
    f_lin_dist_0_pos = f_xk(vel(:,i)+vel_lin0, trim_vals_0, q(i),  theta_f(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk(vel(:,i)-vel_lin0, trim_vals_0, q(i),  theta_f(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Z_w(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0, q(i)+0.01,  theta_f(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0, q(i)-0.01,  theta_f(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Z_q(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    % Row 3 (M-moment --> qdot state)
    vel_lin0 = [lin_step; 0; 0; 0; 0];
    f_lin_dist_0_pos = f_xk(vel(:,i)+vel_lin0, trim_vals_0, q(i),  theta_f(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk(vel(:,i)-vel_lin0, trim_vals_0, q(i),  theta_f(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    M_u(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    vel_lin0 = [0; lin_step; 0; 0; 0];
    f_lin_dist_0_pos = f_xk(vel(:,i)+vel_lin0, trim_vals_0, q(i),  theta_f(i));
    f_lin_dist_0_pos_Mw = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk(vel(:,i)-vel_lin0, trim_vals_0, q(i),  theta_f(i));
    f_lin_dist_0_neg_Mw = f_lin_dist_0_neg(3,1);
    M_w(i) = (f_lin_dist_0_pos_Mw - f_lin_dist_0_neg_Mw)/(2*lin_step);

    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0, q(i)+0.01,  theta_f(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0, q(i)-0.01,  theta_f(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    M_q(i) = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*0.01);

    % Assemble A-matrix
    A_lin = [X_u(i), X_w(i), X_q(i);
            Z_u(i), Z_w(i), Z_q(i);
            M_u(i), M_w(i), M_q(i)]

    %%%%%%%%%%%%%%%%%%% B-MATRIX %%%%%%%%%%%%%%%%%%%
    lin_step = deg2rad(1);

    % Row 1 (X-force --> qdot state)
    b_step = [lin_step; 0; 0; 0; 0; 0];   % input in collective
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i),  theta_f(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i),  theta_f(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_theta_0 = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; lin_step; 0; 0; 0; 0];   % input in forward cyclic
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i),  theta_f(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i),  theta_f(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_theta_c = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; lin_step; 0; 0; 0];
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i),  theta_f(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(1,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i),  theta_f(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(1,1);
    X_theta_p = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    % Row 2 (Z-force --> qdot state)
    b_step = [lin_step; 0; 0; 0; 0; 0];   % input in collective
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i),  theta_f(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i),  theta_f(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Z_theta_0 = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; lin_step; 0; 0; 0; 0];   % input in forward cyclic
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i),  theta_f(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i),  theta_f(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Z_theta_c = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; lin_step; 0; 0; 0];
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i),  theta_f(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(2,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i),  theta_f(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(2,1);
    Z_theta_p = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    % Row 3 (M-moment --> qdot state)
    b_step = [lin_step; 0; 0; 0; 0; 0];   % input in collective
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i),  theta_f(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i),  theta_f(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    M_theta_0 = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; lin_step; 0; 0; 0; 0];   % input in forward cyclic
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i),  theta_f(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i),  theta_f(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    M_theta_c = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    b_step = [0; 0; lin_step; 0; 0; 0];
    f_lin_dist_0_pos = f_xk(vel(:,i), trim_vals_0+b_step, q(i),  theta_f(i));
    f_lin_dist_0_pos = f_lin_dist_0_pos(3,1);
    f_lin_dist_0_neg = f_xk(vel(:,i), trim_vals_0-b_step, q(i),  theta_f(i));
    f_lin_dist_0_neg = f_lin_dist_0_neg(3,1);
    M_theta_p = (f_lin_dist_0_pos - f_lin_dist_0_neg)/(2*lin_step);

    % Assemble B-MATRIX
    B_lin = [X_theta_0, X_theta_c, X_theta_p;
            Z_theta_0, Z_theta_c, Z_theta_p;
            M_theta_0, M_theta_c, M_theta_p]

    % Assemble C-Matrix
    C_lin = eye(size(A_lin))

    % Assemble D-Matrix
    D_lin = zeros(size(B_lin))       % #columns = #inputs
    
    % Analysis with state-space
    % sys_lin = ss(A_lin, B_lin, C_lin, D_lin);
    % pitch_tf_fit = tf([4^2], [1 2*0.707*4 4^2])
    % bode(sys_lin(3,2)); grid; title('q/\theta_{1s} (50 m/s)');
    % hold on;
    % bode(pitch_tf_fit); hold off;
    % pzplot(sys_lin)
    
    % Save Current Trim SS System in lin_systems
    lin_systems{i, 1} = A_lin;
    lin_systems{i, 2} = B_lin;
end

V_vals_kts = V_vals.*1.944;

figure(1)
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


% %% Nonlinear Simulation
% % Time Paramaters
% simN = 400;
% tEnd = 15;
% dt = (tEnd-t(1))/simN;
% tau = 0.1;
% 
% U = zeros(5,simN);
% U(1,1) = x_k(1,1);  % collective at hover trim
% U(2,1) = x_k(2,1);  % longitudinal cyclic at hover trim
% U(3,1) = x_k(3,1);  % prop collective at hover trim
% % U(4,1) = x_k(4,1);
% % U(5,1) = x_k(5,1);
% t(1)=0;
% u(1)=0;
% w(1)=0;
% q(1)=0;
% theta(1)=0;
% 
% for i = 1:simN
% 
%     % Constant Rotor Collective and prop collective
%     U(1,i) = U(1,1);
%     U(3,i) = U(3,1);
%     % Cyclic step of 1deg between t=2s and t=2.5s
%     if t(i) < 2
%         U(2,i) = U(2,1);
%     end
%     if t(i) >= 2 && t(i) <= 2.5 
%        U(2,i) = U(2,1)+deg2rad(1);
%     else 
%        U(2,i) = 0;
%     end
% 
%     % Calculate Angles
%     alfa_sp(i) = atan2(w(i),u(i));
%     alfa_cp(i) = alfa_sp(i) + U(2,i);
% 
%     % Speed variables
%     V(i) = sqrt(u(i)^2+w(i)^2);
%     mu(i) = V(i)/(Omega*R);
%     mu_x(i) = V(i)/(Omega*R) * cos(alfa_cp(i));
%     mu_z(i) = V(i)/(Omega*R) * sin(alfa_cp(i));
%     lambda_c(i) = mu(i) * sin(alfa_cp(i));
% 
%     vel(:,i) = [u(i); w(i); mu(i); mu_x(i); mu_z(i)];
%     inflow(:,i) = [lambda_c(i); U(4,i)];
% 
%     xdot(:,i) = f_xk(vel(:,i), U(:,i), q(i), theta(i));
% 
%     % Euler Integration
%     u(i+1) = u(i)+dt*xdot(1,i);
%     w(i+1) = w(i)+dt*xdot(2,i);
%     q(i+1) = q(i)+dt*xdot(3,i);
%     theta(i+1) = theta(i)+dt*q(i);
%     U(4,i+1) = U(4,i)+dt*xdot(4,i);
%     U(5,i+1) = U(5,i)+dt*xdot(5,i);
%     t(i+1) = t(i)+dt;
% end
% 
% figure(1)
% subplot(3,1,1)
% plot(t, rad2deg(U(1,:)), t, rad2deg(U(2,:)), t, rad2deg(U(3,:))); grid on; legend('theta 0', 'theta 1s', 'theta p');
% subplot(3,1,2)
% plot(t, rad2deg(theta)); grid on; legend('theta')
% subplot(3,1,3)
% plot(t, u, t, w); grid on; legend('u', 'w')
% 
% figure(2)
% plot(t, U(2,:))
% 
% 
% %% Linear Simulation
% simN = 400;
% tEnd = 30;
% dt = (tEnd-t(1))/simN;
% 
% t = 0:dt:tEnd;
% U = zeros(3, simN);
% U(1,:) = x_k(1,1);
% step_index = t <= 1;
% U(2, step_index) = deg2rad(1);
% U(3,:) = x_k(3,1);
% [y, t, x] = lsim(ss(lin_systems{1,1}, lin_systems{1,2}, C_lin, D_lin), U, t(1:end-1));
% figure()
% plot(t, x);


%% ROTOR RPM - MACH PLOT



















    