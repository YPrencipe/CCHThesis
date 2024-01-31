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
lambda_0_u(1)=sqrt(mass*abs(g)/(area*2*rho))/(Omega_init*R);
lambda_p(1) = 0.002;
vi_hov = sqrt(1/2*mass*g/(2*rho*pi*R^2));

% Define trim variable vector
% x_k = zeros(4,length(V_vals));
x_k(:,1) = [theta_0(1); theta_c(1); theta_p(1); lambda_0_u(1); lambda_p(1)];


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
    theta_c(i) = x_k(3,i); lambda_0_u(i) = x_k(4,i);
    lambda_0_p = x_k(5,i);
   
    %%
    % Calculate Angles
    alfa_sp(i) = atan2(w(i),u(i));
    alfa_cp(i) = alfa_sp(i) + x_k(2,i);

    % Speed variables
    mu(i) = V(i)/(Omega*R);
    mu_x(i) = V(i)/(Omega*R) * cos(alfa_cp(i));
    mu_z(i) = V(i)/(Omega*R) * sin(alfa_cp(i));
    vel(:,i) = [u(i); w(i); V(i); mu(i); mu_x(i); mu_z(i)];

    %%
    lambda_c(i) = mu(i) * sin(alfa_cp(i));
   
    inflow(:,i) = [lambda_c(i); x_k(4,i)]; %; lambda_i_u_mean(i); lambda_i_l_mean(i)];

    %% Flapping Angles
    % a1_u = (8/3*mu_x(i)*theta_0_u(i) + 2*mu_x(i)*lambda_u - 16/lock*q(i)/Omega)/(1-1/2*mu_x(i)^2);
    % a1_l = (8/3*mu_x(i)*theta_0_l(i) + 2*mu_x(i)*lambda_l - 16/lock*q(i)/Omega)/(1-1/2*mu_x(i)^2);
    % a0_u = lock/8 * (theta_0_u(i)*(1+mu_x(i)^2)+4/3*lambda_u/2);
    % a0_l = lock/8 * (theta_0_l(i)*(1+mu_x(i)^2)+4/3*lambda_l/2);
    
    % flapping(:,i) = [a1_u; a1_l; a0_u; a0_l];
    flapping(:,i) = [0; 0; 0; 0];


    %%
    count(i) = 1;
    f_sum(1) = 10;


    % Newton-Rhapson Trim Algorithm 
    while f_sum(count(i)) > 0.0005 % convergence said to occur at 0.0005 according to leishman p 97

        count(i) = count(i)+1;
        tistep = deg2rad(0.01);          % small increment for Newton-Rhapson

        % calculate f(:,i) using x_k
        x_ki(:,i) = x_k(:,i);
        f(:,i) = f_xk(vel(:,i), x_ki(:,i), inflow(:,i), q(i), flapping(:,i), theta_f(i));

        % calculate x_k+1 for perturbation in x(1,i)
        ti1 = [tistep; 0; 0; 0; 0];
        x_k1(:,i) = x_ki(:,i) + ti1;
        % calculate f_xk1_ti(:,i) using x_k+1 for FIRST column of jacobian
        f_xk1_ti(:,i) = f_xk(vel(:,i), x_k1(:,i), inflow(:,i), q(i), flapping(:,i), theta_f(i));

        % calculate x_k+1 for perturbation in x(2,i)
        ti2 = [0; tistep; 0; 0; 0];
        x_k2(:,i) = x_ki(:,i) + ti2;
        % calculate f_xk2_ti(:,i) using x_k+1 for SECOND column of jacobian
        f_xk2_ti(:,i) = f_xk(vel(:,i), x_k2(:,i), inflow(:,i), q(i), flapping(:,i), theta_f(i));
        
        % calculate x_k+1 for perturbation in x(3,i)
        ti3 = [0; 0; tistep; 0; 0];
        x_k3(:,i) = x_ki(:,i) + ti3;
        % calculate f_xk_ti(:,i) using x_k+1 for THIRD column of jacobian
        f_xk3_ti(:,i) = f_xk(vel(:,i), x_k3(:,i), inflow(:,i), q(i), flapping(:,i), theta_f(i));

        % calculate x_k+1 for perturbation in x(4,i)
        ti4 = [0; 0; 0; tistep; 0];
        x_k4(:,i) = x_ki(:,i) + ti4;
        % calculate f_xk_ti(:,i) using x_k+1 for THIRD column of jacobian
        f_xk4_ti(:,i) = f_xk(vel(:,i), x_k4(:,i), inflow(:,i), q(i), flapping(:,i), theta_f(i));

        % calculate x_k+1 for perturbation in x(5,i)
        ti5 = [0; 0; 0; 0; tistep];
        x_k5(:,i) = x_ki(:,i) + ti5;
        % calculate f_xk_ti(:,i) using x_k+1 for THIRD column of jacobian
        f_xk5_ti(:,i) = f_xk(vel(:,i), x_k5(:,i), inflow(:,i), q(i), flapping(:,i), theta_f(i));
        
        % assemble all 3 columns for Jacobian
        ti = [tistep; tistep; tistep; tistep; tistep];
        Jac = [(f_xk1_ti(:,i)-f(:,i))./ti, (f_xk2_ti(:,i)-f(:,i))./ti, (f_xk3_ti(:,i)-f(:,i))./ti, ...
            (f_xk4_ti(:,i)-f(:,i))./ti, (f_xk5_ti(:,i)-f(:,i))./ti];
        
        % calculate new input vector according to Jacobian
        x_k(:,i) = x_ki(:,i) - inv(Jac) * f(:,i);

        % f(:,i)
        f(:,i) = f_xk(vel(:,i), x_k(:,i), inflow(:,i), q(i), flapping(:,i), theta_f(i));

        f_sum(count(i)) = sum(abs(f(:,i)),"all");
        
        

        % f(:,i) = f_xk(vel(:,i), x_k(:,i), inflow(:,i), q(i));

    end
    Jac


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
    alfa_dp_u(i) = alfa_cp(i)-a1_u(i);
    alfa_dp_l(i) = alfa_cp(i)-a1_l(i);

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

end

% PLOTTING
%%
figure(1)
subplot(4,1,1)
plot(mu, rad2deg(x_k(1,:)), '-*'), grid;
legend('\theta_0'); title('Mean Collective')
subplot(4,1,2)
plot(mu, rad2deg(x_k(2,:)), '-*'), grid;
legend('\theta_c'); title('Longitudinal Cyclic')
subplot(4,1,3)
plot(mu, rad2deg(x_k(3,:)), '-*'), grid;
legend('\theta p'); title('Pusher Prop Collective')
subplot(4,1,4)
plot(mu, rad2deg(theta_f), '-*'), grid;
legend('\theta f'); xlabel("Advance Ratio \mu [-]"); title('Fuselage Pitch (Scheduled)')


%%
figure(2)
subplot(3,1,1)
plot(mu, Hdp_u, mu, Hdp_l), grid, title('H-Forces');
ylabel('H-Force [N]'), xlabel('Advance Ratio \mu [-]');
legend('H-force U', 'H-force L');
subplot(3,1,2)
plot(mu, T_u, '-*', mu, T_l, '--o', mu, W, ...
    mu, (T_u+T_l), '--o'), grid, title('Thrust Forces and Weight');
ylabel('Rotor Thrust [N]'), xlabel('Advance Ratio \mu [-]');
legend('T u', 'T l', 'Weight', 'T tot')
subplot(3,1,3)
plot(mu, Tp, '-*', mu, D_fus+Hdp_l+Hdp_u, mu, T_u.*sin(alfa_dp_u+i_theta), '-o', mu, T_l.*sin(alfa_dp_l+i_theta), '--o');
grid, title('Propeller Thrust Compared to Drag');
ylabel('Force[N]'), xlabel('Advance Ratio \mu [-]');
legend('T p', 'Drag (Dfus+Hl+Hu) [N]');

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
figure(7)
plot(mu, rad2deg(a1_u), '--', mu, rad2deg(a1_l), '--', ...
    mu, rad2deg(alfa_cp), '-*', mu, rad2deg(alfa_sp), '-*', ...
    mu, rad2deg(alfa_dp_u), '-*', mu, rad2deg(alfa_dp_l), '-*',...
    mu, rad2deg(x_k(2,:)), '--o')
% ylim([-10, 5])
legend("a 1 u", "a 1 l ", ...
    "\alpha_{cp}", "\alpha_{sp}", ...
    "\alpha_{dp u}" , "\alpha_{dp l}", ...
    "\theta_c", 'Location', 'northwest'), grid
xlabel("Advance Ratio \mu [-]")

    
    
    
    
    