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
theta(1)=deg2rad(1);
theta_p(1) = deg2rad(0);
t(1)=0;
u(1)=0;
w(1)=0;
q(1)=0;
lambda_0_u(1)=sqrt(mass*abs(g)/(area*2*rho))/(Omega*R);
vi_hov = sqrt(1/2*mass*g/(2*rho*pi*R^2));

% Define speed vector
V_vals = [0.1, 1,2,3,4, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, ...
    70, 75, 80, 85, 90, 95, 100, 105, 110];
% V_vals = 0.1:80;

% Define trim variable vector
x_k = zeros(4,length(V_vals));
x_k(:,1) = [theta_0(1); theta_p(1); theta_c(1); lambda_0_u(1)];


%% Trim Routine

for i = 1:length(V_vals)
    i
    V(i) = V_vals(i);
    u(i) = V(i)*cos(atan2(w(i),u(i)));
    w(i) = V(i)*sin(atan2(w(i),u(i)));

    if i>1
        x_k(:,i) = x_k(:,i-1);      % initial guess of each trim iteration is the previous trim state - for quick convergence
    end
   
    % Calculate Angles
    alfa_sp(i) = atan2(w(i),u(i));
    alfa_cp(i) = x_k(3,i) - alfa_sp(i);

    % Speed variables
    mu(i) = V(i)/(Omega*R);
    mu_x(i) = V(i)/(Omega*R) * cos(alfa_cp(i));
    mu_z(i) = V(i)/(Omega*R) * sin(alfa_cp(i));
    vel(:,i) = [u(i); w(i); V(i); mu(i); mu_x(i); mu_z(i)];

    %%
    % % Define Polar Coordinates
    % psi_values = 0:deg2rad(1):2*pi;
    % r_values = 0:0.01:R;
    % [r, psi] = meshgrid(r_values, psi_values);
    % % Include Uniform (upper)/Linear (Lower) Inflow Distributions
    % ksi(i) = atan(mu_x(i)/(mu_z(i)+x_k(4,i)));
    % k_x = 15*pi/32 * tan(ksi(i)/2);        % Pitt and Peters estimates for k_x and k_y
    % k_y = 0;
    % 
    % % Interference Effects -- omitted for now due to numerical errors
    % if (0<=mu(i)) && (mu(i)<= 0.316)
    %     delta_l2u = -2.15*mu(i) + 0.68;
    % else
    %     % delta_l2u = -2.15*mu(i) + 0.68;     % this is not correct but 0 gives issues
    %     delta_l2u = 0;
    % end
    % if (0<=mu(i)) && (mu(i)<=0.381)
    %     delta_u2l = -3.81*mu(i) + 1.45;    
    % else
    %     % delta_u2l = -3.81*mu(i) + 1.45;      % not correct but 0 gives issues
    %     delta_u2l = 0;                         % 0 gives jumps in C_T and a_1
    % end
    % 
    % % Inflow interferences and linear distribution n stuff
    % lambda_i_u_prime = x_k(4,i)*(1 + 0*r.*cos(psi) + 0*r.*sin(psi));
    % lambda_i_l_prime = x_k(4,i)*(1 + k_x*r.*cos(psi) + k_y*r.*sin(psi));
    % 
    % lambda_i_u = lambda_i_u_prime; %+ lambda_i_l_prime*delta_l2u;
    % lambda_i_l = lambda_i_l_prime + lambda_i_u_prime*1; %+ lambda_i_u_prime*delta_u2l;
    % 
    % lambda_i_u_mean(i) = mean(lambda_i_u, "all");
    % lambda_i_l_mean(i) = mean(lambda_i_l, "all");

    % Mean inflows and forward flight inflow component
    % lambda_i_u_mean(i) = x_k(4,i);      % Equal to the upper rotor mean inflow lambda_0_u
    % lambda_i_l_mean(i) = 2*x_k(4,i);    % Double the upper rotor inflow
    lambda_c(i) = mu(i) * sin(alfa_cp(i));

    inflow(:,i) = [lambda_c(i); x_k(4,i)]; %; lambda_i_u_mean(i); lambda_i_l_mean(i)];

    %%
    count(i) = 1;
    f_sum(1) = 10;


    % Newton-Rhapson Trim Algorithm 
    while f_sum(count(i)) > 0.0001

        count(i) = count(i)+1;
        tistep = deg2rad(0.1);          % small increment for Newton-Rhapson

        % calculate f(:,i) using x_k
        x_ki(:,i) = x_k(:,i);
        f(:,i) = f_xk(vel(:,i), x_ki(:,i), inflow(:,i), q(i));

        % calculate x_k+1 for perturbation in x(1,i)
        ti1 = [tistep; 0; 0; 0];
        x_k1(:,i) = x_ki(:,i) + ti1;
        % calculate f_xk1_ti(:,i) using x_k+1 for FIRST column of jacobian
        f_xk1_ti(:,i) = f_xk(vel(:,i), x_k1(:,i), inflow(:,i), q(i));

        % calculate x_k+1 for perturbation in x(2,i)
        ti2 = [0; tistep; 0; 0];
        x_k2(:,i) = x_ki(:,i) + ti2;
        % calculate f_xk2_ti(:,i) using x_k+1 for SECOND column of jacobian
        f_xk2_ti(:,i) = f_xk(vel(:,i), x_k2(:,i), inflow(:,i), q(i));
        
        % calculate x_k+1 for perturbation in x(3,i)
        ti3 = [0; 0; tistep; 0];
        x_k3(:,i) = x_ki(:,i) + ti3;
        % calculate f_xk_ti(:,i) using x_k+1 for THIRD column of jacobian
        f_xk3_ti(:,i) = f_xk(vel(:,i), x_k3(:,i), inflow(:,i), q(i));

        % calculate x_k+1 for perturbation in x(3,i)
        ti4 = [0; 0; 0; tistep];
        x_k4(:,i) = x_ki(:,i) + ti4;
        % calculate f_xk_ti(:,i) using x_k+1 for THIRD column of jacobian
        f_xk4_ti(:,i) = f_xk(vel(:,i), x_k4(:,i), inflow(:,i), q(i));
        
        % assemble all 3 columns for Jacobian
        ti = [tistep; tistep; tistep; 0.001];
        Jac = [(f_xk1_ti(:,i)-f(:,i))./ti, (f_xk2_ti(:,i)-f(:,i))./ti, (f_xk3_ti(:,i)-f(:,i))./ti ...
            , (f_xk4_ti(:,i)-f(:,i))./ti];
        
        % calculate new input vector according to Jacobian
        x_k(:,i) = x_ki(:,i) - inv(Jac) * f(:,i);

        % f(:,i)
        f(:,i) = f_xk(vel(:,i), x_k(:,i), inflow(:,i), q(i));

        f_sum(count(i)) = sum(abs(f(:,i)),"all");
        
        

        % f(:,i) = f_xk(vel(:,i), x_k(:,i), inflow(:,i), q(i));

    end
    Jac

    % Defining the final calculated values, not necessary for trim but nice
    % for analysis and debugging purposes
    lambda_i_l(i) = ( 1 + (-3.81*mu(i) + 1.45) ) * x_k(4,i);
    C_T_Glau_u(i) = x_k(4,i) * sqrt(mu(i)^2+(mu_z(i)+x_k(4,i))^2);
    C_T_Glau_l(i) = 2*x_k(4,i) * sqrt(mu(i)^2+(mu_z(i)-lambda_i_l(i)+x_k(4,i))^2);
    C_T_BEM_u(i) = 1/4*c_l_a*sigma*(2/3*x_k(1,i)*(1+3/2*mu(i)^2)-(lambda_c(i)+x_k(4,i)));
    
    C_T_BEM_l(i) = 1/4*c_l_a*sigma*(2/3*x_k(1,i)*(1+3/2*mu(i)^2)-(lambda_c(i)+lambda_i_l(i)));

    T_u(i) = C_T_BEM_u(i) * rho * (Omega*R)^2 * pi*R^2;
    T_l(i) = ( 1 + (-3.81*mu_x(i) + 1.45) ) * T_u(i);

    a_1_u(i) = (8/3*mu_x(i)*x_k(1,i) - 2*mu_x(i)*(lambda_c(i)+x_k(4,i)) - 16/lock*q(i)/Omega) / (1-1/2*mu_x(i)^2);
    a_1_l(i) = (8/3*mu_x(i)*x_k(1,i) - 2*mu_x(i)*(lambda_c(i)+2*x_k(4,i)) - 16/lock*q(i)/Omega) / (1-1/2*mu_x(i)^2);

    udot(i) = f(1,i);
    wdot(i) = f(2,i);
    qdot(i) = f(3,i);
    lambda_0_u_dot(i) = f(4,i);

    % Euler Integration
    u(i+1) = V(i)*cos(atan2(w(i),u(i)));
    w(i+1) = V(i)*sin(atan2(w(i),u(i)));
    q(i+1) = 0;

    Drag(i) = 1/2*rho*V(i)^2*CDS;

    Tp_des(i) = 0.5 * Drag(i);   % prop accounts for 50% of drag produced by heli
    CTp(i) = Tp_des(i) / (rho * (Omega*R)^2 * pi*R^2);
    lambda_p(i) = sqrt(CTp(i)/2);

    CTp(i) = 1/2*sigma*C_l_p*(x_k(2,i)/3 - lambda_p(i)/2);
    T_p(i) = CTp(i) * (rho * (Omega*R)^2 * pi*R^2);

end

theta_0_mean = 2*x_k(1,:);  % Multiply the collective of upper by 2 because
                            % the lower one also has the same theta

% PLOTTING
%%
figure(1)
plot(mu, rad2deg(theta_0_mean), '-*', mu, rad2deg(x_k(2,:)), '-*', mu, rad2deg(x_k(3,:)), '-*'),grid;
% ylim([-5, 30])
legend('\theta_0', '\theta f', '\theta_c', 'Location', 'southwest'); xlabel("Advance Ratio \mu [-]")

%%
figure(2)
plot(mu, x_k(4,:), mu, 2*x_k(4,:),mu,lambda_c, ...
    mu, sqrt(-V.^2/2 + sqrt((V.^2/2).^2 + vi_hov^4))/(Omega*R), ...
    mu, 2*sqrt(-V.^2/2 + sqrt((V.^2/2).^2 + vi_hov^4))/(Omega*R) ),grid
legend('lambda i u mean', 'lambda i l mean', 'lambda c', '^4 upper', '^4 lower'); xlabel("Advance Ratio \mu [-]")

%%
% figure(3)
% plot(mu, C_T_Glau_u, '-*', mu, C_T_Glau_l, mu, C_T_BEM_u, mu, C_T_BEM_l),grid
% legend('CT Glau lower', 'CT Glau lower', 'CT BEM upper', 'CT BEM lower'); xlabel("Advance Ratio \mu [-]")

figure(4)
plot(V, C_T_Glau_u, '-*', V, C_T_Glau_l, V, C_T_BEM_u, V, C_T_BEM_l, V, CTp),grid
legend('CT Glau upper', 'CT Glau lower', 'CT BEM upper', 'CT BEM lower', 'C T Prop'); xlabel("Advance Ratio \mu [-]")

%%
% figure(5)
% plot(mu, T_u/g, mu, T_l/g),grid
% legend('T upper [kg]', 'T lower [kg]'); xlabel("Advance Ratio \mu [-]")

figure(6)
plot(V, T_u/g, V, T_l/g),grid
legend('T upper [kg]', 'T lower [kg]'); xlabel("Speed V [-]")

% IMPORTANT --- SPEED AND ADVANCE RATIO GIVE DIFFERENT RESULTS FOR LOWER
% ROTOR????

%%
figure(7)
plot(mu, rad2deg(a_1_u), '--', mu, rad2deg(a_1_l), '--', ...
    mu, rad2deg(alfa_cp), '-*', mu, rad2deg(alfa_sp), '-*', ...
    mu, rad2deg(x_k(3,:)), '-o')
legend("a 1 u", "a 1 l ", "\alpha_{dp}", "\alpha_{sp}", "\theta_c", 'Location', 'northwest'), grid
% xlabel("Advance Ratio \mu [-]")

    
    
    
    
    