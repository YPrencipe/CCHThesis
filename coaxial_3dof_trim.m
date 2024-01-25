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
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Set-up
% Load coaxial helicopter parameters
coaxial_heli_parameters;

% Initial Values
theta_0_u(1) = deg2rad(1);
theta_0_l(1) = deg2rad(1);
theta_c(1) = deg2rad(1);
t(1)=0;
u(1)=0;
w(1)=0;
q(1)=0;
theta(1)=deg2rad(0);
lambda_0_u(1)=sqrt(mass*abs(g)/(area*2*rho))/(Omega*R);

% Define speed vector
u_vals = [0.1, 0.2, 0.3,0.4,0.5,0.6, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80];

% Define trim variable vector
x_k = zeros(4,length(u_vals));
x_k(:,1) = [theta_0_u(1); theta_0_l(1); theta_c(1); lambda_0_u(1)];


%% Trim Routine

for i = 1:length(u_vals)
    i
    u(i) = u_vals(i);   % cycling through each speed defined in the speed range vector u_vals
    w(i) = 0;           % w set to 0 during trim procedure
    theta(i) = 0;       % fix theta_fuselage to 0 for all flightspeeds for trimming purposes

    if i>1
        x_k(:,i) = x_k(:,i-1);      % initial guess of each trim iteration is the previous trim state - for quick convergence
    end
   
    % Calculate Angles
    if u(i)==0 	
        if w(i)>0 	
            phi(i)=pi/2;
        else 
            phi(i)=-pi/2;
        end
    else
        phi(i)=atan(w(i)/u(i));
    end
    
    if u(i)<0
        phi(i)=phi(i)+pi;
    end
    
    alfa_c(i) = x_k(3,i) - phi(i);
    % alfa_c(i) = 0;

    % Speed variables
    V(i) = sqrt(u(i)^2+w(i)^2);
    mu(i) = V(i)/(Omega*R) * cos(alfa_c(i));
    mu_x(i) = u(i)/(Omega*R);
    mu_z(i) = w(i)/(Omega*R);
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
    lambda_i_u_mean(i) = x_k(4,i);      % Equal to the upper rotor mean inflow lambda_0_u
    lambda_i_l_mean(i) = 2*x_k(4,i);    % Double the upper rotor inflow
    lambda_c(i) = V(i)/(Omega*R) * sin(alfa_c(i));

    inflow(:,i) = [lambda_c(i); x_k(4,i); lambda_i_u_mean(i); lambda_i_l_mean(i)];

    %%
    count(i) = 1;
    f_sum(1) = 10;


    % Newton-Rhapson Trim Algorithm 
    while f_sum(count(i)) > 0.001

        count(i) = count(i)+1;
        tistep = deg2rad(0.1);          % small increment for Newton-Rhapson

        % calculate f(:,i) using x_k
        x_ki(:,i) = x_k(:,i);
        f(:,i) = f_xk(vel(:,i), x_ki(:,i), inflow(:,i), q(i), theta(i));

        % calculate x_k+1 for perturbation in x(1,i)
        ti1 = [tistep; 0; 0; 0];
        x_k1(:,i) = x_ki(:,i) + ti1;
        % calculate f_xk1_ti(:,i) using x_k+1 for FIRST column of jacobian
        f_xk1_ti(:,i) = f_xk(vel(:,i), x_k1(:,i), inflow(:,i), q(i), theta(i));

        % calculate x_k+1 for perturbation in x(2,i)
        ti2 = [0; tistep; 0; 0];
        x_k2(:,i) = x_ki(:,i) + ti2;
        % calculate f_xk2_ti(:,i) using x_k+1 for SECOND column of jacobian
        f_xk2_ti(:,i) = f_xk(vel(:,i), x_k2(:,i), inflow(:,i), q(i), theta(i));
        
        % calculate x_k+1 for perturbation in x(3,i)
        ti3 = [0; 0; tistep; 0];
        x_k3(:,i) = x_ki(:,i) + ti3;
        % calculate f_xk_ti(:,i) using x_k+1 for THIRD column of jacobian
        f_xk3_ti(:,i) = f_xk(vel(:,i), x_k3(:,i), inflow(:,i), q(i), theta(i));

        % calculate x_k+1 for perturbation in x(3,i)
        ti4 = [0; 0; 0; tistep];
        x_k4(:,i) = x_ki(:,i) + ti4;
        % calculate f_xk_ti(:,i) using x_k+1 for THIRD column of jacobian
        f_xk4_ti(:,i) = f_xk(vel(:,i), x_k4(:,i), inflow(:,i), q(i), theta(i));
        
        % assemble all 3 columns for Jacobian
        ti = [tistep; tistep; tistep; tistep];
        Jac = [(f_xk1_ti(:,i)-f(:,i))./ti, (f_xk2_ti(:,i)-f(:,i))./ti, (f_xk3_ti(:,i)-f(:,i))./ti ...
            , (f_xk4_ti(:,i)-f(:,i))./ti];
        
        % calculate new input vector according to Jacobian
        x_k(:,i) = x_ki(:,i) - inv(Jac) * f(:,i);

        % f(:,i)
        f(:,i) = f_xk(vel(:,i), x_k(:,i), inflow(:,i), q(i), theta(i));

        f_sum(count(i)) = sum(abs(f(:,i)),"all");
        
    end
    Jac

    % Defining the final calculated values, not necessary for trim but nice
    % for analysis and debugging purposes
    lambda_i_u_mean(i) = x_k(4,i);
    lambda_i_l_mean(i) = 2*x_k(4,i);

    C_T_Glau_u(i) = lambda_i_u_mean(i) * sqrt(mu(i)^2+(mu_z(i)+lambda_i_u_mean(i))^2);
    C_T_Glau_l(i) = 2*lambda_i_l_mean(i) * sqrt(mu(i)^2+(mu_z(i)-lambda_i_l_mean(i)+lambda_i_u_mean(i))^2);
    C_T_BEM_u(i) = 1/4*c_l_a*sigma*(2/3*x_k(1,i)*(1+3/2*mu(i)^2)-(lambda_c(i)+x_k(4,i)));
    C_T_BEM_l(i) = 1/4*c_l_a*sigma*(2/3*x_k(2,i)*(1+3/2*mu(i)^2)-(lambda_c(i)+x_k(4,i)));

    T_u(i) = C_T_BEM_u(i) * rho * (Omega*R)^2 * pi*R^2;
    T_l(i) = C_T_BEM_l(i) * rho * (Omega*R)^2 * pi*R^2;

    a_1_u(i) = (8/3*mu_x(i)*x_k(1,i) - 2*mu_x(i)*(lambda_c(i)+lambda_i_u_mean(i)) - 16/lock*q(i)/Omega) / (1-1/2*mu_x(i)^2);
    a_1_l(i) = (8/3*mu_x(i)*x_k(2,i) - 2*mu_x(i)*(lambda_c(i)+lambda_i_l_mean(i)) - 16/lock*q(i)/Omega) / (1-1/2*mu_x(i)^2);

    udot(i) = -g*sin(theta(i)) - CDS/mass*0.5*rho*u(i)*V(i) + (T_u(i)*sin(x_k(3,i)-a_1_u(i))+T_l(i)*sin(x_k(3,i)-a_1_l(i)))/mass - q(i)*w(i);
    wdot(i) = g*cos(theta(i)) - CDS/mass*0.5*rho*w(i)*V(i) - (T_u(i)*cos(x_k(3,i)-a_1_u(i))+T_l(i)*cos(x_k(3,i)-a_1_l(i)))/mass + q(i)*u(i);
    qdot(i) = -(T_u(i)*sin(x_k(3,i)-a_1_u(i))*0.89+T_l(i)*sin(x_k(3,i)-a_1_l(i))*(0.89+0.77))/Iyy;
    thetadot(i) = q(i);

    lambda_0_u_dot(i) = (C_T_BEM_u(i)-C_T_Glau_u(i))/0.1; % Quasi-Dynamic Inflow

    % Euler Integration
    w(i+1) = 0;
    q(i+1) = 0;
    theta(i+1) = 0;
    % lambda_0_u(i+1) = lambda_0_u(i)+dt*labidot_u(i);

end

% PLOTTING
%%
figure(1)
plot(mu, rad2deg(x_k(1,:)), '-*', mu, rad2deg(x_k(2,:)), '-*', mu, rad2deg(x_k(3,:)), '-*'),grid;
% ylim([-5, 30])
legend('theta_0_u', 'theta_0_l', 'theta_c'); xlabel("Advance Ratio \mu [-]")

%%
figure(2)
plot(mu, lambda_i_u_mean, mu, lambda_i_l_mean,mu,lambda_c),grid
legend('lambda i u mean', 'lambda i l mean', 'lambda c'); xlabel("Advance Ratio \mu [-]")

%%
figure(3)
plot(mu, C_T_Glau_u, '-*', mu, C_T_Glau_l, mu, C_T_BEM_u, mu, C_T_BEM_l),grid
legend('CT Glau lower', 'CT Glau lower', 'CT BEM upper', 'CT BEM lower'); xlabel("Advance Ratio \mu [-]")

figure(4)
plot(V, C_T_Glau_u, '-*', V, C_T_Glau_l, V, C_T_BEM_u, V, C_T_BEM_l),grid
legend('CT Glau lower', 'CT Glau lower', 'CT BEM upper', 'CT BEM lower'); xlabel("Advance Ratio \mu [-]")

%%
figure(5)
plot(mu, T_u/g, mu, T_l/g),grid
legend('T upper [kg]', 'T lower [kg]'); xlabel("Advance Ratio \mu [-]")

figure(6)
plot(V, T_u/g, V, T_l/g),grid
legend('T upper [kg]', 'T lower [kg]'); xlabel("Speed V [-]")

% IMPORTANT --- SPEED AND ADVANCE RATIO GIVE DIFFERENT RESULTS FOR LOWER
% ROTOR????

%%
figure(7)
plot(V, rad2deg(a_1_u))
legend("a 1"); xlabel("Advance Ratio \mu [-]")

% vi_lst = sqrt(-V_lst.^2/2 + sqrt(V_lst.^4/4+1));



    
    
    
    
    