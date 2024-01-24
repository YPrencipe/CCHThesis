% Filename: coaxial_model_3dof
% Course: Thesis
% Supervisor: M.D. Pavel
% Author: Ynias Prencipe 
% Student Number: 4777158
% Date of Delivery:
% Description: 
% Assumptions: 

clear; clc; close all;

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
x(1)=0;
lambda_0_u(1)=sqrt(mass*abs(g)/(area*2*rho))/(Omega*R);
z(1)=0;

Vmax = 100;

% Time Paramaters
simN = 400;
tEnd = 40;
dt = (tEnd-t(1))/simN;
tau = 0.1;

u_vals = [0.1, 0.2, 0.3,0.4,0.5,0.6, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80];

x_k = zeros(4,length(u_vals));
x_k(:,1) = [theta_0_u(1); theta_0_l(1); theta_c(1); lambda_0_u(1)];

for i = 1:length(u_vals)
    i
    u(i) = u_vals(i);
    w(i) = 0;        % w set to 0 during trim procedure
    theta(i) = 0;    % fix theta_fuselage to 0 for all flightspeeds for trimming purposes

    if i>1
        x_k(:,i) = x_k(:,i-1);
    end
   
    % Calculate Parameters
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

    V(i) = sqrt(u(i)^2+w(i)^2);

    mu(i) = V(i)/(Omega*R) * cos(alfa_c(i));
    lambda_c(i) = V(i)/(Omega*R) * sin(alfa_c(i));
    mu_x(i) = u(i)/(Omega*R);
    mu_z(i) = w(i)/(Omega*R);

    vel(:,i) = [u(i); w(i); V(i); mu(i); mu_x(i); mu_z(i)];

    %%
    % Define Polar Coordinates
    psi_values = 0:deg2rad(1):2*pi;
    r_values = 0:0.01:R;
    [r, psi] = meshgrid(r_values, psi_values);
    % Include Uniform (upper)/Linear (Lower) Inflow Distributions
    ksi(i) = atan(mu_x(i)/(mu_z(i)+x_k(4,i)));
    k_x = 15*pi/32 * tan(ksi(i)/2);        % Pitt and Peters estimates for k_x and k_y
    k_y = 0;

    % Interference Effects
    if (0<=mu(i)) && (mu(i)<= 0.316)
        delta_l2u = -2.15*mu(i) + 0.68;
    else
        % delta_l2u = -2.15*mu(i) + 0.68;     % this is not correct but 0 gives issues
        delta_l2u = 0;
    end
    if (0<=mu(i)) && (mu(i)<=0.381)
        delta_u2l = -3.81*mu(i) + 1.45;    
    else
        % delta_u2l = -3.81*mu(i) + 1.45;      % not correct but 0 gives issues
        delta_u2l = 0;                         % 0 gives jumps in C_T and a_1
    end

    % Inflow interferences and linear distribution n stuff
    lambda_i_u_prime = x_k(4,i)*(1 + 0*r.*cos(psi) + 0*r.*sin(psi));
    lambda_i_l_prime = x_k(4,i)*(1 + k_x*r.*cos(psi) + k_y*r.*sin(psi));

    lambda_i_u = lambda_i_u_prime; %+ lambda_i_l_prime*delta_l2u;
    lambda_i_l = lambda_i_l_prime + lambda_i_u_prime*1; %+ lambda_i_u_prime*delta_u2l;

    % lambda_i_u_mean(i) = mean(lambda_i_u, "all");
    % lambda_i_l_mean(i) = mean(lambda_i_l, "all");
    lambda_i_u_mean(i) = x_k(4,i);
    lambda_i_l_mean(i) = 2*x_k(4,i);

    inflow(:,i) = [lambda_c(i); x_k(4,i); lambda_i_u_mean(i); lambda_i_l_mean(i)];

    %%
    count(i) = 1;
    f_sum(1) = 10;

    while f_sum(count(i)) > 0.001

        count(i) = count(i)+1;
        tistep = deg2rad(0.1);

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

    lambda_i_u_mean(i) = x_k(4,i);
    lambda_i_l_mean(i) = 2*x_k(4,i);

    C_T_Glau_u(i) = lambda_i_u_mean(i) * sqrt(mu(i)^2+(mu_z(i)+lambda_i_u_mean(i))^2);
    C_T_Glau_l(i) = 2*lambda_i_l_mean(i) * sqrt(mu(i)^2+(mu_z(i)-lambda_i_l_mean(i)+lambda_i_u_mean(i))^2);
    C_T_BEM_u(i) = 1/4*c_l_a*sigma*(2/3*x_k(1,i)*(1+3/2*mu(i)^2)-(lambda_c(i)+x_k(4,i)));
    C_T_BEM_l(i) = 1/4*c_l_a*sigma*(2/3*x_k(2,i)*(1+3/2*mu(i)^2)-(lambda_c(i)+2*x_k(4,i)));

    T_u(i) = C_T_BEM_u(i) * rho * (Omega*R)^2 * pi*R^2;
    T_l(i) = C_T_BEM_l(i) * rho * (Omega*R)^2 * pi*R^2;

    a_1_u(i) = (8/3*mu_x(i)*x_k(1,i) - 2*mu_x(i)*(lambda_c(i)+lambda_i_u_mean(i)) - 16/lock*q(i)/Omega) / (1-1/2*mu_x(i)^2);
    a_1_l(i) = (8/3*mu_x(i)*x_k(2,i) - 2*mu_x(i)*(lambda_c(i)+lambda_i_l_mean(i)) - 16/lock*q(i)/Omega) / (1-1/2*mu_x(i)^2);

    udot(i) = -g*sin(theta(i)) - CDS/mass*0.5*rho*u(i)*V(i) + (T_u(i)*sin(x_k(3,i)-a_1_u(i))+T_l(i)*sin(x_k(3,i)-a_1_l(i)))/mass - q(i)*w(i);
    wdot(i) = g*cos(theta(i)) - CDS/mass*0.5*rho*w(i)*V(i) - (T_u(i)*cos(x_k(3,i)-a_1_u(i))+T_l(i)*cos(x_k(3,i)-a_1_l(i)))/mass + q(i)*u(i);
    qdot(i) = -(T_u(i)*sin(x_k(3,i)-a_1_u(i))*0.89+T_l(i)*sin(x_k(3,i)-a_1_l(i))*(0.89+0.77))/Iyy;
    thetadot(i) = q(i);

    % Quasi-Dynamic Inflow
    lambda_0_u_dot(i) = (C_T_BEM_u(i)-C_T_Glau_u(i))/0.1;
    

    % Euler Integration
    w(i+1) = 0;
    q(i+1) = 0;
    theta(i+1) = 0;
    % lambda_0_u(i+1) = lambda_0_u(i)+dt*labidot_u(i);



end

%%
figure(1)
plot(V, rad2deg(x_k(1,:)), '-*', V, rad2deg(x_k(2,:)), '-*', V, rad2deg(x_k(3,:)), '-*'),grid;
% ylim([-5, 30])
legend('theta_0_u', 'theta_0_l', 'theta_c')

%%
figure(2)
plot(V, lambda_i_u_mean, V, lambda_i_l_mean,V,lambda_c),grid
legend('lambda i u mean', 'lambda i l mean', 'lambda c')

%%
figure(3)
plot(V, C_T_Glau_u, '-*', V, C_T_Glau_l, V, C_T_BEM_u, V, C_T_BEM_l),grid
legend('CT Glau lower', 'CT Glau lower', 'CT BEM upper', 'CT BEM lower')

%%
figure(4)
plot(V, T_u/g, V, T_l/g),grid
legend('T upper [kg]', 'T lower [kg]')

%%
figure(5)
plot(V, lambda_0_u_dot)


    
    
    
    
    