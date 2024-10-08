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
Omega = 40;

% Initial Values
theta_0(1) = deg2rad(6.1);
theta_c(1) = deg2rad(0);
t(1)=0;
u(1)=0;
w(1)=0;
q(1)=0;
theta(1)=deg2rad(0.01);
x(1)=0;
lambda_0_u(1)=sqrt(mass/2*abs(g)/(area*2*rho))/(Omega*R);
z(1)=0;

% Time Paramaters
simN = 400;
tEnd = 15;
dt = (tEnd-t(1))/simN;
tau = 0.1;

z_des = 8; 
KaltP = 0.4;
% KaltD = 0.2;
Ki_c = 0.4;

for i = 1:simN

    % Cyclic step input of 1deg between t=0.5 and t=1s
    if t(i) >= 2 & t(i) <= 3 
       theta_c(i)=deg2rad(1);%longit cyclic input
    else 
       theta_c(i)=deg2rad(0);
    end
    
    % Control law for collective
    c(i)=u(i)*sin(theta(i))-w(i)*cos(theta(i));
    height(i)=-z(i);

    e_s(i) = z_des - height(i);

    theta_0(i)=theta_0(1);

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

    alfa_c(i) = theta_c(i) - phi(i);

    V(i) = sqrt(u(i)^2+w(i)^2);

    mu(i) = V(i)/(Omega*R);
    lambda_c(i) = V(i)/(Omega*R) * sin(alfa_c(i));
    mu_x(i) = V(i)/(Omega*R) * cos(alfa_c(i));
    mu_z(i) = V(i)/(Omega*R) * sin(alfa_c(i));

    %%
    % [lambdai_u_mean, lambdai_l_mean] = calculateRotorInflow(mu(i), lambda_0_u(i), mu_x(i), mu_z(i));
    % lambda_i_u_mean(i) = lambdai_u_mean;
    % lambda_i_l_mean(i) = lambdai_l_mean;
    
    % Define Polar Coordinates
    psi_values = 0:deg2rad(1):2*pi;
    r_values = 0:0.01:R;
    [r, psi] = meshgrid(r_values, psi_values);
    % Include Uniform (upper)/Linear (Lower) Inflow Distributions
    ksi(i) = atan(mu_x(i)/(mu_z(i)+lambda_0_u(i)));
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

    lambda_i_u_prime = lambda_0_u(i)*(1 + 0*r.*cos(psi) + 0*r.*sin(psi));
    lambda_i_l_prime = lambda_0_u(i)*(1 + k_x*r.*cos(psi) + k_y*r.*sin(psi));

    lambda_i_u = lambda_i_u_prime + lambda_i_l_prime*delta_l2u;
    lambda_i_l = lambda_i_l_prime + lambda_i_u_prime*delta_u2l;

    lambda_i_u_mean(i) = mean(lambda_i_u, "all");
    lambda_i_l_mean(i) = mean(lambda_i_l, "all");

    %%
    % Flapping Motion
    a_1_u(i) = (8/3*mu(i)*theta_0(i) - 2*mu(i)*(lambda_c(i)+lambda_i_u_mean(i)) - 16/lock*q(i)/Omega) / (1-1/2*mu(i)^2);
    a_1_l(i) = (8/3*mu(i)*theta_0(i) - 2*mu(i)*(lambda_c(i)+lambda_i_l_mean(i)) - 16/lock*q(i)/Omega) / (1-1/2*mu(i)^2);
 
    a1r_u(i) = theta_c(i) - a_1_u(i);
    a1r_l(i) = theta_c(i) - a_1_l(i);

    % Inflow and Thrust Calculation
    C_T_BEM_u(i) = 1/4*c_l_a*sigma*(2/3*theta_0(i)*(1+3/2*mu(i)^2)-(mu_z(i)+lambda_0_u(i)));
    C_T_BEM_l(i) = 1/4*c_l_a*sigma*(2/3*theta_0(i)*(1+3/2*mu(i)^2)-(mu_z(i)+lambda_0_u(i)));

    %%
    C_T_Gl_u(i) = lambda_i_u_mean(i) * sqrt(mu_x(i)^2+(mu_z(i)+lambda_i_u_mean(i))^2);

    T_u(i) = C_T_BEM_u(i) * rho * (Omega*R)^2 * pi*R^2;
    T_l(i) = C_T_BEM_l(i) * rho * (Omega*R)^2 * pi*R^2;

    %%% Horizontal Stabilizer
    alpha_ht = -360:360;
    Cl_ht = zeros(size(alpha_ht));
    range_indices = (alpha_ht >= -20 & alpha_ht <= 20);
    Cl_ht(range_indices) = interp1([-20, 20], [-1.2, 1.2], alpha_ht(range_indices), 'linear');
    alfa_ht(i) = atan2((w(i)+q(i)*l_h), u(i)) + alfa_ht_0;
    Cl_ht_interp(i) = interp1(alpha_ht, Cl_ht, rad2deg(alfa_ht(i)), 'linear', 'extrap');

    V_ht(i) = sqrt(u(i)^2 + (w(i)+q(i)*l_h)^2);
    L_h(i) = 0.5*rho*V_ht(i)^2*S_ht*Cl_ht_interp(i);

    %%% Fuselage
    alfa_fus(i) = atan2(w(i),u(i));
    C_M_fus(i) = (V(i)/(Omega*R))^2 * 1/(pi*R^2*R) * K_fus * Vol_fus * (alfa_fus(i) - alfa_fus_m_0);

    %% Forces and Moments
    X_u(i) = T_u(i)*sin(a1r_u(i));
    X_l(i) = T_l(i)*sin(a1r_l(i));
    X_fus(i) = -CDS*0.5*rho*u(i)*V(i);

    Z_u(i) = -T_u(i)*cos(a1r_u(i));
    Z_l(i) = -T_l(i)*cos(a1r_l(i));
    Z_fus(i) = -CDS/mass*0.5*rho*w(i)*V(i);
    Z_ht(i) = -L_h(i);  

    M_u(i) = -X_u(i)*h_u;
    M_l(i) = -X_l(i)*h_l;
    M_fus(i) = rho*pi*R^2*(Omega*R^2)*R*C_M_fus(i);
    M_ht(i) = l_h*Z_ht(i);


    X(i) = X_u(i) + X_l(i) + X_fus(i);
    Z(i) = Z_u(i) + Z_l(i) + Z_fus(i) + Z_ht(i);
    M(i) = M_u(i) + M_l(i) + M_ht(i) + M_fus(i);

    % Equations of Motion
    udot(i) = -g*sin(theta(i)) + X(i)/mass - q(i)*w(i);
    wdot(i) = g*cos(theta(i)) + Z(i)/mass - q(i)*u(i);
    qdot(i) = M(i)/Iyy;
    thetadot(i) = q(i);
    xdot(i) = u(i)*cos(theta(i))+w(i)*sin(theta(i));
    zdot(i) = -c(i);
    labidot_u(i) = (C_T_BEM_u(i)-C_T_Gl_u(i))/tau;


    % Euler Integration
    u(i+1) = u(i)+dt*udot(i);
    w(i+1) = w(i)+dt*wdot(i);
    q(i+1) = q(i)+dt*qdot(i);
    theta(i+1) = theta(i)+dt*thetadot(i);

    lambda_0_u(i+1) = lambda_0_u(i)+dt*labidot_u(i);

    x(i+1) = x(i)+dt*xdot(i);
    z(i+1) = z(i)+dt*zdot(i);

    t(i+1) = t(i)+dt;

end

%%
figure(1)
subplot(3,1,1)
plot(t(1:end-1),rad2deg(theta_c)), grid, legend('theta c');
subplot(3,1,2)
plot(t,rad2deg(theta),t,rad2deg(q)), grid, legend('theta f', 'q');
subplot(3,1,3)
plot(t,u,t(1:end-1),udot,t,w), grid, legend('u', 'udot','w');
figure(2)
plot(t,x,t(1:end-1),height); grid on; legend('x', 'height')

figure(2)
plot(t(1:end-1),M_ht,t(1:end-1),M_u,t(1:end-1),M_l,t(1:end-1),M_fus); 
grid on; legend('M_ht','M_u','M_l','M_{fus}')

figure(5)
plot(t(1:end-1),rad2deg(a_1_u)); hold on; plot(t(1:end-1), rad2deg(a_1_l)); 
plot(t(1:end-1), rad2deg(theta_c)); plot(t(1:end-1), rad2deg(alfa_ht));
plot(t(1:end-1), rad2deg(alfa_c)); plot(t(1:end-1), rad2deg(phi)); hold off; grid on;
legend('a_1_u','a_1_l','theta_c','alfa_ht','alfa_c','phi')

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

