%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename:     f_xk
% Project:      MSc Thesis
% Supervisor:   M.D. Pavel
% Author:       Ynias Prencipe 
% Student Nr.:  4777158
% 
% Description:  Function to determine the state vector using the small 
% increments of the trim vector during the Newton-Rhapson trim algorithm
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

function f_xki_ti = f_xk(vel, x_k, inflow, q, theta)

    % Load helicopter parameters
    coaxial_heli_parameters;
    sigma = N_u*c*R/(pi*R^2);       % for some reason it doesnt want to read sigma from the parameters file, perhaps because of the variable name 'sigma' which is a greek letter

    % Define variables for readability
    lambda_c = inflow(1); lambda_0_u = inflow(2); %lambda_i_u_mean = inflow(3); lambda_i_l_mean = inflow(4);
    mu = vel(4); mu_x = vel(5); mu_z = vel(6); 
    lambda_i_u_mean = x_k(4);
    lambda_i_l_mean = 2*x_k(4);

    % Longitudinal flapping angle
    a_1_u = (8/3*mu_x*x_k(1) - 2*mu_x*(lambda_c+lambda_i_u_mean) - 16/lock*q/Omega) / (1-1/2*mu_x^2);
    a_1_l = (8/3*mu_x*x_k(2) - 2*mu_x*(lambda_c+lambda_i_l_mean) - 16/lock*q/Omega) / (1-1/2*mu_x^2);

    % Thrust coefficients
    C_T_BEM_ui = 1/4*c_l_a*sigma*(2/3*x_k(1)*(1+3/2*mu^2)-(lambda_c+x_k(4)));
    C_T_BEM_li = 1/4*c_l_a*sigma*(2/3*x_k(2)*(1+3/2*mu^2)-(lambda_c+2*x_k(4)));
    C_T_Glau_ui = lambda_i_u_mean * sqrt(mu^2+(mu_z+lambda_i_u_mean)^2);      % has to be changed, see paper, maybe to 2*?
    C_T_Glau_li = 2*lambda_i_l_mean * sqrt(mu^2+(mu_z+lambda_i_l_mean)^2);      % has to be changed, see paper

    % Thrust using BEM 
    T_u = C_T_BEM_ui * rho * (Omega*R)^2 * pi*R^2;
    T_l = C_T_BEM_li * rho * (Omega*R)^2 * pi*R^2;

    % State Equations
    udot = -g*sin(theta) - CDS/mass*0.5*rho*vel(1)*vel(3) + (T_u*sin(x_k(3)-a_1_u)+T_l*sin(x_k(3)-a_1_l))/mass - q*vel(2);
    wdot = g*cos(theta) - CDS/mass*0.5*rho*vel(2)*vel(3) - (T_u*cos(x_k(3)-a_1_u)+T_l*cos(x_k(3)-a_1_l))/mass + q*vel(1);
    qdot = -(T_u*sin(x_k(3)-a_1_u)*0.89+T_l*sin(x_k(3)-a_1_l)*(0.89+0.77))/Iyy;
    thetadot = q;
    lambda_0_u_dot = (C_T_BEM_ui-C_T_Glau_ui)/0.1;
    % lambda_0_l_dot = (C_T_BEM_li-C_T_Glau_li)/0.1;

    f_xki_ti = [udot; wdot; qdot; lambda_0_u_dot];
end