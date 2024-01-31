% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Filename:     f_xk
% % Project:      MSc Thesis
% % Supervisor:   M.D. Pavel
% % Author:       Ynias Prencipe 
% % Student Nr.:  4777158
% % 
% % Description:  Function to determine the state vector using the small 
% % increments of the trim vector during the Newton-Rhapson trim algorithm
% %
% % Assumptions:  - Uniform inflows on both upper and lower rotors
% %               - Upper rotor inflow is fully added to lower rotor inflow
% %               - No rotor interference from lower to upper
% %               - No wake contraction
% % 
% % Issues to be solved:  - CT lower drops with speed , and becomes negative
% %                       - Strange behaviour of inflow at low speeds
% %                       - BEM thrust coefficient (CT) in wrong plane?
% % 
% % To Add:       - Pusher propeller
% %               - Elevator + Horizontal Stabilizer / Stabilator
% %               ? Wake contraction effects
% %               ? Rotor interference factor in function of advance ratio
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function f_xki_ti = f_xk(vel, x_k, inflow, q, theta)
% 
%     % Load helicopter parameters
%     coaxial_heli_parameters;
%     sigma = N_u*c*R/(pi*R^2);       % for some reason it doesnt want to read sigma from the parameters file, perhaps because of the variable name 'sigma' which is a greek letter
% 
%     % Define variables for readability
%     lambda_c = inflow(1); %lambda_0_u = inflow(2); %lambda_i_u_mean = inflow(3); lambda_i_l_mean = inflow(4);
%     mu = vel(4); mu_x = vel(5); mu_z = vel(6); 
%     % lambda_i_u_mean = lambda_0_u;
%     % lambda_i_l_mean = 2*lambda_0_u;
% 
%     % Longitudinal flapping angle
%     a1_u = (8/3*mu_x*x_k(1) - 2*mu_x*(lambda_c+lambda_0_u) - 16/lock*q/Omega) / (1-1/2*mu_x^2);
%     a1_l = (8/3*mu_x*x_k(2) - 2*mu_x*(lambda_c+2*lambda_0_u) - 16/lock*q/Omega) / (1-1/2*mu_x^2);
%     % a1_l=a1_u;
%     % a1_u = x_k(3);
%     % a1_l = x_k(3);
% 
%     % Thrust coefficients
%     C_T_Glau_ui = lambda_0_u * sqrt(mu^2+(mu_z+lambda_0_u)^2);      % has to be changed, see paper, maybe to 2*?
%     C_T_Glau_li = 2*2*lambda_0_u * sqrt(mu^2+(mu_z-2*lambda_0_u+lambda_0_u)^2);      % has to be changed, see paper
%     C_T_BEM_ui = 1/4*c_l_a*sigma*(2/3*x_k(1)*(1+3/2*mu^2)-(lambda_c+lambda_0_u));
%     C_T_BEM_li = 1/4*c_l_a*sigma*(2/3*x_k(2)*(1+3/2*mu^2)-(lambda_c+lambda_0_u));
%     % C_T_BEM_ui = 1/4*c_l_a*sigma*(2/3*x_k(1)*(1+3/2*mu^2)+(mu_z-lambda_0_u)/2);
%     % C_T_BEM_li = 1/4*c_l_a*sigma*(2/3*x_k(2)*(1+3/2*mu^2)+(mu_z-2*lambda_0_u)/2);
% 
%     % Thrust using BEM 
%     T_u = C_T_BEM_ui * rho * (Omega*R)^2 * pi*R^2;
%     T_l = C_T_BEM_li * rho * (Omega*R)^2 * pi*R^2;
% 
%     % State Equations
%     udot = -g*sin(theta) - CDS/mass*0.5*rho*vel(1)*vel(3) + (T_u*sin(x_k(3)-a1_u)+T_l*sin(x_k(3)-a1_l))/mass - q*vel(2);
%     wdot = g*cos(theta) - CDS/mass*0.5*rho*vel(2)*vel(3) - (T_u*cos(x_k(3)-a1_u)+T_l*cos(x_k(3)-a1_l))/mass + q*vel(1);
%     qdot = -(T_u*sin(x_k(3)-a1_u)*0.89+T_l*sin(x_k(3)-a1_l)*(0.89+0.77))/Iyy;
%     thetadot = q;
%     lambda_0_u_dot = (C_T_BEM_ui-C_T_Glau_ui)/0.1;
%     % lambda_0_l_dot = (C_T_BEM_li-C_T_Glau_li)/0.1;
% 
%     f_xki_ti = [udot; wdot; qdot; lambda_0_u_dot];
% end








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

% function f_xki_ti = f_xk(vel, x_k, inflow, q)
% 
%     % Load helicopter parameters
%     coaxial_heli_parameters;
%     sigma = N_u*c*R/(pi*R^2);       % for some reason it doesnt want to read sigma from the parameters file, perhaps because of the variable name 'sigma' which is a greek letter
% 
%     % Define variables for readability
%     lambda_c = inflow(1); %lambda_0_u = inflow(2); %lambda_i_u_mean = inflow(3); lambda_i_l_mean = inflow(4);
%     u = vel(1); w = vel(2); V = vel(3); 
%     mu = vel(4); mu_x = vel(5); mu_z = vel(6); 
% 
%     theta_0 = x_k(1); theta_f = x_k(2); theta_c = x_k(3); 
%     lambda_0_u = lambda_0_u;
% 
%     % Longitudinal flapping angle
%     a1_u = (8/3*mu_x*x_k(1) - 2*mu_x*(lambda_c+lambda_0_u) - 16/lock*q/Omega) / (1-1/2*mu_x^2);
%     a1_l = (8/3*mu_x*x_k(1) - 2*mu_x*(lambda_c+2*lambda_0_u) - 16/lock*q/Omega) / (1-1/2*mu_x^2);
% 
%     % Thrust coefficients
%     lambda_i_l = ( 1 + (-3.81*mu + 1.45) ) * lambda_0_u;
%     C_T_Glau_ui = lambda_0_u * sqrt(mu_x^2+(mu_z+lambda_0_u)^2);      % has to be changed, see paper, maybe to 2*?
%     % C_T_Glau_li = 2*lambda_0_u * sqrt(mu^2+(mu_z-lambda_i_l+lambda_0_u)^2);      % has to be changed, see paper
%     C_T_BEM_ui = 1/4*c_l_a*sigma*(2/3*x_k(1)*(1+3/2*mu_x^2)-(lambda_c+lambda_0_u));
%     % C_T_BEM_li = 1/4*c_l_a*sigma*(2/3*x_k(1)*(1+3/2*mu^2)-(lambda_c+lambda_i_l));
% 
%     % Thrust using BEM 
%     T_u = C_T_BEM_ui * rho * (Omega*R)^2 * pi*R^2;
%     T_l = ( 1 + (-3.81*mu_x + 1.45) ) * T_u;
% 
%     % Pusher Propeller Model
%     % Assumptions: - rigid blades - no flapping - no interference between
%     % blades - uniform induced velocity - inertial forces ignored - 
%     % spanwise flow ignored - thrust force acting straight through C.G. 
%     % on X-axis - no blade twist YET
%     % Momentum theory for vertical flight as described in Johnson can be
%     % used for the inflow modelling, since the pusher propeller is in fact
%     % in 'vertical' flight w.r.t. its own plane
% 
%     % theta_p = 0;
%     % diffTp = 10;
%     % while diffTp > 1
%     %     theta_p = theta_p + deg2rad(0.01);
%     % 
%     %     Tp_des = 0.5 * 1/2 * rho * V^2 * CDS;   % prop accounts for 50% of drag produced by heli
%     %     CTp = Tp_des / (rho * (Omega*R)^2 * pi*R^2);
%     %     lambda_p = sqrt(CTp/2);
%     % 
%     %     v_n_p = (lambda_p + mu_x) * omega_p*R_p;     % simplified for trimmed fwd flight
%     %     v_tan_p = omega_p*r_p;      % simplified for trimmed forward flight
%     %     v_res_p = sqrt(v_n_p.^2+v_tan_p.^2);
%     %     phi_i_p = atan2(v_n_p, v_tan_p);
%     %     alfa_p = theta_p + atan2(v_n_p, v_tan_p);
%     % 
%     %     l_p = 1/2*rho*v_res_p.^2*c_p*C_l_p.*alfa_p;
%     %     d_p = 1/2*rho*v_res_p.^2*c_p*C_d_p.*alfa_p;
%     %     dTp = l_p.*cos(phi_i_p) - d_p.*sin(phi_i_p);
%     %     Tp = N_p * sum(dTp, "all");
%     % 
%     %     diffTp = Tp_des - Tp;
%     % end
% 
%     Tp_des = 0.5 * 1/2 * rho * V^2 * CDS;   % prop accounts for 50% of drag produced by heli
%     CTp_des = Tp_des / (rho * (Omega*R)^2 * pi*R^2);
%     lambda_p = sqrt(CTp_des/2);
% 
%     diff_p = 10;
%     theta_p = 0;
%     while diff_p > 0.1
%         theta_p = theta_p + deg2rad(0.01);
%         CTp = 1/2*sigma*C_l_p*(theta_p/3 - lambda_p/2);
%         T_p = CTp * (rho * (Omega*R)^2 * pi*R^2);
%         diff_p = Tp_des - T_p;
%     end
% 
% 
%     % State Equations
%     udot = -g*sin(theta_f) - CDS/mass*0.5*rho*V*cos(theta_f) + (T_u*sin(a1r_u+i_theta)+T_l*sin((a1r_l+i_theta)))/mass - q*w ;%+ T_p/mass;
%     wdot = g*cos(theta_f) - CDS/mass*0.5*rho*V*sin(-theta_f) - (T_u*cos(a1r_u+i_theta)+T_l*cos((a1r_l+i_theta)))/mass + q*u;
%     qdot = -(T_u*sin(a1r_u+i_theta)*0.89+T_l*sin((a1r_l+i_theta))*(0.89+0.77))/Iyy;
%     lambda_0_u_dot = (C_T_BEM_ui-C_T_Glau_ui)/0.1;
%     % lambda_0_l_dot = (C_T_BEM_li-C_T_Glau_li)/0.1;
% 
%     f_xki_ti = [udot; wdot; qdot; lambda_0_u_dot];
% end


%% TRY 3

function f_xki_ti = f_xk(vel, x_k, inflow, q, flapping, theta_f)

    % Load helicopter parameters
    coaxial_heli_parameters;
    sigma = N_u*c*R/(pi*R^2);       % for some reason it doesnt want to read sigma from the parameters file, perhaps because of the variable name 'sigma' which is a greek letter

    % Define variables for readability
    lambda_c = inflow(1); %lambda_0_u = inflow(2); %lambda_i_u_mean = inflow(3); lambda_i_l_mean = inflow(4);
    u = vel(1); w = vel(2); V = vel(3); 
    mu = vel(4); mu_x = vel(5); mu_z = vel(6); 
    a1_u = flapping(1); a1_l = flapping(2); a0_u = flapping(3); a0_l = flapping(4);
    p = 0; r = 0; v = 0;

    % Rotor RPM Scheduling
    if V < 70
        Omega = 40;
    else
        Omega = 40 - 5/30*(V-70);
    end

    % theta_f = deg2rad(0); 
    delta_e = 0; 

    theta_0 = x_k(1); theta_c = x_k(2); theta_p = x_k(3);
    lambda_0_u = x_k(4); lambda_0_p = x_k(5);

    lambda_u = (mu_z - lambda_0_u)/2;
    lambda_l = (mu_z - 2*lambda_0_u)/2;
    a0_u = lock/(8*nu_b2) * (theta_0/2*(1+mu_x^2) + 4/3*lambda_u + twist_r*(4/5+2/3*mu_x) - 4/3*mu_x*theta_c );
    a0_l = lock/(8*nu_b2) * (theta_0/2*(1+mu_x^2) + 4/3*lambda_l + twist_r*(4/5+2/3*mu_x) - 4/3*mu_x*theta_c );
    a1_u = (8/3*mu_x*theta_0/2 + 2*mu_x*lambda_u - 16/lock*q/Omega + 2*twist_r*mu_x - (1+3/2*mu_x^2)*theta_c) / (1-1/2*mu_x^2); 
    a1_l = (8/3*mu_x*theta_0/2 + 2*mu_x*lambda_l - 16/lock*q/Omega + 2*twist_r*mu_x - (1+3/2*mu_x^2)*theta_c) / (1-1/2*mu_x^2); 
    
    % H-forces
    CHdp_u = sigma*0.015*mu_x/4 + sigma*c_l_a/4*((a1_u*mu_x^2/2 + mu_x*lambda_u)*theta_0/2 ...
        + q/Omega*(-a0_u/3) + (a0_u^2+a1_u^2)*mu_x/2);
    CHdp_l = sigma*0.015*mu_x/4 + sigma*c_l_a/4*((a1_l*mu_x^2/2 + mu_x*lambda_l)*theta_0/2 ...
        + q/Omega*(-a0_l/3) + (a0_l^2+a1_l^2)*mu_x/2);

    Hdp_u = CHdp_u * rho * (Omega*R)^2 * pi*R^2;
    Hdp_l = CHdp_l * rho * (Omega*R)^2 * pi*R^2;

    % Rotor Angles
    a1r_u = theta_c - a1_u;
    a1r_l = theta_c - a1_l;

    % Thrust coefficients
    C_T_Glau_ui = 2*lambda_0_u * sqrt( (mu*cos(a1r_u))^2 + (mu*sin(a1r_u)+lambda_0_u)^2 );
    % C_T_BEM_ui = 1/4*c_l_a*sigma*(2/3*theta_0 *(1+3/2*mu_x^2)-(lambda_c+lambda_0_u));
    C_T_BEM_ui = sigma*c_l_a/2*( (1/3+mu_x^2/2)*theta_0/2 + (1+mu_x^2)/8*twist_r + lambda_u/2 );
    C_T_BEM_li = sigma*c_l_a/2*( (1/3+mu_x^2/2)*theta_0/2 + (1+mu_x^2)/8*twist_r + lambda_l/2 );

    % Thrust using BEM 
    T_u = C_T_BEM_ui * rho * (Omega*R)^2 * pi*R^2;
    % T_l = C_T_BEM_li * rho * (Omega*R)^2 * pi*R^2;
    T_l = 2*T_u;

    % Pusher Propeller Model
    alfa_sp_p = atan2(w,u);

    mu_xp = sqrt(v^2 + (w + Kp*Omega*R*(x_k(4)+2*x_k(4)) + q*l_p)^2)/(omega_p*R_p);
    lambda_p = -(u-q*h_p-r*d_p)/(omega_p*R_p) - lambda_0_p;

    C_T_Glau_p = 2*lambda_0_p*sqrt( (mu*sin(alfa_sp_p))^2 + ...
        (mu*cos(alfa_sp_p)+lambda_0_p)^2 );       % this is like climbing flight for a rotor
    C_T_BEM_p = sigma_p*c_l_a_p/2*( (1/3+mu_xp^2/2)*theta_p + (1+mu_xp^2)/8*twist_p + lambda_p/2 );
    Tp = C_T_BEM_p * rho * (omega_p*R_p)^2 * pi*R_p^2;

    % All-moving horizontal tail (elevator)
    L_h = 1/2*rho*V^2*S_h*C_l_h*(delta_e+theta_f);
    L_h = 0;  

    %% Equations of Motion
    udot = -g*sin(theta_f) - 0.5*rho*u*V*CDS/mass ...
        + (T_u*sin(a1r_u+i_theta)+T_l*sin(a1r_l+i_theta))/mass - q*w + Tp/mass ...
        - (Hdp_u*cos(a1r_u+i_theta)+Hdp_l*cos(a1r_l+i_theta))/mass;
    wdot = g*cos(theta_f) - 0.5*rho*w*V*CDS/mass ...
        - (T_u*cos(a1r_u+i_theta)+T_l*cos(a1r_l+i_theta))/mass + q*u - L_h/mass ...
        - (Hdp_u*sin(a1r_u+i_theta)+Hdp_l*sin(a1r_l+i_theta))/mass;
    qdot = (-(T_u*sin(a1r_u+i_theta)*0.89+T_l*sin(a1r_l+i_theta)*(0.89+0.77)) - L_h*l_h - ...
        (T_u*cos(a1r_u+i_theta)*0.0+T_l*cos(a1r_l+i_theta)*0.0)) / Iyy ;
    
    lambda_0_u_dot = (C_T_BEM_ui-C_T_Glau_ui)/0.1;
    lambda_p_dot = (C_T_BEM_p-C_T_Glau_p)/0.1;

    f_xki_ti = [udot; wdot; qdot; lambda_0_u_dot; lambda_p_dot];
end





