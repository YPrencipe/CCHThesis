function [f_xki_ti, other_output] = f_xk(vel, x_k, q, theta_f)

    % Load helicopter parameters
    coaxial_heli_parameters;
    sigma = N_u*c*R/(pi*R^2);       % for some reason it doesnt want to read sigma from the parameters file, perhaps because of the variable name 'sigma' which is a greek letter

    %% Define variables for readability
    u = vel(1); w = vel(2);
    
    V = sqrt(u^2+w^2);
    p = 0; r = 0; v = 0;
    theta_0 = x_k(1); theta_c = x_k(2); theta_p = x_k(3);
    lambda_0_u = x_k(4); lambda_0_l = x_k(5); lambda_0_p = x_k(6);

    %% Rotor Model
    % Rotor RPM Scheduling
    if V > 70
        Omega = 40 - 5/30*(V-70);
        % Omega = 40;
    else
        Omega = 40;
    end

    alfa_sp = atan2(w,u);
    alfa_cp = alfa_sp + x_k(2);

    mu = V/(Omega*R);
    mu_x = V/(Omega*R) * cos(alfa_cp);
    mu_z = V/(Omega*R) * sin(alfa_cp);

    lambda_u = (mu_z - lambda_0_u)/2;
    lambda_l = (mu_z - lambda_0_l)/2;
    a0_u = lock/(8*nu_b2) * (theta_0*(1+mu_x^2) + 4/3*lambda_u + twist_r*(4/5+2/3*mu_x) - 4/3*mu_x*theta_c );
    a0_l = lock/(8*nu_b2) * (theta_0*(1+mu_x^2) + 4/3*lambda_l + twist_r*(4/5+2/3*mu_x) - 4/3*mu_x*theta_c );
    a1_u = (8/3*mu_x*theta_0 + 2*mu_x*lambda_u - 16/lock*q/Omega + 2*twist_r*mu_x - (1+3/2*mu_x^2)*theta_c) / (1-1/2*mu_x^2); 
    a1_l = (8/3*mu_x*theta_0 + 2*mu_x*lambda_l - 16/lock*q/Omega + 2*twist_r*mu_x - (1+3/2*mu_x^2)*theta_c) / (1-1/2*mu_x^2); 
    
    % H-forces
    CHdp_u = sigma*0.015*mu_x/4 + sigma*c_l_a/4*((a1_u*mu_x^2/2 + mu_x*lambda_u)*theta_0 ...
        + q/Omega*(-a0_u/3) + (a0_u^2+a1_u^2)*mu_x/2);
    CHdp_l = sigma*0.015*mu_x/4 + sigma*c_l_a/4*((a1_l*mu_x^2/2 + mu_x*lambda_l)*theta_0 ...
        + q/Omega*(-a0_l/3) + (a0_l^2+a1_l^2)*mu_x/2);

    Hdp_u = CHdp_u * rho * (Omega*R)^2 * pi*R^2;
    Hdp_l = CHdp_l * rho * (Omega*R)^2 * pi*R^2;

    % Rotor Angles
    a1r_u = -theta_c + a1_u;
    a1r_l = -theta_c + a1_l;

    % Thrust coefficients
        % Upper Rotor
        C_T_Glau_u = 2*lambda_0_u * sqrt( (mu*cos(alfa_cp))^2 + (mu*sin(alfa_cp)-lambda_0_u)^2 );       %% NOTE THE 1 IN THE FRONT OF THE EQ.!!
        C_T_BEM_u = sigma*c_l_a/2*( (1/3+mu_x^2/2)*theta_0 + (1+mu_x^2)/8*twist_r + lambda_u/2 );
        
        % Lower Rotor
        C_T_Glau_l = 2*lambda_0_l * sqrt( (mu*cos(alfa_cp))^2 + (mu*sin(alfa_cp)-lambda_0_l+lambda_0_u)^2 );
        C_T_BEM_l = sigma*c_l_a/2*( (1/3+mu_x^2/2)*theta_0 + (1+mu_x^2)/8*twist_r + lambda_l/2 );

    % Thrust using BEM
    T_u = C_T_BEM_u * rho * (Omega*R)^2 * pi*R^2;
    T_l = C_T_BEM_l * rho * (Omega*R)^2 * pi*R^2;
    % T_l = T_u;

    %% Pusher Propeller Model
    alfa_sp_p = atan2(w,u);

    mu_xp = sqrt(v^2 + (w + Kp*Omega*R*(x_k(4)+2*x_k(4)) + q*l_p)^2)/(omega_p*R_p);
    lambda_p = -(u-q*h_p-r*d_p)/(omega_p*R_p) - lambda_0_p;

    C_T_Glau_p = 2*lambda_0_p*sqrt( (mu*sin(alfa_sp_p))^2 + ...
        (mu*cos(alfa_sp_p)+lambda_0_p)^2 );       % this is like climbing flight for a rotor
    C_T_BEM_p = sigma_p*c_l_a_p/2*( (1/3+mu_xp^2/2)*theta_p + (1+mu_xp^2)/8*twist_p + lambda_p/2 );
    Tp = C_T_BEM_p * rho * (omega_p*R_p)^2 * pi*R_p^2;

    %% Horizontal tail model
    V_ht = sqrt(u^2 + (w+q*l_h)^2);
    alfa_ht = atan2((w+q*l_h), u) + alfa_ht_0;
    L_h = 1/2*rho*V_ht^2*S_ht*C_l_ht*alfa_ht; 

    %% Elevator Model

    %% Forces and Moments
    % X Forces
    X_MR_u_H = -Hdp_u*cos(a1r_u+gamma_s);
    X_MR_u_T = T_u*sin(a1r_u+gamma_s);
    X_MR_u = X_MR_u_H + X_MR_u_T;

    X_MR_l_H = -Hdp_l*cos(a1r_l+gamma_s);
    X_MR_l_T = T_l*sin(a1r_l+gamma_s);
    X_MR_l = X_MR_l_H + X_MR_l_T;

    X_MR = X_MR_u + X_MR_l;
    alfa_fus = atan2(w,u);
    X_fus = -0.5*rho*V^2*CDS*cos(alfa_fus);
    X_prop = Tp;
    X_ht = -L_h*sin(alfa_ht_0);
    X = X_MR + X_fus + X_prop + X_ht;

    % Y-forces
    Y = 0;

    % Z forces
    Z_MR_u_T = -T_u*cos(a1r_u+gamma_s);
    Z_MR_u_H = -Hdp_u*sin(a1r_u+gamma_s);
    Z_MR_u = Z_MR_u_T + Z_MR_u_H;

    Z_MR_l_T = -T_l*cos(a1r_l+gamma_s);
    Z_MR_l_H = -Hdp_l*sin(a1r_l+gamma_s);
    Z_MR_l = Z_MR_l_T + Z_MR_l_H;

    Z_MR = Z_MR_u + Z_MR_l;

    Z_fus = -0.5*rho*V^2*CDS*sin(alfa_fus);

    Z_ht = -L_h*cos(alfa_ht_0);

    Z = Z_MR + Z_fus + Z_ht;

    % L Moment
    L = 0;

    % M Moment
    M_MR = -h_u*X_MR_u - h_l*X_MR_l;
    M_hinge = (Omega*R)^2*hinge_offset_ratio*m_bl*sin(a1r_u+gamma_s);    % moment due to hinge offset
    M_flap = -N_b*K_b/2*a1_u;
    M_ht = l_h*Z_ht;       %% CHECK SIGN FOR M_w
    C_M_fus = (V/(Omega*R))^2 * 1/(pi*R^2*R) * K_fus * V_fus * (alfa_fus - alfa_fus_m_0);
    M_fus = rho*pi*R^2*(Omega*R^2)*R*C_M_fus;

    % M = M_MR + M_flap + M_hinge + M_ht + M_fus;
    M = (N_u/2*K_b)*a1_u + (N_l/2*K_b)*a1_l + h_l*T_l*a1_l + h_u*T_u*a1_u - x_cg*(T_u+T_l) - gamma_s*(h_u*T_u + h_l*T_u) + M_fus + (l_h+x_cg)*Z_ht;
 
    % N Moment
    N = 0;

    %% Equations of Motion
    udot = -q*w + X/mass -g*sin(theta_f);
    wdot = q*u + Z/mass + g*cos(theta_f);
    qdot = M/Iyy;
    
    lambda_0_u_dot = (C_T_BEM_u-C_T_Glau_u)/0.1;
    lambda_0_l_dot = (C_T_BEM_l-C_T_Glau_l)/0.1;
    lambda_p_dot = (C_T_BEM_p-C_T_Glau_p)/0.1;

    f_xki_ti = [udot; wdot; qdot; lambda_0_u_dot; lambda_0_l_dot; lambda_p_dot];

    % Other Calcs
    a_sound = 343;
    V_tip_adv = Omega*R + u;
    M_tip_adv = V_tip_adv/a_sound;




    other_output = [M_tip_adv];

    % F = [X; Y; Z];
    % Sigma = [L; M; N];  % big Sigma, moments
end






