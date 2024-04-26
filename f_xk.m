function [f_xki_ti, M_components, Mw_components, angles, inflow] = f_xk(vel, x_k, q, theta_f, delta_e)

    % Load helicopter parameters
    coaxial_heli_parameters;
    sigma = N_u*c*R/(pi*R^2);       % for some reason it doesnt want to read sigma from the parameters file, perhaps because of the variable name 'sigma' which is a greek letter

    %% Define variables for readability
    u = vel(1); w = vel(2);
    V = sqrt(u^2+w^2);
    p = 0; r = 0; v = 0;
    theta_0 = x_k(1); theta_c = x_k(2); theta_p = x_k(3); 
    lambda_0_u = x_k(4); lambda_0_l = x_k(5); lambda_0_p = x_k(6);
    alfa_fus = atan2(w,u);

    %% Rotor Model
    % Rotor RPM Scheduling
    if V > 70
        % Omega = 40 - 5/30*(V-70);
        Omega = 40;
    else
        Omega = 40;
    end

    alfa_sp = atan2(w,u);
    alfa_cp = alfa_sp + theta_c;

    mu = V/(Omega*R);
    mu_x = V/(Omega*R) * cos(alfa_cp);
    mu_z = V/(Omega*R) * sin(alfa_cp);

    %%%%%%% INFLOW
    % Interference Factors
    if (0<=mu) && (mu<= 0.316)
        delta_l2u = -2.15*mu + 0.68;
    else
        % delta_l2u = -2.15*mu + 0.68;     % this is not correct but 0 gives issues
        delta_l2u = 0;
    end
    if (0<=mu) && (mu<=0.381)
        delta_u2l = -3.81*mu + 1.45;    
    else
        % delta_u2l = -3.81*mu + 1.45;      % not correct but 0 gives issues
        delta_u2l = 0;                         % 0 gives jumps in C_T and a_1
    end
    
    % lambda_0_u = lambda_0_u + lambda_0_l*delta_l2u;
    % lambda_0_l = lambda_0_l + lambda_0_u*delta_u2l;

    lambda_u = (mu_z-lambda_0_u)/2;
    lambda_l = (mu_z-lambda_0_l)/2;

    %%%%%%% FLAPPING ANGLES
    a0_u = lock/(8*nu_b2) * (theta_0*(1+mu_x^2) + 4/3*lambda_u + twist_r*(4/5+2/3*mu_x) - 4/3*mu_x*theta_c );
    a0_l = lock/(8*nu_b2) * (theta_0*(1+mu_x^2) + 4/3*lambda_l + twist_r*(4/5+2/3*mu_x) - 4/3*mu_x*theta_c );
    a1_u = (8/3*mu_x*theta_0 + 2*mu_x*lambda_u - 16/lock*q/Omega + 2*twist_r*mu_x - (1+3/2*mu_x^2)*theta_c) / (1-1/2*mu_x^2); 
    a1_l = (8/3*mu_x*theta_0 + 2*mu_x*lambda_l - 16/lock*q/Omega + 2*twist_r*mu_x - (1+3/2*mu_x^2)*theta_c) / (1-1/2*mu_x^2); 

    %%%%%%% H-FORCES
    CHdp_u = sigma*0.015*mu_x/4 + sigma*c_l_a/4*((a1_u*mu_x^2/2 + mu_x*lambda_u)*theta_0 ...
        + q/Omega*(-a0_u/3) + (a0_u^2+a1_u^2)*mu_x/2);
    CHdp_l = sigma*0.015*mu_x/4 + sigma*c_l_a/4*((a1_l*mu_x^2/2 + mu_x*lambda_l)*theta_0 ...
        + q/Omega*(-a0_l/3) + (a0_l^2+a1_l^2)*mu_x/2);

    Hdp_u = CHdp_u * rho * (Omega*R)^2 * pi*R^2;
    Hdp_l = CHdp_l * rho * (Omega*R)^2 * pi*R^2;

    % Rotor Angles
    a1r_u = theta_c - a1_u;
    a1r_l = theta_c - a1_l;

    alfa_dp_u = alfa_cp + a1_u;
    alfa_dp_l = alfa_cp + a1_l;

    % Thrust coefficients
        % Upper Rotor
        C_T_Glau_u = 2*lambda_0_u * sqrt( (mu*cos(alfa_dp_u))^2 + (mu*sin(alfa_dp_u) + lambda_0_u)^2 );       %% NOTE THE 1 IN THE FRONT OF THE EQ.!!
        C_T_BEM_u = sigma*c_l_a/2*( (1/3+mu_x^2/2)*theta_0 + (1+mu_x^2)/8*twist_r + lambda_u/2 );

        % Lower Rotor
        C_T_Glau_l = 2*lambda_0_l * sqrt( (mu*cos(alfa_dp_l))^2 + (mu*sin(alfa_dp_l) + lambda_0_l + lambda_0_u)^2 );
        C_T_BEM_l = sigma*c_l_a/2*( (1/3+mu_x^2/2)*theta_0 + (1+mu_x^2)/8*twist_r + lambda_l/2 );

    % Thrust using BEM
    T_u = C_T_BEM_u * rho * (Omega*R)^2 * pi*R^2;
    T_l = C_T_BEM_l * rho * (Omega*R)^2 * pi*R^2;

    %% Pusher Propeller Model
    alfa_sp_p = atan2(w,u);

    mu_xp = sqrt(v^2 + (w + Kp*Omega*R*(lambda_0_u+2*lambda_0_u) + q*l_p)^2)/(omega_p*R_p);
    lambda_p = -(u-q*h_p-r*d_p)/(omega_p*R_p) - lambda_0_p;

    C_T_Glau_p = 2*lambda_0_p*sqrt( (mu*sin(alfa_sp_p))^2 + ...
        (mu*cos(alfa_sp_p)+lambda_0_p)^2 );       % this is like climbing flight for a rotor
    C_T_BEM_p = sigma_p*c_l_a_p/2*( (1/3+mu_xp^2/2)*theta_p + (1+mu_xp^2)/8*twist_p + lambda_p/2 );
    Tp = C_T_BEM_p * rho * (omega_p*R_p)^2 * pi*R_p^2;

    %% HORIZONTAL STABILIZER MODEL
    alpha_ht = -360:360;
    Cl_ht = zeros(size(alpha_ht));
    range_indices = (alpha_ht >= -20 & alpha_ht <= 20);
    Cl_ht(range_indices) = interp1([-20, 20], [-1.2, 1.2], alpha_ht(range_indices), 'linear');
    alfa_ht = atan2((w+q*l_h), u) + alfa_ht_0;
    
    Cl_ht_interp = interp1(alpha_ht, Cl_ht, rad2deg(alfa_ht), 'linear', 'extrap');
    V_ht = sqrt(u^2 + (w+q*l_h)^2);
    L_h = 0.5*rho*V_ht^2*S_ht*(Cl_ht_interp+(0.7*rad2deg(delta_e)));           % 0.859 comes from george's thesis

    %% ELEVATOR MODEL
    alpha_e = -360:360;
    Cl_e = zeros(size(alpha_e));
    range_indices = (alpha_e >= -20 & alpha_e <= 20);
    Cl_e(range_indices) = interp1([-20, 20], [-1.2, 1.2], alpha_e(range_indices), 'linear');
    alfa_e = alfa_ht + delta_e;

    Cl_e_interp = interp1(alpha_e, Cl_e, rad2deg(alfa_e), 'linear', 'extrap');
    V_e = V_ht;
    L_e = 0.5*rho*V_e^2*S_e*Cl_e_interp;
    % L_e = 0;

    %% Forces and Moments
    %%%%% X FORCE COMPONENTS
    X_MR_u = T_u*sin(a1r_u+gamma_s) - Hdp_u*cos(a1r_u+gamma_s);
    X_MR_l = T_l*sin(a1r_l+gamma_s) - Hdp_l*cos(a1r_l+gamma_s);
    X_MR = X_MR_u + X_MR_l;
    X_fus = -0.5*rho*V^2*CDS*cos(alfa_fus);
    X_prop = Tp;
    X_ht = -L_h*sin(alfa_ht_0);     % alfa_ht_0 is correct!! look at the correct reference frame
    X_e = -L_e*sin(alfa_ht_0+delta_e);

    %%%%% Z FORCE COMPONENTS
    Z_MR_u = -T_u*cos(a1r_u+gamma_s) + Hdp_u*sin(a1r_u+gamma_s);
    Z_MR_l = -T_l*cos(a1r_l+gamma_s) + Hdp_l*sin(a1r_l+gamma_s);
    Z_MR = Z_MR_u + Z_MR_l;
    Z_fus = -0.5*rho*V^2*CDS*sin(alfa_fus);
    Z_ht = -L_h*cos(alfa_ht_0);
    Z_e = -L_e*cos(alfa_ht_0+delta_e);

    %%%%% M MOMENT COMPONENTS
    M_MR_u = -h_u*X_MR_u - Z_MR_u*x_MR;
    M_MR_l = -h_l*X_MR_l - Z_MR_l*x_MR;
    M_hinge_u = -(Omega)^2*R*hinge_offset_ratio*R*m_bl*sin(a1r_u+gamma_s);    % moment due to hinge offset
    M_hinge_l = -(-Omega)^2*R*hinge_offset_ratio*R*m_bl*sin(a1r_l+gamma_s);    % moment due to hinge offset
    M_hinge = M_hinge_u + M_hinge_l;
    % M_hinge = 0;

    M_flap = N_u*K_b/2*a1_u + N_l*K_b/2*a1_l;
    % M_flap = 0;

    C_M_fus = (V/(Omega*R))^2 * 1/(pi*R^2*R) * K_fus * Vol_fus * (alfa_fus - alfa_fus_m_0);
    M_fus = rho*pi*R^2*(Omega*R^2)*R*C_M_fus;
    M_ht = (l_h+x_cg)*Z_ht;
    M_e = (l_h+x_cg)*Z_e;

    %%%%% TOTAL FORCES AND MOMENTS
    Forces.X = X_MR + X_fus + X_prop + X_ht + X_e;
    Forces.Z = Z_MR + Z_fus + Z_ht + Z_e;
    Moments.M = M_MR_u + M_MR_l + M_fus + M_ht + M_hinge + M_flap + M_e;


    %% Equations of Motion
    udot = Forces.X/mass - g*sin(theta_f) - q*w ;
    wdot = q*u + Forces.Z/mass + g*cos(theta_f);
    qdot = Moments.M/Iyy;
        qdot_MR_u = M_MR_u/Iyy;
        qdot_MR_l = M_MR_l/Iyy;
        qdot_hinge = M_hinge/Iyy;
        qdot_flap = M_flap/Iyy;
        qdot_fus = M_fus/Iyy;
        qdot_ht = M_ht/Iyy;
        qdot_e = M_e/Iyy;


    lambda_0_u_dot = (C_T_BEM_u-C_T_Glau_u)/0.1;
    lambda_0_l_dot = (C_T_BEM_l-C_T_Glau_l)/0.1;
    lambda_p_dot = (C_T_BEM_p-C_T_Glau_p)/0.1;

    f_xki_ti = [udot; wdot; qdot; lambda_0_u_dot; lambda_0_l_dot; lambda_p_dot];
    M_components = [M_MR_u, M_MR_l, M_hinge, M_flap, M_fus, M_ht, M_e];
    Mw_components = [qdot_MR_u, qdot_MR_l, qdot_hinge, qdot_flap, qdot_fus, qdot_ht, qdot_e];
    angles = rad2deg([a0_u, a0_l, a1_u, a1_l, a1r_u, a1r_l, alfa_sp, alfa_cp, alfa_dp_u, alfa_dp_l]);
    inflow = [lambda_0_u, lambda_0_l];

    
end




