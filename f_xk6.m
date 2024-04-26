function [f_xki_ti, M_components, Mw_components, angles, inflow, latAngles] = ...
    f_xk6(vel, x_k, p, q, r, theta_f, delta_e, delta_r, theta_cdiff)

    % Load helicopter parameters
    coaxial_heli_parameters;
    sigma = N_u*c*R/(pi*R^2);       % for some reason it doesnt want to read sigma from the parameters file, perhaps because of the variable name 'sigma' which is a greek letter

    %% Define variables for readability
    u = vel(1); v = vel(2); w = vel(3);
    V = sqrt(u^2+v^2+w^2);

    theta_0 = x_k(1); theta_d = x_k(2); theta_s=x_k(3);
    theta_c = x_k(4); phi_f = x_k(5); theta_p = x_k(6); 
    lambda_0_u = x_k(7); lambda_0_l = x_k(8); lambda_0_p = x_k(9);

    alfa_fus = atan2(w,u);
    beta_fus = atan2(v,V);

    %% Rotor Model
    alfa_sp = atan2(w,u);
    alfa_cp = alfa_sp + theta_s;

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

    lambda_u = (mu_z-lambda_0_u)/2;
    lambda_l = (mu_z-lambda_0_l)/2;

    %%%%%%% FLAPPING ANGLES
    a0_u = lock/(8*nu_b2) * ((theta_0+theta_d/2)*(1+mu_x^2) + 4/3*lambda_u + 2/3*mu_x*p/Omega + twist_r*(4/5+2/3*mu_x) - 4/3*mu_x*theta_s );
    a0_l = lock/(8*nu_b2) * ((theta_0-theta_d/2)*(1+mu_x^2) + 4/3*lambda_l + 2/3*mu_x*p/Omega + twist_r*(4/5+2/3*mu_x) - 4/3*mu_x*theta_s );
    a1_u = (8/3*mu_x*(theta_0+theta_d/2) + 2*mu_x*lambda_u + p/Omega - 16/lock*q/Omega + 2*twist_r*mu_x - (1+3/2*mu_x^2)*theta_s) / (1-1/2*mu_x^2); 
    a1_l = (8/3*mu_x*(theta_0-theta_d/2) + 2*mu_x*lambda_l + p/Omega - 16/lock*q/Omega + 2*twist_r*mu_x - (1+3/2*mu_x^2)*theta_s) / (1-1/2*mu_x^2); 
    b1_u = -8/lock*(nu_b2-1)/(1+1/2*mu_x^2)*a1_u + (4/3*mu_x*a0_u + q/Omega - 16/lock*p/Omega + (1+1/2*mu_x^2)*(theta_c+theta_cdiff/2))/(1+1/2*mu_x^2);
    b1_l = -8/lock*(nu_b2-1)/(1+1/2*mu_x^2)*a1_l + (4/3*mu_x*a0_l + q/Omega - 16/lock*p/Omega + (1+1/2*mu_x^2)*(theta_c-theta_cdiff/2))/(1+1/2*mu_x^2);

    %%%%%%% H-FORCES
    CHdp_u = sigma*0.015*mu_x/4 + sigma*c_l_a/4*((a1_u*mu_x^2/2 + mu_x*lambda_u)*(theta_0+theta_d/2) ...
        + q/Omega*(-a0_u/3) + (a0_u^2+a1_u^2)*mu_x/2);
    CHdp_l = sigma*0.015*mu_x/4 + sigma*c_l_a/4*((a1_l*mu_x^2/2 + mu_x*lambda_l)*(theta_0-theta_d/2) ...
        + q/Omega*(-a0_l/3) + (a0_l^2+a1_l^2)*mu_x/2);

    Hdp_u = CHdp_u * rho * (Omega*R)^2 * pi*R^2;
    Hdp_l = CHdp_l * rho * (Omega*R)^2 * pi*R^2;

    %%%%%%% S-FORCES
    CSdp_u = sigma*c_l_a/4 * (-1/2*mu_x*a0_u*(theta_0+theta_d/2) + (-a0_u*mu_x/3+b1_u*mu_x^2/4-q/(Omega*4))*twist_r ...
        - 3*a0_u*mu_x*(mu_x*a1_u-lambda_u) + b1_u*(mu_x*a1_u-lambda_u)/2 + a0_u*a1_u*(mu_x^2+1)/3);
    CSdp_l = sigma*c_l_a/4 * (-1/2*mu_x*a0_l*(theta_0-theta_d/2) + (-a0_l*mu_x/3+b1_l*mu_x^2/4-q/(Omega*4))*twist_r ...
        - 3*a0_l*mu_x*(mu_x*a1_l-lambda_l) + b1_l*(mu_x*a1_l-lambda_l)/2 + a0_l*a1_l*(mu_x^2+1)/3);
    
    Sdp_u = CSdp_u * rho * (Omega*R)^2 * pi*R^2;
    Sdp_l = CSdp_l * rho * (Omega*R)^2 * pi*R^2;

    % Rotor Angles
    a1r_u = theta_s - a1_u;
    a1r_l = theta_s - a1_l;
    b1r_u = b1_u + theta_c + theta_cdiff/2;
    b1r_l = b1_l + theta_c - theta_cdiff/2;

    alfa_dp_u = alfa_cp + a1_u;
    alfa_dp_l = alfa_cp + a1_l;

    %%%%%%%% Thrust coefficients
        % Upper Rotor
        C_T_Glau_u = 2*lambda_0_u * sqrt( (mu*cos(alfa_dp_u))^2 + (mu*sin(alfa_dp_u) + lambda_0_u)^2 );       %% NOTE THE 1 IN THE FRONT OF THE EQ.!!
        C_T_BEM_u = sigma*c_l_a/2*( (1/3+mu_x^2/2)*(theta_0+theta_d/2) + (1+mu_x^2)/8*twist_r + lambda_u/2 );

        % Lower Rotor
        C_T_Glau_l = 2*lambda_0_l * sqrt( (mu*cos(alfa_dp_l))^2 + (mu*sin(alfa_dp_l) + lambda_0_l + lambda_0_u)^2 );
        C_T_BEM_l = sigma*c_l_a/2*( (1/3+mu_x^2/2)*(theta_0-theta_d/2) + (1+mu_x^2)/8*twist_r + lambda_l/2 );

    % Thrust using BEM
    T_u = C_T_BEM_u * rho * (Omega*R)^2 * pi*R^2;
    T_l = C_T_BEM_l * rho * (Omega*R)^2 * pi*R^2;

    %%%%%%% ROTOR TORQUES
    lambda_dp_u = V/(Omega*R)*sin(alfa_dp_u) - lambda_0_u;
    lambda_dp_l = V/(Omega*R)*sin(alfa_dp_l) - lambda_0_l;
    CQdp_u = sigma*(CD_r/8*(1+4.7*mu_x^2) - C_T_BEM_u*lambda_dp_u - CHdp_u*mu_x);
    CQdp_l = -sigma*(CD_r/8*(1+4.7*mu_x^2) - C_T_BEM_l*lambda_dp_l - CHdp_l*mu_x);

    Qdp_u = CQdp_u * rho * (Omega*R)^2 * R * pi*R^2;
    Qdp_l = CQdp_l * rho * (Omega*R)^2 * R * pi*R^2;

    %% Pusher Propeller Model
    alfa_sp_p = atan2(w,u);

    mu_xp = sqrt(v^2 + (w + Kp*Omega*R*(lambda_0_u+2*lambda_0_u) + q*l_p)^2)/(omega_p*R_p);
    lambda_p = -(u-q*h_p-r*d_p)/(omega_p*R_p) - lambda_0_p;

    C_T_Glau_p = 2*lambda_0_p*sqrt( (mu*sin(alfa_sp_p))^2 + ...
        (mu*cos(alfa_sp_p)+lambda_0_p)^2 );       % this is like climbing flight for a rotor
    C_T_BEM_p = sigma_p*c_l_a_p/2*( (1/3+mu_xp^2/2)*theta_p + (1+mu_xp^2)/8*twist_p + lambda_p/2 );
    Tp = C_T_BEM_p * rho * (omega_p*R_p)^2 * pi*R_p^2;

    %% HORIZONTAL STABILIZER AND ELEVATOR
    % horizontal stabilizer
    alpha_ht = -360:360;
    Cl_ht = zeros(size(alpha_ht));
    range_indices = (alpha_ht >= -20 & alpha_ht <= 20);
    Cl_ht(range_indices) = interp1([-20, 20], [-1.2, 1.2], alpha_ht(range_indices), 'linear');
    alfa_ht = atan2((w+q*l_h), u) + alfa_ht_0;
    
    Cl_ht_interp = interp1(alpha_ht, Cl_ht, rad2deg(alfa_ht), 'linear', 'extrap');
    V_ht = sqrt(u^2 + (w+q*l_h)^2);
    L_h = 0.5*rho*V_ht^2*S_ht*(Cl_ht_interp+(0.7*rad2deg(delta_e)));           % 0.859 comes from george's thesis

    %% VERTICAL STABILIZER AND RUDDER
    % vertical stabilizer
    beta_vt = -360:360;
    Cl_vt = zeros(size(beta_vt));
    range_indices = (beta_vt >= -20 & beta_vt <= 20);
    Cl_vt(range_indices) = interp1([-20, 20], [-1.2, 1.2], beta_vt(range_indices), 'linear');
    beta_act_vt = atan2((v + p*h_v - r*l_v), u) + beta_vt_0;
    
    Cl_vt_interp = interp1(beta_vt, Cl_vt, rad2deg(beta_act_vt), 'linear', 'extrap');
    V_vt = sqrt(u^2 + (v + p*h_v - r*l_v)^2);
    L_v = 0.5*rho*V_vt^2*S_vt*(Cl_vt_interp+(0.7*rad2deg(delta_r)));           % 0.859 comes from george's thesis

    %% Forces and Moments
    %%%%% X FORCE COMPONENTS
    X_MR_u = T_u*sin(a1r_u+gamma_s) - Hdp_u*cos(a1r_u+gamma_s);
    X_MR_l = T_l*sin(a1r_l+gamma_s) - Hdp_l*cos(a1r_l+gamma_s);
    X_MR = X_MR_u + X_MR_l;
    X_fus = -0.5*rho*V^2*CDS*cos(alfa_fus);
    X_prop = Tp;

    %%%%% Y FORCE COMPONENTS
    Y_MR_u = T_u*sin(b1_u) + Sdp_u*cos(b1_u);
    Y_MR_l = T_l*sin(b1_l) + Sdp_l*cos(b1_l);
    Y_vt = -L_v;

    %%%%% Z FORCE COMPONENTS
    Z_MR_u = -T_u*cos(a1r_u+gamma_s) + Hdp_u*sin(a1r_u+gamma_s);
    Z_MR_l = -T_l*cos(a1r_l+gamma_s) + Hdp_l*sin(a1r_l+gamma_s);
    Z_MR = Z_MR_u + Z_MR_l;
    Z_fus = -0.5*rho*V^2*CDS*sin(alfa_fus);
    Z_ht = -L_h*cos(alfa_ht_0);

    %%%%% L MOMENT COMPONENTS
    L_MR_u = Y_MR_u*h_u - Z_MR_u*y_MR;
    L_MR_l = Y_MR_l*h_l - Z_MR_l*y_MR;
    L_hinge_u = (Omega*R)^2*hinge_offset_ratio*R*m_bl*sin(b1r_u);
    L_hinge_l = (-Omega*R)^2*hinge_offset_ratio*R*m_bl*sin(b1r_l);
    L_vt = h_v*Y_vt;

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

    %%%%% N MOMENT COMPONENTS
    N_MR_u = Qdp_u + X_MR_u*y_MR - Y_MR_u*x_MR;
    N_MR_l = Qdp_l + X_MR_l*y_MR - Y_MR_l*x_MR;
    N_vt = -l_v*Y_vt;
    C_N_fus = (V/(Omega*R))^2 * 1/(pi*R^2*R) * K_fus * Vol_fus * beta_fus;
    N_fus = rho*pi*R^2*(Omega*R^2)*R*C_N_fus;

    %%%%% TOTAL FORCES AND MOMENTS
    X = X_MR + X_fus + X_prop;
    Y = Y_MR_u + Y_MR_l + Y_vt;
    Z = Z_MR + Z_fus + Z_ht;
    L = L_MR_u + L_MR_l + L_hinge_u + L_hinge_l + L_vt;
    M = M_MR_u + M_MR_l + M_fus + M_ht + M_hinge + M_flap;
    N = N_MR_u + N_MR_l + N_vt + N_fus;


    %% Equations of Motion
    udot = X/mass - g*sin(theta_f) - q*w + r*v;
    vdot = Y/mass + g*cos(theta_f)*sin(phi_f) - r*u + p*w ;
    wdot = Z/mass + g*cos(theta_f)*cos(phi_f) - p*v + q*u;

    pdot = ((Iyy*Izz-Izz^2-Ixz^2)*r*q + (Ixx-Iyy+Izz)*Ixz*p*q + Izz*L + Ixz*N)/(Ixx*Izz-Ixz^2);
    qdot = (M + (Izz-Ixx)*p*r - Ixz*(r^2-p^2)) / Iyy;
        qdot_MR_u = M_MR_u/Iyy;
        qdot_MR_l = M_MR_l/Iyy;
        qdot_hinge = M_hinge/Iyy;
        qdot_flap = M_flap/Iyy;
        qdot_fus = M_fus/Iyy;
        qdot_ht = M_ht/Iyy;
    rdot = ((Ixx^2-Ixx*Iyy+Ixz^2)*p*q - (Ixx-Iyy+Izz)*Ixz*q*r + Ixz*L + Ixx*N)/(Ixx*Izz - Ixz^2);

    lambda_0_u_dot = (C_T_BEM_u-C_T_Glau_u)/0.1;
    lambda_0_l_dot = (C_T_BEM_l-C_T_Glau_l)/0.1;
    lambda_p_dot = (C_T_BEM_p-C_T_Glau_p)/0.1;

    f_xki_ti = [udot; vdot; wdot; pdot; qdot; rdot; lambda_0_u_dot; lambda_0_l_dot; lambda_p_dot];
    
    M_components = [M_MR_u, M_MR_l, M_hinge, M_flap, M_fus, M_ht];
    Mw_components = [qdot_MR_u, qdot_MR_l, qdot_hinge, qdot_flap, qdot_fus, qdot_ht];
    angles = rad2deg([a0_u, a0_l, a1_u, a1_l, a1r_u, a1r_l, alfa_sp, alfa_cp, alfa_dp_u, alfa_dp_l]);
    latAngles = rad2deg([b1_u, b1_l, b1r_u, b1r_l]);
    inflow = [lambda_0_u, lambda_0_l];

    
end




