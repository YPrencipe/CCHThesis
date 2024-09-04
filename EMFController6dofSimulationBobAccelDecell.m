%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename:     EMFController6dofSimulationBobAccellDecell
% Project:      MSc Thesis
% Supervisor:   M.D. Pavel
% Author:       Ynias Prencipe 
% Student Nr.:  4777158
% 
% Description:  controller with several manoeuvres including bob-up/down
%               with accell/decell
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc; %close all;

% Load Newest Trim States and Helicopter Parameters
load('trim_saved_6dof.mat'); load('TrimLinSave_6dof.mat'); 
load('V_transit_save.mat'); 
load('ref_3211'); 
load('theta_cmd_3211'); load('phi_cmd_3211'); load('r_cmd_3211');
load('V_x_cmd'); load('V_z_cmd');
load('t'); 
load('X_theta_x.mat'); load('Z_theta_z.mat');
coaxial_heli_parameters;

%% Simulation Parameters

sys6DoF = ss(TrimLinSave_6dof{9, 1}, TrimLinSave_6dof{9, 2}, ...
    TrimLinSave_6dof{9, 3}, TrimLinSave_6dof{9, 4});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ManoeuvreSelection = 2; % 1. 3-2-1-1 
                        % 2. Bob up/down w accel decell
                        % 3. Step input for attitude quickness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ManoeuvreSelection == 1
    subSelection3211 = input('Choose attitude \n1 = theta, 2 = phi, 3 = psi : ');
end


if ManoeuvreSelection == 1
    tEnd = 10;
    simN = 400;
    u(1)=40;
    v(1)=0;
    w(1)=0;
    PitchDropBack = 0;
end
if ManoeuvreSelection == 2
    tEnd = 70;
    simN = 1600;
    u(1)=40;
    v(1)=0;
    w(1)=0;
    theta_cmd(1) = 0;
    phi_cmd(1) = 0;
    PitchDropBack = 0;
end
if ManoeuvreSelection == 3
    tEnd = 5;
    simN = 400;
    u(1)=40;
    v(1)=0;
    w(1)=0;
    theta_cmd(1) = 0;
    phi_cmd(1) = 0;
    attitude = input('Choose attitude, 1 pitch, 2 roll, 3 yaw: ');
    distVal = input('Choose step value in degreees: ');
    desired_Attitude = deg2rad(distVal);
    PitchDropBack = 0;
end

% TRIM SPEED
[~, idx_trimspeed] = min(abs(V_vals - u(1)));


% U = zeros(6,simN);
U(1,1) = trim_saved_6dof(1,idx_trimspeed);  % Coll trim
U(2,1) = trim_saved_6dof(2,idx_trimspeed);  % DiffColl trim
U(3,1) = trim_saved_6dof(3,idx_trimspeed);  % LonCyc trim
U(4,1) = trim_saved_6dof(4,idx_trimspeed);  % LatCyc trim
U(5,1) = trim_saved_6dof(5,idx_trimspeed);  % roll trim
U(6,1) = trim_saved_6dof(6,idx_trimspeed);  % PropColl trim
U(7,1) = trim_saved_6dof(7,idx_trimspeed);  % quasidynamic inflow state trim upper rotor
U(8,1) = trim_saved_6dof(8,idx_trimspeed);  % quasidynamic inflow state trim lower rotor
U(9,1) = trim_saved_6dof(9,idx_trimspeed);  % quasidynamic inflow state trim pusher propeller
t(1)=0;

p(1)=0;
q(1)=0;
r(1)=0;
phi(1)=0;
theta(1)=0;
psi(1)=0;

% Time Paramaters
dt = (tEnd-t(1))/simN;

% CA ON or OFF
CA = 1;
V_l_transit = V_transit_save(1)/1.94384;
V_u_transit = V_transit_save(2)/1.94384-20;

% Control Parameters Pitch
if CA == true
    generalFactor = 1;
    % if ManoeuvreSelection == 1
    %     if u(1) == 0
    %         Kp_theta = 80;
    %         Ki_theta = 0.1;
    %         Kp_q = -10;
    %         r_theta = 1;
    %         Kd_theta = 0;
    % 
    %         Kp_phi = 120;
    %         Ki_phi = 0.8;
    %         Kp_p = -2.8;
    %         Kd_phi = 0.5;
    %         r_phi = 1 * generalFactor;
    % 
    %         Kp_r = 25;
    %         Ki_r = 0;
    %         r_r = 1;
    %     end
    %     if u(1) == 40
    %         Kp_theta = 250;
    %         Ki_theta = 10;
    %         Kp_q = -25;
    %         r_theta = 1;
    %         Kd_theta = 0;
    % 
    %         Kp_phi = 120;
    %         Ki_phi = 0.8;
    %         Kp_p = -2.8;
    %         Kd_phi = 0.5;
    %         r_phi = 1 * generalFactor;
    % 
    %         Kp_r = 25;
    %         Ki_r = 25;
    %         r_r = 1;
    %     end
    % end

    if ManoeuvreSelection == 1 || 2 || 3
        Kp_theta = 80;
        Ki_theta = 15;
        Kp_q = -20;
        r_theta = 1;
        Kd_theta = 0;

        Kp_phi = 120;
        Ki_phi = 0.8;
        Kp_p = -2.8;
        Kd_phi = 0.5;
        r_phi = 1 * generalFactor;

        Kp_r = 1;
        Ki_r = 1;
        r_r = 1;

        if ManoeuvreSelection == 1
            Kp_psi = 50;
            Ki_psi = 0.1;
            Kp_r = -8;
        elseif ManoeuvreSelection == 3
            Kp_psi = 50;
            Ki_psi = 0.1;
            Kp_r = -8;
        end
    end

    % factor = 1;
    factor = generalFactor;

    Kp_x = -0.9 * factor;       % for manoeuvre 3
    % Kp_x = -0.1 * factor;    %for manoeuvre 2
    Ki_x = -0 * factor;

    Kp_z = 5 * factor;
    % Kp_z = 8 * factor;
    Ki_z = 0.1 * factor;
    Kd_z = 0.0 * factor;

    Kp_y = -deg2rad(1.5) * factor;
    % Kp_y = -deg2rad(1.5) * 0;
    Ki_y = 0 * factor;
end
if CA == false
    if ManoeuvreSelection == 1
        Kp_theta = -0.9;
        Ki_theta = -0.9;
        Kp_q = 0.6;
        r_theta = 1.4;
    end
    if ManoeuvreSelection == 2 || 3
        Kp_theta = -0.3;
        Ki_theta = -0.1;
        Kp_q = 0.07;
        r_theta = 1;
    end

    Kp_x = -deg2rad(4);
    Ki_x = -deg2rad(0.1);

    Kp_z = 0;
    Ki_z = 0;
end

% Initialize list to store Q values
Q_values = [];
r_values = [];

% Iterate over different values of r
% for r = 0.2:0.1:4

% EMF Control Scaling Parameter
r_x = 1;
r_y = 1;
r_z = 1;
integral_phi = 0;
integral_theta = 0;
integral_r = 0;
integral_x = 0;
integral_y = 0;
integral_z = 0;

% Reset t
t = 0;

%% Simulation
for i = 1:simN

    
    % Speed definitions
    V(i) = sqrt(u(i)^2 + v(i)^2 + w(i)^2);
    vel(:,i) = [u(i); v(i); w(i); 0; 0; 0];  

    if ManoeuvreSelection == 2 || 3
        if PitchDropBack == true
            if t(i)>4
                Kp_theta = 0;
                Ki_theta = 0;
                Kp_q = 0;
                r_theta = 0;
                Kd_theta = 0;
        
                Kp_phi = 0;
                Ki_phi = 0;
                Kp_p = 0;
                Kd_phi = 0;
                r_phi = 0;
        
                Kp_r = 0;
                Ki_r = 0;
                r_r = 0;
            
                Kp_x = 0;
                Ki_x = 0;
            
                Kp_z = 0;
                Ki_z = 0;
                Kd_z = 0;
            
                Kp_y = 0;
                Ki_y = 0;
            end  
        end
    end

    % Keep certain inputs constant for t
    % U(2,i) = U(2,1);    % DiffColl
    % U(5,i) = U(5,1);    % Roll
    % delta_r(i) = 0;
    theta_cdiff(i) = 0;

    % Select Derivatives for Correct Speed
    [~, idx_speed] = min(abs(V_vals - V(i)));

    L_theta_c = TrimLinSave_6dof{idx_speed,2}(4,4);
    M_theta_s = TrimLinSave_6dof{idx_speed,2}(5,3);
    M_delta_e = TrimLinSave_6dof{idx_speed,2}(5,7);
    N_theta_d = TrimLinSave_6dof{idx_speed,2}(6,2);
    N_delta_r = TrimLinSave_6dof{idx_speed,2}(6,8);
    % X_theta_x = X_theta_x(idx_speed);
    X_theta_p = TrimLinSave_6dof{idx_speed,2}(1,6);
    Z_theta_0 = TrimLinSave_6dof{idx_speed,2}(3,1);
    % Z_theta_z = Z_theta_z(idx_speed);

    L_v = TrimLinSave_6dof{idx_speed,1}(4,2);

    if ManoeuvreSelection == 3
        if t(i) < 1
            ref_step(i) = 0;
        else t(i) >= 1
            ref_step(i) = desired_Attitude;
        end
    end

    % Inner Loop ACAH Controllers

    if ManoeuvreSelection == 1
        if subSelection3211 == 1
            error_theta(i) = theta_cmd_3211(i) - theta(i);
            phi_cmd_3211(i) = 0;
            error_phi(i) = phi_cmd_3211(i) - phi(i);
            r_cmd_3211(i) = 0;
            error_r(i) = r_cmd_3211(i) - r(i);
        elseif subSelection3211 == 2
            error_theta(i) = 0 - theta(i);
            error_phi(i) = theta_cmd_3211(i) - phi(i);
            r_cmd_3211(i) = 0;
            error_r(i) = r_cmd_3211(i) - r(i);
        else
            error_theta(i) = 0 - theta(i);
            phi_cmd_3211(i) = 0;
            error_phi(i) = phi_cmd_3211(i) - phi(i);
            error_r(i) = theta_cmd_3211(i) - psi(i);
        end
    end
    if ManoeuvreSelection == 2
        % phi_cmd(i) = 0;
        error_phi(i) = phi_cmd(i) - phi(i);
        error_theta(i) = theta_cmd(i) - theta(i);
        psi_cmd(i) = 0;
        error_psi(i) = psi_cmd(i) - psi(i);
        r_cmd(i) = 0;
        error_r(i) = r_cmd(i) - r(i);
    end
    if ManoeuvreSelection == 3
        if attitude == 1
            phi_cmd(i) = 0;
            error_phi(i) = phi_cmd(i) - phi(i);
            theta_cmd(i) = ref_step(i);
            error_theta(i) = theta_cmd(i) - theta(i);
            r_cmd(i) = 0;
            error_r(i) = r_cmd(i) - r(i);
        elseif attitude == 2
            phi_cmd(i) = ref_step(i);
            error_phi(i) = phi_cmd(i) - phi(i);
            theta_cmd(i) = 0;
            error_theta(i) = theta_cmd(i) - theta(i);
            r_cmd(i) = 0;
            error_r(i) = r_cmd(i) - r(i);
        else
            phi_cmd(i) = 0;
            error_phi(i) = phi_cmd(i) - phi(i);
            theta_cmd(i) = 0;
            error_theta(i) = theta_cmd(i) - theta(i);
            psi_cmd(i) = ref_step(i);
            error_r(i) = psi_cmd(i) - psi(i);
        end
    end

    integral_phi = integral_phi + error_phi(i) * dt;
    integral_theta = integral_theta + error_theta(i) * dt;


    if i>1
        gamma_phi(i) = r_phi * (Kp_phi*error_phi(i) + Ki_phi * integral_phi + ... 
            Kp_p*p(i) + Kd_phi * (error_phi(i) - error_phi(i-1)) / dt);
        gamma_theta(i) = r_theta * (Kp_theta*error_theta(i) + Ki_theta ... 
            * integral_theta + Kp_q*q(i) + Kd_theta*(error_theta(i)-error_theta(i-1))/dt);
    else
        gamma_phi(i) = r_phi * (Kp_phi*error_phi(i) + Ki_phi * integral_phi + ... 
            Kp_p*p(i));
        gamma_theta(i) = r_theta * (Kp_theta*error_theta(i) + Ki_theta ... 
            * integral_theta + Kp_q*q(i) );

    end

   
    
    if ManoeuvreSelection  == 2
        integral_r = integral_r + error_r(i) * dt;
        gamma_r(i) = r_r * (Kp_r*error_r(i) + Ki_r * integral_r );
    elseif ManoeuvreSelection == 1
        integral_r = integral_r + error_r(i) * dt;
        gamma_r(i) = r_r * (Kp_psi*error_r(i) + Ki_psi * integral_r + Kp_r*r(i));
    else
        integral_r = integral_r + error_r(i) * dt;
        gamma_r(i) = r_r * (Kp_psi*error_r(i) + Ki_psi * integral_r + Kp_r*r(i));
    end

    % Outer Loop Controllers
    error_x(i) = V_x_cmd(i) - u(i);
    integral_x = integral_x + error_x(i) * dt;
    gamma_x(i) = r_x * (Kp_x*error_x(i) + Ki_x*integral_x);

    error_z(i) = V_z_cmd(i) - w(i);
    integral_z = integral_z + error_z(i) * dt;

    if i > 1
        gamma_z(i) = r_z * (Kp_z*error_z(i) + Ki_z*integral_z + Kd_z * (error_z(i) - error_z(i-1))/dt );
    else
        gamma_z(i) = r_z * (Kp_z*error_z(i) + Ki_z*integral_z);
    end

    V_y_cmd(i) = 0;
    error_y = V_y_cmd(i) - v(i);
    integral_y = integral_y + error_y * dt;
    gamma_y(i) = r_y * (Kp_y*error_y + Ki_y*integral_y);

    phi_cmd(i+1) = gamma_y(i);
    % phi_cmd(i+1) = deg2rad(1);

    % If not using Control Allocation
    if CA == false
        % Pitch
        U(3,i) = gamma_theta(i);
        delta_e(i) = 0;

        % LonVel
        theta_x_cmd(i) = gamma_x(i);
        theta_p(i) = U(6,1);

        % Height
        theta_z_cmd(i) = 0;
        U(1,i) = U(1,1);

    end

    % Define Reference 3-2-1-1 Signal
    if ManoeuvreSelection == 1
        if t(i) < 1
            ref_3211(i) = 0;
        elseif t(i) >= 1 && t(i) < 4 
            ref_3211(i) = deg2rad(3);
        elseif t(i) >= 4 && t(i) < 6 
            ref_3211(i) = deg2rad(-3);
        elseif t(i) >= 6 && t(i) < 7 
            ref_3211(i) = deg2rad(3);
        elseif t(i) >= 7 && t(i) < 8 
            ref_3211(i) = deg2rad(-3);
        else t(i) >= 8;
            ref_3211(i) = 0;
        end
    end

    

    % Define Bob-Up Bob-Down with Accel/Decell Manoeuvre
    if t(i) < 5
        V_x_ref(i) = u(1);
    elseif t(i) >= 5 && t(i) < 25
        V_x_ref(i) = 30;
    elseif t(i) >= 25 && t(i) < 40
        V_x_ref(i) = 40;
    elseif t(i) >= 40
        V_x_ref(i) = 30;
    end

    if t(i) < 10
        V_z_ref(i) = w(1);
    elseif t(i) >= 10 && t(i) < 20
        V_z_ref(i) = -2;
    elseif t(i) >= 20 && t(i) < 45
        V_z_ref(i) = 0;
    elseif t(i) >= 45 && t(i) < 55
        V_z_ref(i) = 2;
    elseif t(i) >= 55
        V_z_ref(i) = 0;
    end

    % CA Inner Loop
    if V(i)<V_l_transit
        rho_theta_c(i) = 1;
        rho_delta_e(i) = epsilon;
    elseif V(i)>V_l_transit && V(i)<V_u_transit
        rho_theta_c(i) = 1 - 1/(V_u_transit-V_l_transit)*(V(i)-V_l_transit);
        rho_delta_e(i) = 1/(V_u_transit-V_l_transit)*(V(i)-V_l_transit);
    else
        rho_theta_c(i) = epsilon;
        rho_delta_e(i) = 1;
    end

    % CA Outer Loop
    if V(i)<V_l_transit
        rho_theta_x_cmd(i) = 1;
        rho_theta_p(i) = epsilon;

        rho_theta_0(i) = 1;
        rho_theta_z_cmd(i) = epsilon;

        rho_theta_d(i) = 1;
        rho_delta_r(i) = epsilon;
    elseif V(i)>V_l_transit && V(i)<V_u_transit
        rho_theta_x_cmd(i) = 1 - 1/(V_u_transit-V_l_transit)*(V(i)-V_l_transit);
        rho_theta_p(i) = 1/(V_u_transit-V_l_transit)*(V(i)-V_l_transit);

        rho_theta_0(i) = 1 - 1/(V_u_transit-V_l_transit)*(V(i)-V_l_transit);
        rho_theta_z_cmd(i) = 1/(V_u_transit-V_l_transit)*(V(i)-V_l_transit);

        rho_theta_d(i) = 1 - 1/(V_u_transit-V_l_transit)*(V(i)-V_l_transit);
        rho_delta_r(i) = 1/(V_u_transit-V_l_transit)*(V(i)-V_l_transit);
    else
        rho_theta_x_cmd(i) = epsilon;
        rho_theta_p(i) = 1;

        rho_theta_0(i) = epsilon;
        rho_theta_z_cmd(i) = 1;

        rho_theta_d(i) = epsilon;
        rho_delta_r(i) = 1;
    end


    if CA == true

        [~, idx_trimspeed] = min(abs(V_vals - V(i)));

        % LonVel CA
        B_x = [9.8066, -X_theta_p];
        W_x = [1/rho_theta_x_cmd(i), 0; 0, 1/rho_theta_p(i)];
        U_x(:,i) = inv(W_x)*transpose(B_x)* ...
            inv(B_x*inv(W_x)*transpose(B_x))*gamma_x(i);
        theta_x_cmd(i) = U_x(1,i);
        U(6,i) = U_x(2,i) + trim_saved_6dof(6,idx_trimspeed);

        % Height CA
        B_z = [Z_theta_0, 9.8066];
        W_z = [1/rho_theta_0(i), 0; 0, 1/rho_theta_z_cmd(i)];
        U_z(:,i) = inv(W_z)*transpose(B_z)* ...
            inv(B_z*inv(W_z)*transpose(B_z))*gamma_z(i);
        U(1,i) = U_z(1,i) + trim_saved_6dof(1,idx_trimspeed);
        theta_z_cmd(i) = U_z(2,i);

        % Pitch CA
        B_theta = [M_theta_s, M_delta_e];
        W_theta = [1/rho_theta_c(i), 0; 0, 1/rho_delta_e(i)];
        U_theta(:,i) = inv(W_theta)*transpose(B_theta)* ...
            inv(B_theta*inv(W_theta)*transpose(B_theta))*gamma_theta(i);
        U(3,i) = U_theta(1,i) + trim_saved_6dof(3,idx_trimspeed);
        delta_e(i) = U_theta(2,i); 

        % Roll CA
        U(4,i) = gamma_phi(i)/L_theta_c + trim_saved_6dof(4,idx_trimspeed);

        % Yaw CA
        B_psi = [N_theta_d, N_delta_r];
        W_psi = [1/rho_theta_d(i), 0; 0, 1/rho_delta_r(i)];
        U_psi(:,i) = inv(W_psi)*transpose(B_psi)* ...
            inv(B_psi*inv(W_psi)*transpose(B_psi))*gamma_r(i);
        U(2,i) = U_psi(1,i) + trim_saved_6dof(2,idx_trimspeed);
        delta_r(i) = U_psi(2,i); 

    end

    

    % Actuator Rate Constraints
    if i>1
        % LonCyc
        cyclic_rate(i) = (U(3,i)-U(3,i-1))/dt;
        if abs(cyclic_rate(i)) > 28.8
            U(3,i) = sign(cyclic_rate(i))*28.8*dt + U(3,i-1);
        end
        cyclic_rate(i) = (U(3,i)-U(3,i-1))/dt;

        % LatCyc
        cyclic_rate(i) = (U(4,i)-U(4,i-1))/dt;
        if abs(cyclic_rate(i)) > 28.8
            U(4,i) = sign(cyclic_rate(i))*28.8*dt + U(4,i-1);
        end
        cyclic_rate(i) = (U(4,i)-U(4,i-1))/dt;
        
        % Elev
        elev_rate(i) = (delta_e(i)-delta_e(i-1))/dt;
        if abs(elev_rate(i)) > 28.8
            delta_e(i) = sign(elev_rate(i))*28.8*dt + delta_e(i-1);
        end
        elev_rate(i) = (delta_e(i)-delta_e(i-1))/dt;

        % Rudd
        rudd_rate(i) = (delta_r(i)-delta_r(i-1))/dt;
        if abs(rudd_rate(i)) > 28.8
            delta_r(i) = sign(rudd_rate(i))*28.8*dt + delta_r(i-1);
        end
        rudd_rate(i) = (delta_r(i)-delta_r(i-1))/dt;

        % Coll
        coll_rate(i) = (U(1,i)-U(1,i-1))/dt;
        if abs(coll_rate(i)) > 16
            U(1,i) = sign(coll_rate(i))*16*dt + U(1,i-1);
        end
        coll_rate(i) = (U(1,i)-U(1,i-1))/dt;

        % PropColl
        propcoll_rate(i) = (U(6,i)-U(6,i-1))/dt;
        if abs(propcoll_rate(i)) > 16
            U(6,i) = sign(propcoll_rate(i))*16*dt + U(6,i-1);
        end
        propcoll_rate(i) = (U(6,i)-U(6,i-1))/dt;

    end

    % Actuator Constraints
        % LonCyc
        if abs(U(3,i)) > deg2rad(12)
            U(3,i) = sign(U(3,i))*deg2rad(12);
        end
        % Elev
        if abs(delta_e(i)) > deg2rad(20)
            delta_e(i) = sign(delta_e(i))*deg2rad(20);
        end
        % Rudd
        if abs(delta_r(i)) > deg2rad(20)
            delta_r(i) = sign(delta_r(i))*deg2rad(20);
        end
        % Coll
        if U(1,i) > deg2rad(20)
            U(1,i) = deg2rad(20);
        end
        if U(1,i) < deg2rad(0)
            U(1,i) = deg2rad(0);
        end
        % PropColl
        if U(6,i) > deg2rad(100)
            U(6,i) = deg2rad(100);
        end
        if U(6,i) < deg2rad(-10)
            U(6,i) = deg2rad(-10);
        end
        if abs(U(4,i)) > deg2rad(7)
            U(4,i) = sign(U(4,i))*deg2rad(7);
        end

    % Calc Total Required Pitch
    theta_cmd(i+1) = theta_x_cmd(i) + theta_z_cmd(i);

    

    % State Update
    [xdot(:,i)] = f_xk6(vel(:,i), U(:,i), p(i), q(i), r(i), theta(i), delta_e(i), delta_r(i), theta_cdiff(i));

    % U(1,i) = trim_saved_6dof(1,idx_trimspeed);

    % Euler Integration
    u(i+1) = u(i)+dt*xdot(1,i);
    v(i+1) = v(i)+dt*xdot(2,i);
    w(i+1) = w(i)+dt*xdot(3,i);
    p(i+1) = p(i)+dt*xdot(4,i);
    % p(i+1) = 0;
    q(i+1) = q(i)+dt*xdot(5,i);
    % q(i+1) = 0;
    r(i+1) = r(i)+dt*xdot(6,i);
    % r(i+1) = 0;
    phi(i+1) = phi(i)+dt*p(i);
    theta(i+1) = theta(i)+dt*q(i);
    psi(i+1) = psi(i)+dt*r(i);
    U(7,i+1) = U(7,i)+dt*xdot(7,i);
    U(8,i+1) = U(8,i)+dt*xdot(8,i);
    U(9,i+1) = U(9,i)+dt*xdot(9,i);
    t(i+1) = t(i)+dt;

    

end

%%
for i=1:simN
    t(i+1) = t(i)+dt;
    if t(i) < 5
        V_x_ref(i) = u(1);
    elseif t(i) >= 5 && t(i) < 25
        V_x_ref(i) = 30;
    elseif t(i) >= 25 && t(i) < 40
        V_x_ref(i) = 40;
    elseif t(i) >= 40
        V_x_ref(i) = 30;
    end

    if t(i) < 10
        V_z_ref(i) = w(1);
    elseif t(i) >= 10 && t(i) < 20
        V_z_ref(i) = -2;
    elseif t(i) >= 20 && t(i) < 45
        V_z_ref(i) = 0;
    elseif t(i) >= 45 && t(i) < 55
        V_z_ref(i) = 2;
    elseif t(i) >= 55
        V_z_ref(i) = 0;
    end
end


%%
% EMF Controller Performance (Q) Calculation
if ManoeuvreSelection == 1
    err_theta = sum(abs(theta(1:end-1) - theta_cmd_3211'));
    err_phi = sum(abs(phi(1:end-1) - theta_cmd_3211'));
    err_r = sum(abs(psi(1:end-1) - theta_cmd_3211'));
else
    err_theta = sum(abs(theta - theta_cmd));
    err_phi = sum(abs(phi - phi_cmd));
    err_r = sum(abs(r - phi_cmd));
end

err_Vx = sum(abs(V_x_ref - u(1:end-1)));
err_Vz = sum(abs(V_z_ref - w(1:end-1)));


V = 1; %deg/in
err_rel_theta = 1/(V*tEnd)*err_theta;
Q_theta = 1-err_rel_theta
err_rel_phi = 1/(V*tEnd)*err_phi;
Q_phi = 1-err_rel_phi
err_rel_r = 1/(V*tEnd)*err_r;
Q_r = 1-err_rel_r



%%
% EMF Pre-defining
save('ref_3211', 'ref_3211');

if ManoeuvreSelection == 1
    M_theta_tf = tf([4^2], [1, 2*0.707*4, 4^2]);
    theta_cmd_3211 = lsim(M_theta_tf, ref_3211, t(1:end-1));
    save('theta_cmd_3211', 'theta_cmd_3211');
    
    M_phi_tf = tf([4^2], [1, 2*0.707*4, 4^2]);
    phi_cmd_3211 = lsim(M_phi_tf, ref_3211, t(1:end-1));
    save('phi_cmd_3211', 'phi_cmd_3211');

    M_r_tf = tf([1], [1/4, 1]);
    r_cmd_3211 = lsim(M_r_tf, ref_3211, t(1:end-1));
    save('r_cmd_3211', 'r_cmd_3211');
end

if ManoeuvreSelection == 2
    M_V_z_tf = tf([3^2], [1, 2*0.7*3, 3^2]);
    V_z_cmd = lsim(M_V_z_tf, V_z_ref, t(1:end-1), V_z_ref(1));
    save('V_z_cmd', 'V_z_cmd');

    A = tf2ss([1^2], [1, 2*0.7*1, 1^2]);
    B = [0; -1];
    C = [1 0];
    D = 0;
    sys = ss(A,B,C,D);
    V_x_cmd = lsim(sys, V_x_ref, t(1:end-1), [40,-40*0.7/0.5]);
    V_x_cmd = V_x_cmd';
    save('V_x_cmd', 'V_x_cmd');
end

save('t', 't');




%% PLOTTING
% PLOT ACAH Response
% figure(71)
% subplot(3,1,1)
% plot(t(1:end-1), rad2deg(ref_3211), t(1:end-1), rad2deg(theta_cmd), ...
%      t(1:end-1), rad2deg(theta(1:end-1)));
% grid on; legend('3-2-1-1 ref', '\theta_{cmd} model', '\theta_f'); 
% ylabel('\theta [deg]'); xlabel('Time [s]');
% subplot(3,1,2)
% plot(t, rad2deg(U(2,:)), t(1:end-1), rad2deg(delta_e(:))); 
% grid on;
% legend('LonCyc', 'Elev');
% ylabel('Actuator Deflection [deg]')
% xlabel('Time [s]')
% subplot(3,1,3)
% plot(t,rad2deg(q)); grid on; legend('q [deg/s]');

% Model-Following Performance Criterion
% figure(72)
% scatter(r_values, Q_values*100, 100, 'Marker', 's', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'blue', 'LineWidth', 2);
% hold on;
% scatter(r_values, Q_values*100, 100, 'Marker', '.', 'MarkerEdgeColor', 'blue');
% plot(r_values, Q_values*100, '-b'); hold off;
% grid on;
% xlabel('r');
% ylabel('Q [%]');  
% title('Optimization of Control System Activity');

% Plot Bob-Up/Down accell decell response
% if ManoeuvreSelection == 2
%     figure(81)
%     subplot(4,2,1)
%     plot(t(1:end-1), V_x_ref, t(1:end-1), V_x_cmd, t, u); grid on;
%     ylabel('V_x [m/s]'); %ylim([-10;10])
%     legend('Reference Signal', 'Model Signal')
%     subplot(4,2,3)
%     plot(t(1:end-1), V_y_cmd, t, v); grid on;
%     ylabel('V_y [m/s]')
%     subplot(4,2,5)
%     plot(t(1:end-1), V_z_ref, t(1:end-1), V_z_cmd, t, w); grid on;
%     ylabel('V_z [m/s]')
%     xlabel('Time [s]')
% 
%     subplot(4,2,2)
%     plot(t, rad2deg(phi), t, rad2deg(phi_cmd)); grid on;
%     ylabel('\phi [deg]')
%     subplot(4,2,4)
%     if ManoeuvreSelection == 1
%         plot(t(1:end-1), rad2deg(theta_cmd_3211), t, rad2deg(theta)); grid on;
%     else
%         plot(t, rad2deg(theta_cmd), t, rad2deg(theta)); grid on;
%     end
%     ylabel('\theta [deg]')
%     legend('\theta [deg]', '\theta_{cmd} [deg]')
%     subplot(4,2,6)
%     plot(t(1:end-1), psi_cmd, t, rad2deg(psi)); grid on;
%     ylabel('\psi [deg]')
%     subplot(4,2,8)
%     plot(t(1:end-1), r_cmd, t, rad2deg(r)); grid on;
%     ylabel('r [deg/s]')
%     xlabel('Time [s]')
%     subplot(4,2,7)
%     plot(t, rad2deg(p))
%     % 
%     figure(82)
%     subplot(4,2,1)
%     plot(t(1:end-1), rad2deg(U(1,1:end-1))); legend('Coll'); grid on;
%     subplot(4,2,3)
%     plot(t(1:end-1), rad2deg(U(2,1:end-1))); legend('DiffColl'); grid on;
% 
%     subplot(4,2,5)
%     plot(t(1:end-1), rad2deg(U(3,1:end-1))); legend('LonCyc'); grid on;
%     subplot(4,2,7)
%     plot(t(1:end-1), rad2deg(U(4,1:end-1))); legend('LatCyc'); grid on;
% 
%     subplot(4,2,2)
%     plot(t(1:end-1), rad2deg(U(6,1:end-1))); legend('PropColl');grid on;
% 
%     subplot(4,2,4)
%     plot(t(1:end-2), rad2deg(delta_e(1:end-1))); legend('Elev'); grid on;
%     subplot(4,2,6)
%     plot(t(1:end-2), rad2deg(delta_r(1:end-1))); legend('Rudder'); grid on;
% end

%%
% if ManoeuvreSelection == 1
%     figure(86)
%     subplot(3,1,1)
%     plot(t(1:end-1), rad2deg(ref_3211), t(1:end-1), rad2deg(theta_cmd_3211), ...
%         t, rad2deg(theta)); grid on; ylabel('\theta')
%     subplot(3,1,2)
%     plot(t(1:end-1), rad2deg(phi_cmd_3211), t, rad2deg(phi)); grid on; ylabel('\phi')
%     subplot(3,1,3)
%     plot(t(1:end-1), rad2deg(r_cmd_3211), t, rad2deg(r)); grid on; ylabel('r')
% 
% 
%     figure(87)
%     subplot(3,1,1)
%     plot(t,rad2deg(U(3,:)), t(1:end-1), rad2deg(delta_e)); legend('LonCyc', 'Elev'); grid on;
%     subplot(3,1,2)
%     plot(t,rad2deg(U(4,:))); legend('LatCyc'); grid on;
%     subplot(3,1,3)
%     plot(t,rad2deg(U(2,:)), t(1:end-1),rad2deg(delta_r)); legend('Differential Collective', 'Rudder'); grid on;
% 
% end

%% HQA
if ManoeuvreSelection == 3
    load('qpitchhq1.mat');load('qpitchhq2.mat');

    if attitude == 1

        figure(20)
    
        subplot(3,1,1)
        plot(t(1:end-1), rad2deg(U(3,1:end-1)), 'LineWidth', 2); hold on;
        plot(t(1:end-1), rad2deg(delta_e), 'LineWidth', 2); hold off;
        grid on; 
        legend('LonCyc [deg]', 'Elev [deg]', 'Interpreter', ...
            'latex', 'FontSize', 10, 'Location', 'SouthEast')
        ylabel('Input [deg]', 'Interpreter', 'latex', 'FontSize', 14)
        
        subplot(3,1,2)
        plot(t, rad2deg(q), 'LineWidth', 2); 
        grid on; 
        max_q_value = max(q);
        q_pk_deg = rad2deg(max_q_value)
        disp(['Maximum Value of q: ', num2str(max_q_value), 'Interpreter', 'latex']);
        hold on;
        yline(rad2deg(max_q_value), '--r');
        hold off;
        ylabel('q [deg/s]', 'Interpreter', 'latex', 'FontSize', 14)
        ylim([rad2deg(min(q))*1.2, rad2deg(max_q_value)*1.4]);
        xlim([0,t(end)]);
        
        subplot(3,1,3)
        plot(t, rad2deg(theta), 'LineWidth', 2); 
        grid on;
        hold on;
        plot(t(1:end-1), rad2deg(ref_step), 'm:', 'LineWidth', 1.5); 
        ylim([rad2deg(min(theta))*1.2, rad2deg(max(theta)*1.4)])
        hold off;
        ylabel('Attitude [deg]', 'Interpreter', 'latex', 'FontSize', 14)
        xlim([0,t(end)]);
        
        % Calculate rise time for theta
        percent_10 = 0.1; % 10%
        percent_90 = 0.9; % 90%
        steady_state_value_theta = mean(rad2deg(theta(end-100:end))); % Assuming the steady state is the mean of last 100 samples
        ten_percent_value_theta = percent_10 * steady_state_value_theta;
        ninety_percent_value_theta = percent_90 * steady_state_value_theta;
        idx_10_theta = find(rad2deg(theta) > ten_percent_value_theta, 1, 'first');
        idx_90_theta = find(rad2deg(theta) > ninety_percent_value_theta, 1, 'first');
        rise_time_theta = t(idx_90_theta) - t(idx_10_theta);
        
        % Find maximum value of theta
        max_value_theta = max(theta);
        theta_pk_deg = rad2deg(max_value_theta)
        
        % Find minimum value after the rise time for theta
        min_value_after_rise_theta = min(theta(idx_90_theta:end));
        rad2deg(min_value_after_rise_theta)
        
        % Plot horizontal red dashed lines for theta
        hold on;
        yline(rad2deg(max_value_theta), '--r');
        yline(rad2deg(min_value_after_rise_theta), '--r');
        hold off;
        legend('$\theta$ [deg]', 'Interpreter', 'latex', 'FontSize', 10, ...
            'Location', 'SouthEast')
        
        % Display results for theta
        disp(['Rise Time for \theta: ', num2str(rise_time_theta)]);
        disp(['Maximum Value for \theta: ', num2str(max_value_theta)]);
        disp(['Minimum Value after Rise Time for \theta: ', num2str(min_value_after_rise_theta)]);
        
        % %% PLOT PITCH ATTITUDE QUICKNESS
        Qpitch_5deg = max_q_value / (max_value_theta - theta(1))
        figure(21)
        plot(rad2deg(min_value_after_rise_theta), Qpitch_5deg, 'p', 'MarkerSize', 10)
        grid on; hold on;
        plot(qpitchhq1(:,1), qpitchhq1(:,2), 'k--', 'LineWidth', 2);
        plot(qpitchhq2(:,1), qpitchhq2(:,2), 'k--', 'LineWidth', 2);
        % Add text labels for the different levels
        text(10, 1.8, 'Level 1', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
        text(10, 0.8, 'Level 2', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
        text(10, 0.2, 'Level 3', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
        
            otherQvals = [  2.9358, 2.9714, 2.9358, 2.9057, 2.8670, 2.8309];
            QvalsX = [4.6056 9.2817 13.8512 18.5795 23.1041 27.9526];
            rad2deg(min_value_after_rise_theta)
            p = plot(QvalsX(1:6), otherQvals(1:6), 'o', 'MarkerSize', 10);
            p.MarkerFaceColor = [1 0.5 0];
    
        hold off;
        xlim([0, 30]); ylim([0, 3.5]);
        xlabel('Minimum attitude change, $\Delta\theta_{min}$ [deg]', ...
            'Interpreter', 'latex', 'FontSize', 14);
        ylabel('$\frac{q_{pk}}{\Delta\theta_{pk}}$ [1/sec]', ...
            'Interpreter', 'latex', 'FontSize', 14); 
        
        savedStates = [t; U(3,:); [delta_e, delta_e(end)]; q; theta];
        nameFile = [num2str(distVal) 'degPitch.mat']
        filePath = fullfile('attitudeQuicknessSaves', nameFile);
        save(filePath, 'savedStates');

    elseif attitude == 2
        figure(20)
    
        subplot(3,1,1)
        plot(t(1:end-1), rad2deg(U(4,1:end-1)), 'LineWidth', 2);
        grid on; 
        legend('LatCyc [deg]', 'Interpreter', ...
            'latex', 'FontSize', 10, 'Location', 'SouthEast')
        ylabel('Input [deg]', 'Interpreter', 'latex', 'FontSize', 14)
        
        subplot(3,1,2)
        plot(t, rad2deg(p), 'LineWidth', 2); 
        grid on; 
        max_p_value = max(p);
        p_pk_deg = rad2deg(max_p_value)
        disp(['Maximum Value of p: ', num2str(max_p_value), 'Interpreter', 'latex']);
        hold on;
        yline(rad2deg(max_p_value), '--r');
        hold off;
        ylabel('p [deg/s]', 'Interpreter', 'latex', 'FontSize', 14)
        ylim([rad2deg(min(p))*1.2, rad2deg(max_p_value)*1.4]);
        xlim([0,t(end)]);
        
        subplot(3,1,3)
        plot(t, rad2deg(phi), 'LineWidth', 2); 
        grid on;
        hold on;
        plot(t(1:end-1), rad2deg(ref_step), 'm:', 'LineWidth', 1.5); 
        ylim([rad2deg(min(phi))*1.2, rad2deg(max(phi)*1.4)])
        hold off;
        ylabel('Attitude [deg]', 'Interpreter', 'latex', 'FontSize', 14)
        xlim([0,t(end)]);
        
        
        % Calculate rise time for phi
        percent_10 = 0.1; % 10%
        percent_90 = 0.9; % 90%
        steady_state_value_phi = mean(rad2deg(phi(end-100:end))); % Assuming the steady state is the mean of last 100 samples
        ten_percent_value_phi = percent_10 * steady_state_value_phi;
        ninety_percent_value_phi = percent_90 * steady_state_value_phi;
        idx_10_phi = find(rad2deg(phi) > ten_percent_value_phi, 1, 'first');
        idx_90_phi = find(rad2deg(phi) > ninety_percent_value_phi, 1, 'first');
        rise_time_phi = t(idx_90_phi) - t(idx_10_phi);
        
        % Find maximum value of theta
        max_value_phi = max(phi);
        phi_pk_deg = rad2deg(max_value_phi)
        
        % Find minimum value after the rise time for theta
        min_value_after_rise_phi = min(phi(idx_90_phi:end));
        rad2deg(min_value_after_rise_phi)
        
        % Plot horizontal red dashed lines for theta
        hold on;
        yline(rad2deg(max_value_phi), '--r');
        yline(rad2deg(min_value_after_rise_phi), '--r');
        hold off;
        legend('$\phi$ [deg]', 'Interpreter', 'latex', 'FontSize', 10, ...
            'Location', 'SouthEast')        
        % Display results for theta
        disp(['Rise Time for \phi: ', num2str(rise_time_phi)]);
        disp(['Maximum Value for \phi: ', num2str(max_value_phi)]);
        disp(['Minimum Value after Rise Time for \phi: ', num2str(min_value_after_rise_phi)]);
        
        % %% PLOT PITCH ATTITUDE QUICKNESS
        Qroll_5deg = max_p_value / (max_value_phi - phi(1))
        figure(21)
        % plot(rad2deg(min_value_after_rise_phi), Qroll_5deg, 'p', 'MarkerSize', 10)
        grid on; hold on;
        load('prollhq1.mat'); load('prollhq2.mat');
        plot(prollhq1(:,1), prollhq1(:,2), 'k--', 'LineWidth', 2);
        plot(prollhq2(:,1), prollhq2(:,2), 'k--', 'LineWidth', 2);
        % Add text labels for the different levels
        text(20, 2.5, 'Level 1', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
        text(10, 1.9, 'Level 2', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
        text(20, 1, 'Level 3', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
        
        
            otherQvals = [ 4.3211 3.8373 3.0694 2.4665 2.0404]';
            QvalsX = [ 9.1511 18.2739 27.2275 36.8204 45.4554];
            rad2deg(min_value_after_rise_phi)
            pp = plot(QvalsX(1:5), otherQvals(1:5), 'o', 'MarkerSize', 10);
            pp.MarkerFaceColor = [1 0.5 0];
    
        hold off;
        xlim([0, 50]); ylim([0, 5])
        xlabel('Minimum attitude change, $\Delta\phi_{min}$ [deg]', ...
            'Interpreter', 'latex', 'FontSize', 14);
        ylabel('$\frac{p_{pk}}{\Delta\phi_{pk}}$ [1/sec]', ...
            'Interpreter', 'latex', 'FontSize', 14); 
    else
        figure(20)

        subplot(3,1,1)
        plot(t(1:end-1), rad2deg(U(2,1:end-1)), 'LineWidth', 2); hold on;
        plot(t(1:end-1), rad2deg(delta_r), 'LineWidth', 2); hold off;
        grid on; 
        legend('DiffColl [deg]','Rudd [deg]', 'Interpreter', ...
            'latex', 'FontSize', 10, 'Location', 'SouthEast')
        ylabel('Input [deg]', 'Interpreter', 'latex', 'FontSize', 14)
        
        subplot(3,1,2)
        plot(t, rad2deg(r), 'LineWidth', 2); 
        grid on; 
        max_r_value = max(r);
        r_pk_deg = rad2deg(max_r_value)
        disp(['Maximum Value of r: ', num2str(max_r_value), 'Interpreter', 'latex']);
        hold on;
        yline(rad2deg(max_r_value), '--r');
        hold off;
        ylabel('r [deg/s]', 'Interpreter', 'latex', 'FontSize', 14)
        ylim([rad2deg(min(r))*1.2, rad2deg(max_r_value)*1.4]);
        xlim([0,t(end)]);
        
        subplot(3,1,3)
        plot(t, rad2deg(psi), 'LineWidth', 2); 
        grid on;
        hold on;
        plot(t(1:end-1), rad2deg(ref_step), 'm:', 'LineWidth', 1.5); 
        ylim([rad2deg(min(psi))*1.2, rad2deg(max(psi)*1.4)])
        hold off;
        ylabel('Attitude [deg]', 'Interpreter', 'latex', 'FontSize', 14)
        xlim([0,t(end)]);
        
        
        % Calculate rise time for theta
        percent_10 = 0.1; % 10%
        percent_90 = 0.9; % 90%
        steady_state_value_psi = mean(rad2deg(psi(end-100:end))); % Assuming the steady state is the mean of last 100 samples
        ten_percent_value_psi = percent_10 * steady_state_value_psi;
        ninety_percent_value_psi = percent_90 * steady_state_value_psi;
        idx_10_psi = find(rad2deg(phi) > ten_percent_value_psi, 1, 'first');
        idx_90_psi = find(rad2deg(psi) > ninety_percent_value_psi, 1, 'first');
        rise_time_psi = t(idx_90_psi) - t(idx_10_psi);
        
        % Find maximum value of theta
        max_value_psi = max(psi);
        psi_pk_deg = rad2deg(max_value_psi)
        
        % Find minimum value after the rise time for theta
        min_value_after_rise_psi = min(psi(idx_90_psi:end));
        rad2deg(min_value_after_rise_psi)
        
        % Plot horizontal red dashed lines for theta
        hold on;
        yline(rad2deg(max_value_psi), '--r');
        yline(rad2deg(min_value_after_rise_psi));
        hold off;
        legend('$\psi$ [deg]', 'Interpreter', 'latex', 'FontSize', 10, ...
            'Location', 'SouthEast') 
        
        % Display results for theta
        disp(['Rise Time for \psi: ', num2str(rise_time_psi)]);
        disp(['Maximum Value for \psi: ', num2str(max_value_psi)]);
        disp(['Minimum Value after Rise Time for \psi: ', num2str(min_value_after_rise_psi)]);
        
        % %% PLOT PITCH ATTITUDE QUICKNESS
        Qyaw_5deg = max_r_value / (max_value_psi - psi(1))
        figure(21)
        plot(rad2deg(min_value_after_rise_psi), Qyaw_5deg, 'p', 'MarkerSize', 10)
        grid on; hold on;
        load('ryawhq1.mat'); load('ryawhq2.mat'); 
        plot(ryawhq1(:,1), ryawhq1(:,2), 'k--', 'LineWidth', 2);
        plot(ryawhq2(:,1), ryawhq2(:,2), 'k--', 'LineWidth', 2);
        % Add text labels for the different levels
        text(20, 2.25, 'Level 1', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
        text(20, 1, 'Level 2', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
        text(20, 0.25, 'Level 3', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
        
            rad2deg(min_value_after_rise_psi)

            otherQvals = [ 2.4645 1.9070 1.6454 1.4781 1.3672 1.2865]';
            QvalsX = [ 8.7957 17.8957 27.1347 36.2954 44.7940 52.5691];
            pp = plot(QvalsX(1:6), otherQvals(1:6), 'o', 'MarkerSize', 10);
            pp.MarkerFaceColor = [1 0.5 0];
    
    
        hold off;
        xlim([0, 60]); ylim([0, 3])
        xlabel('Minimum attitude change, $\Delta\psi_{min}$ [deg]', ...
            'Interpreter', 'latex', 'FontSize', 14);
        ylabel('$\frac{r_{pk}}{\Delta\psi_{pk}}$ [1/sec]', ...
            'Interpreter', 'latex', 'FontSize', 14); 
    end


end

%%
%%%%%%%%%%%%%%%%%%%%%%
% figure(60)
% subplot(4,1,1)
% plot(t(1:end-1), rad2deg(U(3,1:end-1)), t(1:end-1), rad2deg(delta_e)); 
% grid on; legend('LonCyc [deg]', 'Elev [deg]', 'Interpreter', 'latex')
% 
% subplot(4,1,2)
% plot(t, u, t, v, t, w); grid on; 
% legend('u [m/s]', 'v [m/s]', 'w [m/s]', 'Interpreter', 'latex')
% 
% subplot(4,1,3)
% plot(t, p, t, rad2deg(q), t, r); grid on; 
% legend('p [rad/s]', 'q [deg/s]', 'r [rad/s]', 'Interpreter', 'latex')
% 
% subplot(4,1,4)
% plot(t, rad2deg(phi), t, rad2deg(theta), t, rad2deg(psi), 'LineWidth', 1.5); grid on;
% legend('$\phi_f$ [deg]', '$\theta_f$ [deg]', '$\psi_f$ [deg]', 'Interpreter', 'latex')

%%
% theta_deg = rad2deg(theta);
% custom_color = [0.93, 0.69, 0.13];  % A color that complements yellow and orange
% 
% theta_cmd_3211_deg = rad2deg([theta_cmd_3211; theta_cmd_3211(end)]);  % Padding to match lengths
% figure(61);
% plot(t(1:end-1), rad2deg(ref_3211), '--', 'LineWidth', 1.5);  hold on;
% plot(t(1:end-1), rad2deg(theta_cmd_3211), '-.', 'LineWidth', 2);
% plot(t, theta_deg, 'LineWidth', 2);
% 
% 
% fill([t fliplr(t)], [theta_deg fliplr(theta_cmd_3211_deg')], custom_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% 
% hold off;
% grid on; 
% ylabel('$\theta$ [deg]', 'Interpreter', 'latex', 'FontSize', 14)
% xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14)
% legend('$\theta_{cmd}$','$\theta_{c}$','$\theta$', 'Interpreter', 'latex', 'FontSize', 13)
% 

%%
% figure(62)
% subplot(2,1,1)
% plot(t(1:end-1), V_x_ref, '--', 'LineWidth', 1.5); grid on; hold on;
% plot(t(1:end-1), V_x_cmd, '-.', 'LineWidth', 1.5); hold off
% ylabel('$V_x$', 'Interpreter', 'latex', 'FontSize', 16);
% % xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
% legend('$V_{x_{cmd}}$', '$V_{x_c}$', 'Interpreter', 'latex', ...
%     'FontSize', 16, 'Location', 'Best')
% subplot(2,1,2)
% plot(t(1:end-1), V_z_ref, '--', 'LineWidth', 1.5); grid on; hold on;
% plot(t(1:end-1), V_z_cmd, '-.', 'LineWidth', 1.5); hold off
% ylabel('$V_z$', 'Interpreter', 'latex', 'FontSize', 16);
% xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
% legend('$V_{z_{cmd}}$', '$V_{z_c}$', 'Interpreter', 'latex', ...
%     'FontSize', 16, 'Location', 'SouthEast')

%%
% f = figure(81);
% 
% subplot(3,2,1)
% set(gca, 'FontSize', 10);
% plot(t(1:end-1), V_x_cmd, '--', 'LineWidth', 1.5); hold on;
% plot(t, u, 'LineWidth', 1.5); hold off; grid on;
% ylabel('$V_x$ [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
% ylim([25,45]); xlim([0,t(end)])
% 
% subplot(3,2,2)
% set(gca, 'FontSize', 10);
% plot(t, rad2deg(phi_cmd), '--', 'LineWidth', 1.5); hold on;
% plot(t, rad2deg(phi), 'LineWidth', 1.5); hold off; grid on;
% ylabel('$\phi [^{\circ}]$', 'Interpreter', 'latex', 'FontSize', 14);
% ylim([-2,2]); xlim([0,t(end)])
% 
% subplot(3,2,3)
% set(gca, 'FontSize', 10);
% plot(t(1:end-1), V_y_cmd, '--', 'LineWidth', 1.5); hold on;
% plot(t, v, 'LineWidth', 1.5); hold off; grid on;
% ylabel('$V_y$ [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
% ylim([-1,1]);xlim([0,t(end)])
% 
% subplot(3,2,4)
% set(gca, 'FontSize', 10);
% plot(t, rad2deg(theta_cmd), '--', 'LineWidth', 1.5); hold on;
% plot(t, rad2deg(theta), 'LineWidth', 1.5); hold off; grid on;
% ylabel('$\theta [^{\circ}]$', 'Interpreter', 'latex', 'FontSize', 14)
% ylim([-10,10]);xlim([0,t(end)])
% 
% subplot(3,2,5)
% set(gca, 'FontSize', 10);
% plot(t(1:end-1), -V_z_cmd, '--', 'LineWidth', 1.5); hold on;
% plot(t, -w, 'LineWidth', 1.5); hold off; grid on;
% ylabel('$V_z$ [m/s]', 'Interpreter', 'latex', 'FontSize', 14)
% xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14)
% ylim([-2.5,2.5]);xlim([0,t(end)])
% 
% subplot(3,2,6)
% set(gca, 'FontSize', 11);
% plot(t(1:end-1), rad2deg(psi_cmd), '--', 'LineWidth', 1.5); hold on;
% plot(t, rad2deg(psi), 'LineWidth', 1.5); hold off; grid on;
% ylabel('$\psi [^{\circ}]$', 'Interpreter', 'latex', 'FontSize', 14)
% xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14)
% ylim([-2, 2]);xlim([0,t(end)])
% 
% f.Position = [600 300 700 500];
% 
% % Create an invisible axis to place the legend
% hL = legend('Command', 'EMF Controller', 'Orientation', 'horizontal');
% hL.Position = [0.5 0.97 0 0];
% hL.FontSize = 13;
% hL.Interpreter = 'latex';

%% 
% f = figure(82);
% 
% subplot(4,2,1)
% plot(t(1:end-1), rad2deg(U(1,1:end-1)), 'Color', '#D95319', 'LineWidth', 2); 
% grid on; xlim([0,t(end)]); ylim([7, 15])
% ylabel('$\theta_0 [^{\circ}]$', 'Interpreter', 'latex', 'FontSize', 18)
% 
% subplot(4,2,2)
% plot(t(1:end-1), rad2deg(U(6,1:end-1)), 'Color', '#D95319', 'LineWidth', 2); 
% grid on; xlim([0,t(end)]); ylim([10,35])
% ylabel('$\theta_p [^{\circ}]$', 'Interpreter', 'latex', 'FontSize', 18)
% 
% subplot(4,2,3)
% plot(t(1:end-1), rad2deg(U(2,1:end-1)), 'Color', '#D95319', 'LineWidth', 2); 
% grid on; xlim([0,t(end)]); ylim([-1,1])
% ylabel('$\theta_d [^{\circ}]$', 'Interpreter', 'latex', 'FontSize', 18)
% 
% subplot(4,2,6)
% plot(t(1:end-2), rad2deg(delta_e(1:end-1)), 'Color', '#D95319', 'LineWidth', 2); 
% grid on; xlim([0,t(end)]); ylim([-1,1])
% ylabel('$\delta_e [^{\circ}]$', 'Interpreter', 'latex', 'FontSize', 18)
% 
% subplot(4,2,4)
% plot(t(1:end-1), rad2deg(U(3,1:end-1)), 'Color', '#D95319', 'LineWidth', 2); 
% grid on; xlim([0,t(end)]);
% ylabel('$\theta_{1s} [^{\circ}]$', 'Interpreter', 'latex', 'FontSize', 18)
% 
% subplot(4,2,7)
% plot(t(1:end-2), rad2deg(delta_r(1:end-1)), 'Color', '#D95319', 'LineWidth', 2); 
% grid on; xlim([0,t(end)]); ylim([-1,1])
% ylabel('$\delta_r [^{\circ}]$', 'Interpreter', 'latex', 'FontSize', 18)
% xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14)
% 
% subplot(4,2,5)
% plot(t(1:end-1), rad2deg(U(4,1:end-1)), 'Color', '#D95319', 'LineWidth', 2); 
% grid on; xlim([0,t(end)]); ylim([-1,1])
% ylabel('$\theta_{1c} [^{\circ}]$', 'Interpreter', 'latex', 'FontSize', 18)
% xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14)
% 
% f.Position = [600 300 700 700];

% %%
% figure(1)
% plot(t(1:end-1), rad2deg(ref_3211), '--', 'LineWidth', 2); grid on;
% ylim([-4,4]); 
% ylabel('Commanded Attitude Angle [$^{\circ}$]', 'Interpreter', 'latex', 'FontSize', 14)
% xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14')

%%
% figure(2)
% 
% if subSelection3211 == 1
%     subplot(3,1,1)
%     plot(t, zeros(length(phi)), '--', 'LineWidth', 2, 'Color', '#0072BD'); grid on; hold on;
%     plot(t, zeros(length(phi)), '-.', 'LineWidth', 2, 'Color', '#D95319');
%     plot(t, rad2deg(phi), '-', 'LineWidth', 2, 'Color', '#EDB120'); hold off;
%     ylabel('$\phi$ [$^{\circ}$]', 'Interpreter', 'latex', 'FontSize', 14)
%     ylim([-5,5]); xlim([0,10])
%     subplot(3,1,2)
%     plot(t(1:end-1), rad2deg(ref_3211), '--', 'LineWidth', 2); grid on; hold on;
%     plot(t(1:end-1), rad2deg(theta_cmd_3211), '-.', 'LineWidth', 2);
%     plot(t, rad2deg(theta), '-', 'LineWidth', 2); hold off;
%     ylabel('$\theta$ [$^{\circ}$]', 'Interpreter', 'latex', 'FontSize', 14)
%     ylim([-5,5]); xlim([0,10])
%     subplot(3,1,3)
%     plot(t, zeros(length(psi)), '--', 'LineWidth', 2, 'Color', '#0072BD'); grid on; hold on;
%     plot(t, zeros(length(psi)), '-.', 'LineWidth', 2, 'Color', '#D95319');
%     plot(t, rad2deg(psi), '-', 'LineWidth', 2, 'Color', '#EDB120'); hold off;
%     ylabel('$\psi$ [$^{\circ}$]', 'Interpreter', 'latex', 'FontSize', 14)
%     ylim([-5,5]); xlim([0,10])
%     xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14')
% elseif subSelection3211 == 2
%     subplot(3,1,1)
%     plot(t(1:end-1), rad2deg(ref_3211), '--', 'LineWidth', 2, 'Color', '#0072BD'); grid on; hold on;
%     plot(t(1:end-1), rad2deg(theta_cmd_3211), '-.', 'LineWidth', 2, 'Color', '#D95319');
%     plot(t, rad2deg(phi), '-', 'LineWidth', 2, 'Color', '#EDB120'); hold off;
%     ylabel('$\phi$ [$^{\circ}$]', 'Interpreter', 'latex', 'FontSize', 14)
%     ylim([-5,5]); xlim([0,10])
%     subplot(3,1,2)
%     plot(t, zeros(length(theta)), '--', 'LineWidth', 2, 'Color', '#0072BD'); grid on; hold on;
%     plot(t, zeros(length(theta)), '-.', 'LineWidth', 2, 'Color', '#D95319');
%     plot(t, rad2deg(theta), '-', 'LineWidth', 2, 'Color', '#EDB120'); hold off;
%     ylabel('$\theta$ [$^{\circ}$]', 'Interpreter', 'latex', 'FontSize', 14)
%     ylim([-5,5]); xlim([0,10])
%     subplot(3,1,3)
%     plot(t, zeros(length(psi)), '--', 'LineWidth', 2, 'Color', '#0072BD'); grid on; hold on;
%     plot(t, zeros(length(psi)), '-.', 'LineWidth', 2, 'Color', '#D95319');
%     plot(t, rad2deg(psi), '-', 'LineWidth', 2, 'Color', '#EDB120'); hold off;
%     ylabel('$\psi$ [$^{\circ}$]', 'Interpreter', 'latex', 'FontSize', 14)
%     ylim([-5,5]); 
%     xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14')
% else
%     subplot(3,1,1)
%     plot(t, zeros(length(phi)), '--', 'LineWidth', 2, 'Color', '#0072BD'); grid on; hold on;
%     plot(t, zeros(length(phi)), '-.', 'LineWidth', 2, 'Color', '#D95319');
%     plot(t, rad2deg(phi), '-', 'LineWidth', 2, 'Color', '#EDB120'); hold off;
%     ylabel('$\phi$ [$^{\circ}$]', 'Interpreter', 'latex', 'FontSize', 14)
%     ylim([-5,5]); xlim([0,10])
%     subplot(3,1,2)
%     plot(t, zeros(length(theta)), '--', 'LineWidth', 2, 'Color', '#0072BD'); grid on; hold on;
%     plot(t, zeros(length(theta)), '-.', 'LineWidth', 2, 'Color', '#D95319');
%     plot(t, rad2deg(theta), '-', 'LineWidth', 2, 'Color', '#EDB120'); hold off;
%     ylabel('$\theta$ [$^{\circ}$]', 'Interpreter', 'latex', 'FontSize', 14)
%     ylim([-5,5]); xlim([0,10])
%     subplot(3,1,3)
%     plot(t(1:end-1), rad2deg(ref_3211), '--', 'LineWidth', 2, 'Color', '#0072BD'); grid on; hold on;
%     plot(t(1:end-1), rad2deg(theta_cmd_3211), '-.', 'LineWidth', 2, 'Color', '#D95319');
%     plot(t, rad2deg(psi), '-', 'LineWidth', 2, 'Color', '#EDB120'); hold off;
%     ylabel('$\psi$ [$^{\circ}$]', 'Interpreter', 'latex', 'FontSize', 14)
%     ylim([-5,5]); xlim([0,10])
%     xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14')
% end
% 
% 
% 
% f = figure(3);
% 
% subplot(4,1,1)
% plot(t(1:end-1), rad2deg(U(1,1:end-1)), 'Color', '#0072BD', 'LineWidth', 2); 
% grid on; xlim([0,t(end)]); 
% ylabel('Vertical', 'Interpreter', 'latex', 'FontSize', 14)
% legend('$\theta_{0}$', 'Interpreter', 'latex', 'FontSize', 12)
% 
% subplot(4,1,2)
% plot(t(1:end-1), rad2deg(U(4,1:end-1)), 'Color', '#0072BD', 'LineWidth', 2); 
% grid on; xlim([0,t(end)]); 
% xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel('Lateral', 'Interpreter', 'latex', 'FontSize', 14)
% legend('$\theta_{1c}$', 'Interpreter', 'latex', 'FontSize', 12)
% 
% subplot(4,1,3)
% plot(t(1:end-1), rad2deg(U(3,1:end-1)), '--', 'Color', '#D95319', 'LineWidth', 2); hold on;
% plot(t(1:end-2), rad2deg(delta_e(1:end-1)), 'Color', '#0072BD', 'LineWidth', 2); hold off;
% grid on; xlim([0,t(end)]); 
% ylabel('Longitudinal', 'Interpreter', 'latex', 'FontSize', 14)
% legend('$\theta_{1s}$', '$\delta_e$', 'Interpreter', 'latex', 'FontSize', 12)
% 
% 
% subplot(4,1,4)
% plot(t(1:end-1), rad2deg(U(2,1:end-1)), 'Color', '#0072BD', 'LineWidth', 2); hold on;
% plot(t(1:end-2), rad2deg(delta_r(1:end-1)), '--', 'Color', '#D95319', 'LineWidth', 2); hold off;
% grid on; xlim([0,t(end)]); 
% ylabel('Directional', 'Interpreter', 'latex', 'FontSize', 14)
% legend('$\theta_{d}$', '$\delta_r$', 'Interpreter', 'latex', 'FontSize', 12)
% 
% xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14)
% 
% 
% % f.Position = [600 300 600 600];




























