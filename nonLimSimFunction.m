%% Nonlinear Simulation

clear;clc;close all;

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


tEnd = 14;
simN = 400;
u(1)=0;
v(1)=0;
w(1)=0;


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
V_l_transit = V_transit_save(1);
V_u_transit = V_transit_save(2)-20;

Kp_theta = 80;
Ki_theta = 0;
Kp_q = -25;
Kd_theta = 0;
r_theta = 1;

Kp_phi = 120;
Ki_phi = 0.8;
Kp_p = -2.8;
Kd_phi = 0;
r_phi = 1;

Kp_r = 0.01;
Ki_r = 0;
r_r = 1;

Kp_x = -0.9 ;
Ki_x = -0;

Kp_z = 5;
Ki_z = 0.1;
Kd_z = 0.0;

Kp_y = -deg2rad(1.5);
Ki_y = 0;

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

    % Keep certain inputs constant for t
    theta_cdiff(i) = 0;

    % Select Derivatives for Correct Speed
    [~, idx_speed] = min(abs(V_vals - V(i)));

    L_theta_c = TrimLinSave_6dof{idx_speed,2}(4,4);
    M_theta_s = TrimLinSave_6dof{idx_speed,2}(5,3);
    M_delta_e = TrimLinSave_6dof{idx_speed,2}(5,7);
    N_theta_d = TrimLinSave_6dof{idx_speed,2}(6,2);
    N_delta_r = TrimLinSave_6dof{idx_speed,2}(6,8);
    X_theta_p = TrimLinSave_6dof{idx_speed,2}(1,6);
    Z_theta_0 = TrimLinSave_6dof{idx_speed,2}(3,1);

    % Inner Loop ACAH Controllers
    phi_cmd_3211(i) = 0;
    error_phi(i) = phi_cmd_3211(i) - phi(i);
    error_theta(i) = theta_cmd_3211(i) - theta(i);
    r_cmd_3211(i) = 0;
    error_r(i) = r_cmd_3211(i) - r(i);



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

    integral_r = integral_r + error_r(i) * dt;
    gamma_r(i) = r_r * (Kp_r*error_r(i) + Ki_r * integral_r );

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
            ref_3211(i) = deg2rad(1);
        elseif t(i) >= 4 && t(i) < 6 
            ref_3211(i) = deg2rad(-1);
        elseif t(i) >= 6 && t(i) < 7 
            ref_3211(i) = deg2rad(1);
        elseif t(i) >= 7 && t(i) < 8 
            ref_3211(i) = deg2rad(-1);
        else t(i) >= 8;
            ref_3211(i) = 0;
        end
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
        B_x = [1, -X_theta_p];
        W_x = [1/rho_theta_x_cmd(i), 0; 0, 1/rho_theta_p(i)];
        U_x(:,i) = inv(W_x)*transpose(B_x)* ...
            inv(B_x*inv(W_x)*transpose(B_x))*gamma_x(i);
        theta_x_cmd(i) = U_x(1,i);
        U(6,i) = U_x(2,i) + trim_saved_6dof(6,idx_trimspeed);

        % Height CA
        B_z = [Z_theta_0, 1];
        W_z = [1/rho_theta_0(i), 0; 0, 1/rho_theta_z_cmd(i)];
        U_z(:,i) = inv(W_z)*transpose(B_z)* ...
            inv(B_z*inv(W_z)*transpose(B_z))*gamma_z(i);
        U(1,i) = U_z(1,i) + trim_saved_6dof(1,idx_trimspeed);
        theta_z_cmd(i) = U_z(2,i);

        % Pitch CA
        B_theta = [M_theta_s*0.95, M_delta_e*0.95];
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
        
        % Elev
        elev_rate(i) = (delta_e(i)-delta_e(i-1))/dt;
        if abs(elev_rate(i)) > 28.8
            delta_e(i) = sign(elev_rate(i))*28.8*dt + delta_e(i-1);
        end
        elev_rate(i) = (delta_e(i)-delta_e(i-1))/dt;

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


%%
% EMF Controller Performance (Q) Calculation
if ManoeuvreSelection == 1
    err_theta = sum(abs(theta(1:end-1) - theta_cmd_3211'));
    % err_phi = sum(abs(phi(1:end-1) - phi_cmd_3211'));
    err_phi = sum(abs(phi(1:end-1) - 0));
    err_phi = sum(abs(phi(1:end-1) - 0));
else
    err_theta = sum(abs(theta - theta_cmd));
    err_phi = sum(abs(phi - phi_cmd));
end

V = 1; %deg/in
err_rel_theta = 1/(V*tEnd)*err_theta;
Q_theta = 1-err_rel_theta
err_rel_phi = 1/(V*tEnd)*err_phi;
Q_phi = 1-err_rel_phi

% err_u = sum(abs(u(1:end-1) - V_x_cmd));
% V_u = 1;
% err_rel_u = 1/(V_u*tEnd)*err_u;
% Q_u = 1-err_rel_u

% Save Q value in the list
% Q_values = [Q_values, Q_theta];
% r_values = [r_values, r_theta];

% end

% disp(Q_values)

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
figure(81)
subplot(4,2,1)
plot(t(1:end-1), V_x_ref, t(1:end-1), V_x_cmd, t, u); grid on;
ylabel('V_x [m/s]'); %ylim([-10;10])
legend('Reference Signal', 'Model Signal')
subplot(4,2,3)
plot(t(1:end-1), V_y_cmd, t, v); grid on;
ylabel('V_y [m/s]')
subplot(4,2,5)
plot(t(1:end-1), V_z_ref, t(1:end-1), V_z_cmd, t, w); grid on;
ylabel('V_z [m/s]')
xlabel('Time [s]')

subplot(4,2,2)
plot(t, rad2deg(phi), t, rad2deg(phi_cmd)); grid on;
ylabel('\phi [deg]')
subplot(4,2,4)
if ManoeuvreSelection == 1
    plot(t(1:end-1), rad2deg(theta_cmd_3211), t, rad2deg(theta)); grid on;
else
    plot(t, rad2deg(theta_cmd), t, rad2deg(theta)); grid on;
end
ylabel('\theta [deg]')
legend('\theta [deg]', '\theta_{cmd} [deg]')
subplot(4,2,6)
plot(t(1:end-1), psi_cmd, t, rad2deg(psi)); grid on;
ylabel('\psi [deg]')
subplot(4,2,8)
plot(t(1:end-1), r_cmd, t, rad2deg(r)); grid on;
ylabel('r [deg/s]')
xlabel('Time [s]')
subplot(4,2,7)
plot(t, rad2deg(p))
% 
figure(82)
subplot(4,2,1)
plot(t(1:end-1), rad2deg(U(1,1:end-1))); legend('Coll'); grid on;
subplot(4,2,3)
plot(t(1:end-1), rad2deg(U(2,1:end-1))); legend('DiffColl'); grid on;

subplot(4,2,5)
plot(t(1:end-1), rad2deg(U(3,1:end-1))); legend('LonCyc'); grid on;
subplot(4,2,7)
plot(t(1:end-1), rad2deg(U(4,1:end-1))); legend('LatCyc'); grid on;

subplot(4,2,2)
plot(t(1:end-1), rad2deg(U(6,1:end-1))); legend('PropColl');grid on;

subplot(4,2,4)
plot(t(1:end-2), rad2deg(delta_e(1:end-1))); legend('Elev'); grid on;
subplot(4,2,6)
plot(t(1:end-2), rad2deg(delta_r(1:end-1))); legend('Rudder'); grid on;

% figure(85); 
% if ManoeuvreSelection == 1
%     plot(t, rad2deg(theta), t(1:end-1), rad2deg(theta_cmd_3211)); grid on;
% else
%     plot(t, rad2deg(theta), t, rad2deg(theta_cmd)); grid on;
% end
% ylabel('\theta [deg]')
% xlabel('Time [s]')
% legend('\theta [deg]', '\theta_{cmd} [deg]')

%%
% figure(86)
% if ManoeuvreSelection == 1
%     subplot(3,1,1)
%     plot(t(1:end-1), rad2deg(ref_3211), t(1:end-1), rad2deg(theta_cmd_3211), ...
%         t, rad2deg(theta)); grid on; ylabel('\theta')
%     subplot(3,1,2)
%     plot(t(1:end-1), rad2deg(ref_3211), t(1:end-1), rad2deg(phi_cmd_3211), ...
%         t, rad2deg(phi)); grid on; ylabel('\phi')
%     subplot(3,1,3)
%     plot(t(1:end-1), rad2deg(ref_3211), t(1:end-1), rad2deg(r_cmd_3211), ...
%         t, rad2deg(r)); grid on; ylabel('r')
% end
% 
% figure(87)
% subplot(3,1,1)
% plot(t,rad2deg(U(3,:)), t(1:end-1), rad2deg(delta_e)); legend('LonCyc', 'Elev'); grid on;
% subplot(3,1,2)
% plot(t,rad2deg(U(4,:))); legend('LatCyc'); grid on;
% subplot(3,1,3)
% plot(t,rad2deg(U(2,:)), t(1:end-1),rad2deg(delta_r)); legend('Differential Collective', 'Rudder'); grid on;


