clear;clc;close all;

coaxial_heli_parameters;
load('TrimLinSave_6dof.mat'); 
load('NLR_lin.mat');

for i = 1:length(TrimLinSave_6dof)
    X_u(i) = TrimLinSave_6dof{i, 1}(1, 1);
    
    X_w(i) = TrimLinSave_6dof{i, 1}(1, 3);
    X_q(i) = TrimLinSave_6dof{i, 1}(1, 5);

    Z_u(i) = TrimLinSave_6dof{i, 1}(3, 1);
    Z_w(i) = TrimLinSave_6dof{i, 1}(3, 3);
    Z_q(i) = TrimLinSave_6dof{i, 1}(3, 5);

    M_u(i) = TrimLinSave_6dof{i, 1}(5, 1);
    M_w(i) = TrimLinSave_6dof{i, 1}(5, 3);
    M_q(i) = TrimLinSave_6dof{i, 1}(5, 5);

    X_theta_0(i) = TrimLinSave_6dof{i, 2}(1, 1);
    
    Z_theta_0(i) = TrimLinSave_6dof{i, 2}(3, 1);

end

for i = 1:length(NLR_lin)
    V_vals_kts_NLR(i) = (i-1)*10*1.944;

    X_u_NLR(i) = NLR_lin{i, 1}(1, 1);
    X_w_NLR(i) = NLR_lin{i, 1}(1, 3);
    X_q_NLR(i) = NLR_lin{i, 1}(1, 5);

    Z_u_NLR(i) = NLR_lin{i, 1}(3, 1);
    Z_w_NLR(i) = NLR_lin{i, 1}(3, 3);
    Z_q_NLR(i) = NLR_lin{i, 1}(3, 5);

    M_u_NLR(i) = NLR_lin{i, 1}(5, 1);
    M_w_NLR(i) = NLR_lin{i, 1}(5, 3);
    M_q_NLR(i) = NLR_lin{i, 1}(5, 5);

    X_theta_0_u_NLR(i) = NLR_lin{i, 2}(1, 1);
    X_theta_0_l_NLR(i) = NLR_lin{i, 2}(1, 4);

    Z_theta_0_u_NLR(i) = NLR_lin{i, 2}(3, 1);
    Z_theta_0_l_NLR(i) = NLR_lin{i, 2}(3, 4);

end

figure(2)
% plot(V_vals_kts, X_theta_0, V_vals_kts_NLR, X_theta_0_u_NLR, V_vals_kts_NLR, X_theta_0_l_NLR); 
plot(V_vals_kts, Z_theta_0, V_vals_kts_NLR, Z_theta_0_u_NLR, V_vals_kts_NLR, Z_theta_0_l_NLR); 


% figure(1)
% subplot(3,3,1)
% plot(V_vals_kts, X_u, V_vals_kts_NLR, X_u_NLR); grid on; xlabel('V_f [kts]'); ylabel('X_u'); legend('X_u')
% subplot(3,3,2)
% plot(V_vals_kts, X_w, V_vals_kts_NLR, X_w_NLR); grid on; xlabel('V_f [kts]'); ylabel('X_w'); legend('X_w')
% subplot(3,3,3)
% plot(V_vals_kts, X_q, V_vals_kts_NLR, X_q_NLR); grid on; xlabel('V_f [kts]'); ylabel('X_q'); legend('X_q')
% 
% subplot(3,3,4)
% plot(V_vals_kts, Z_u, V_vals_kts_NLR, Z_u_NLR); grid on; xlabel('V_f [kts]'); ylabel('Z_u'); legend('Z_u')
% subplot(3,3,5)
% plot(V_vals_kts, Z_w, V_vals_kts_NLR, Z_w_NLR); grid on; xlabel('V_f [kts]'); ylabel('Z_w'); legend('Z_w')
% subplot(3,3,6)
% plot(V_vals_kts, Z_q, V_vals_kts_NLR, Z_q_NLR); grid on; xlabel('V_f [kts]'); ylabel('Z_q'); legend('Z_q')
% 
% subplot(3,3,7)
% plot(V_vals_kts, M_u, V_vals_kts_NLR, M_u_NLR); grid on; xlabel('V_f [kts]'); ylabel('M_u'); legend('M_u')
% subplot(3,3,8)
% plot(V_vals_kts, M_w, V_vals_kts_NLR, M_w_NLR); grid on; xlabel('V_f [kts]'); ylabel('M_w'); legend('M_w')
% subplot(3,3,9)
% plot(V_vals_kts, M_q, V_vals_kts_NLR, M_q_NLR); grid on; xlabel('V_f [kts]'); ylabel('M_q'); legend('M_q')
