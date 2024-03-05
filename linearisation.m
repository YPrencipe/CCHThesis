coaxial_heli_parameters;
try3_trim_3dof;

% remember: x_k(:,1) = [theta_0(1); theta_c(1); theta_p(1); lambda_0_u(1); lambda_p(1)];
%           theta_f_trim = deg2rad([4.8, 4.8, 4.5, 4.3, 4, 3.8, 3.7, 3.5, 3.3, 3, 2.5, 2, 1.5, 1, ...
%    0.5, 0 ,-0.5, -1 ,-1.5 ,-2.5, -2.5]); --> pre-scheduled fuselage pitch
%    for trim, according to trim algorithm

% We will make a state-space representation at every trim point
X_0 = -mass*g*sin(theta_f_trim);    % not correct in reader (depending on pitch defintion)
Z_0 = mass*g*cos(theta_f_trim);     % idem
M_0 = 0;


