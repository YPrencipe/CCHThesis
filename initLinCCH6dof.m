% set the rrim speed 'v_kts'  

clc; close all; clear;

trimvelocity = 65; % m/s

% filePath = fullfile('LinSimSaves', 'xxms_xx.mat');
% save(filePath, 'out');

load('TrimLinSave_6dof.mat'); load('V_transit_save.mat');
load('trim_saved_6dof.mat'); load('ytrim.mat');
load('stabDerivs.mat'); load('contrDerivs.mat');

coaxial_heli_parameters;

% ii     [1  2   3   4   5   6   7   8   9  10  11  12  13  14  
%   15  16  17  18  19  20   21   22   23   24   25   26]
kias_s = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, ...
    70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125];

trimvelocity = min(max(trimvelocity, kias_s(1)), kias_s(end));

vv = kias_s(max(find(trimvelocity >= kias_s)));

ii = find(vv == kias_s);

tlr = TrimLinSave_6dof(ii, :);

AA = tlr{1};
BB = tlr{2};
CC = tlr{3};
DD = tlr{4};
ytrim = ytrim(:,ii);

% Find Correct Stability and Control Derivatives according to speed var ii
M_theta_s = BB(5,3);
M_delta_e = BB(5,7);
L_theta_c = BB(4,4);
N_theta_d = BB(6,2);
N_delta_r = BB(6,8);
X_theta_p = BB(1,6);
Z_theta_0 = BB(3,1);

% CA Stuff
V_l_transit = V_transit_save(1);
V_u_transit = V_transit_save(2);

% TrimValues
CollTrim = trim_saved_6dof(1,ii);
DiffLatCycTrim = trim_saved_6dof(5,ii);
PushPropTrim = trim_saved_6dof(6,ii);
















