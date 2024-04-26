% set the rrim speed 'v_kts'  

clc; close all; clear;

v_kts = 55;

load('TrimLinSave.mat'); load('V_transit_save.mat');

coaxial_heli_parameters;

% ii     [1  2   3   4   5   6   7   8   9  10  11  12  13  14  
%   15  16  17  18  19  20   21   22   23   24   25   26]
kias_s = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, ...
    70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125];

v_kts = min(max(v_kts, kias_s(1)), kias_s(end));

vv = kias_s(max(find(v_kts >= kias_s)));

ii = find(vv == kias_s);

tlr = TrimLinSave(ii, :);

AA = tlr{1};
BB = tlr{2};
CC = tlr{3};
DD = tlr{4};

% Find Correct Stability and Control Derivatives according to speed var ii
M_theta_c = BB(3,2);
M_delta_e = BB(3,4);

% CA Stuff
V_l_transit = V_transit_save(1);
V_u_transit = V_transit_save(2);

% Error Compensation Gains
Kp_pitch = -50;
Ki_pitch = -50;
Kd_pitch = 33.3;












