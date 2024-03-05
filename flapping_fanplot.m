% %%%%%%%%%% Flapping Angle Fan Plot %%%%%%%%%%
% 
% coaxial_heli_parameters;
% 
% % Iblade = mblade*R^3*(1-e)^3/3;
% mblade = 30;                    % [kg]
% Iblade = 1921;                  % moment of inertia [kg m^2]
% 
% Kb = 460000;                    % Nm
% relOffset = 0.06;
% 
% w_nr = sqrt(Kb/Iblade);
% v_nr = w_nr / Omega;
% southwell = 1 + 3/2 * relOffset/(1-relOffset);
% 
% OmegaRange = [0:0.1:Omega*2];
% 
% % Southwell Approximation (w_f^2 = w_nr^2 + SouthwellCoefficient*Omega^2)
% w = sqrt(w_nr^2 + southwell.*OmegaRange.^2);
% v = w/Omega;
% 
% % Using equation from yuqing 
% Mblade = mblade*R^2/3;               % mass moment of inertia approximation of a thin rod
% wf = sqrt(1 + relOffset*Mblade/Iblade + Kb./(Iblade.*OmegaRange.^2)).*OmegaRange;
% vf = wf/Omega;
% 
% 
% figure(1)
% plot(OmegaRange, w/(2*pi), OmegaRange, wf/(2*pi))
% ylabel("Flapping Frequency [Hz]"); xlabel("Rotational Speed [rad/s]")
% ylim([0,30]);
% hold on; plot([35, 35], [0, 30]); hold off; 
% legend("Southwell Method", "Literature Method", "Nominal RPM")
% 
% figure(2)
% plot(OmegaRange, v/(2*pi), OmegaRange, vf/(2*pi))
% ylabel("Flapping Frequency [Hz]"); xlabel("Rotational Speed [rad/s]")
% ylim([0,1.2]);
% hold on; plot([35, 35], [0, 30]); hold off; 
% legend("Southwell Method", "Literature Method", "Nominal RPM")

coaxial_heli_parameters;

% Iblade = mblade*R^3*(1-e)^3/3;
mblade = 30;                    % [kg]
Iblade = 1921;                  % moment of inertia [kg m^2]

Kb = 460000;                    % Nm
Omega = 35;                     % Nominal RPM
R = 1.5;                        % Rotor radius [m]

% Define range of relOffset values
relOffset_range = 0.02:0.01:0.10;

OmegaRange = 0.1:0.1:Omega*1.2;

% Preallocate matrix to store results
w_results = zeros(length(OmegaRange), length(relOffset_range));

% Loop over each relOffset value

for i = 1:length(relOffset_range)
    relOffset = relOffset_range(i);

    w_nr = sqrt(Kb/Iblade);
    v_nr = w_nr / Omega;
    southwell = 1 + 3/2 * relOffset/(1-relOffset);

    % Calculate w for current relOffset
    w = sqrt(w_nr^2 + southwell.*OmegaRange.^2);

    % Store results
    w_results(:, i) = w;

    v_results(:, i) = w/Omega;
end

% Plot the results
figure(1);
plot(OmegaRange, v_results/(2*pi)); grid on;
ylabel("Flapping Frequency [Hz]");
xlabel("Rotational Speed [rad/s]");
ylim([0, 0.4]);
hold on;
plot([Omega, Omega], [0, 30], 'k--'); % Plot nominal RPM as vertical line
hold off;
title('Flapping Frequency vs. Rotational Speed for Different relOffset Values');

%%
coaxial_heli_parameters;

% Iblade = mblade*R^3*(1-e)^3/3;
mblade = 30;                    % [kg]
Iblade = 1921;                  % moment of inertia [kg m^2]

Omega = 35;                     % Nominal RPM
R = 1.5;                        % Rotor radius [m]
relOffset = 0.06;               % Single value of relOffset

% Define range of Kb values
Kb_range = 100000:50000:1000000;

% Preallocate matrix to store results
w_results = zeros(length(OmegaRange), length(Kb_range));

% Loop over each Kb value
for i = 1:length(Kb_range)
    Kb = Kb_range(i);

    w_nr = sqrt(Kb/Iblade);
    v_nr = w_nr / Omega;
    southwell = 1 + 3/2 * relOffset/(1-relOffset);

    % Calculate w for current Kb
    w = sqrt(w_nr^2 + southwell.*OmegaRange.^2);

    % Store results
    w_results(:, i) = w;

    v_results(:, i) = w/Omega;
end

% Plot the results
figure(2);
plot(OmegaRange, w_results/(2*pi));
ylabel("Flapping Frequency [Hz]");
xlabel("Rotational Speed [rad/s]");
ylim([0, 80]); grid on;
hold on;
plot([Omega, Omega], [0, 30], 'k--'); % Plot nominal RPM as vertical line
hold off;
title('Flapping Frequency vs. Rotational Speed for Different Kb Values');

