%% verify_dynamics.m

clc;
clear;
close all;

fprintf('--- Verifying dynamics using truth trajectory ---\n');

%% Load truth trajectory
load('Project2_Prob2_truth_traj_50days.mat'); % loads Tt_50 and Xt_50

% Xt_50: each row is [x, y, z, vx, vy, vz, CR, STM (7x7 flattened)]

% Time vector
t_truth = Tt_50; % seconds

% Initial state and STM
X0_truth = Xt_50(1,:)'; % 56x1 initial vector

%% Set up parameters
params = struct();
params.mu_earth = 398600.4415;           % km^3/s^2
params.mu_sun = 132712440017.987;         % km^3/s^2
params.AU = 149597870.7;                  % km
params.solar_flux = 1357;                 % W/m^2
params.speed_of_light = 299792.458;       % km/s
params.area_mass_ratio = 0.01 * 1e-6;             % m^2/kg
params.epoch = 2456296.25;                % Initial epoch JD

%% Integrate using your model
fprintf('Integrating reference trajectory and STM...\n');

dynamics_func = @(t,X) twobodySunSRP(t,X,params);

% ODE options
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

[t_propagated, X_propagated] = ode45(dynamics_func, t_truth, X0_truth, options);

% X_propagated: [length(t_truth) x 56]

%% Error Analysis

N = length(t_truth);

pos_errors = zeros(N,3);
vel_errors = zeros(N,3);
stm_errors = zeros(N,1);

for i = 1:N
    % Your propagated state
    X_num = X_propagated(i,:)';
    
    % Truth state
    X_true = Xt_50(i,:)';
    
    % Position error
    pos_errors(i,:) = (X_num(1:3) - X_true(1:3))';
    
    % Velocity error
    vel_errors(i,:) = (X_num(4:6) - X_true(4:6))';
    
    % STM error
    STM_num = reshape(X_num(8:end),7,7);
    STM_true = reshape(X_true(8:end),7,7);
    
    stm_errors(i) = norm(STM_num - STM_true, 'fro') / norm(STM_true, 'fro'); % relative error
end

%% Compute RMS Errors

pos_rms = sqrt(mean(pos_errors.^2,1)); % km
vel_rms = sqrt(mean(vel_errors.^2,1)); % km/s
stm_rms = sqrt(mean(stm_errors.^2));   % unitless

%% Print Results

fprintf('\n--- Results ---\n');
fprintf('Position RMS Errors (km):\n');
fprintf('  X: %.6f  Y: %.6f  Z: %.6f\n', pos_rms(1), pos_rms(2), pos_rms(3));

fprintf('Velocity RMS Errors (km/s):\n');
fprintf('  VX: %.6f  VY: %.6f  VZ: %.6f\n', vel_rms(1), vel_rms(2), vel_rms(3));

fprintf('STM RMS Error (relative):\n');
fprintf('  STM RMS Error: %.6e\n', stm_rms);

