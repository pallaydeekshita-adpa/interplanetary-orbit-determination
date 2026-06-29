%% Project2_exec.m - For ASEN 6080 Project 2 Part 2
% Clean version: Correct dynamics, measurements, covariance, units

clear; clc; close all;

fprintf("\n--- Starting Project 2 Execution ---\n");

%% 1. Load Truth Trajectory (First State)
load('Project2_Prob2_truth_traj_50days.mat'); % loads Tt_50 and Xt_50

X0_truth = Xt_50(1, :)'; % initial 56x1 state
X0_sc = X0_truth(1:6);    % first 6 entries are position and velocity
CR0 = X0_truth(7);         % 7th entry is CR

% Initial state (7x1)
x_init = [
    -274096790.0;  % X
    -92859240.0;   % Y
    -40199490.0;   % Z
    32.67;         % VX
    -8.94;         % VY
    -3.88;         % VZ
    1.2            % CR
]


%% 2. Setup Parameters
params.mu_earth = 398600.4415;           % km^3/s^2
params.mu_sun = 132712440017.987;         % km^3/s^2
params.AU = 149597870.7;                  % km
params.solar_flux = 1357;                 % W/m^2
params.speed_of_light = 299792.458;       % km/s
params.area_mass_ratio = 0.01 * 1e-6;             % m^2/kg
params.epoch = 2456296.25;                % JD at t = 0

RE = 6378.1363; % km
w_earth = 7.29211585275553e-5; % rad/s

%% 3. DSN Stations (ECEF at theta=0 deg)
stations_ecef = struct();

stations_ecef.Canberra = (RE + 0.691750) * [...
    cosd(-35.398333) * cosd(148.981944);
    cosd(-35.398333) * sind(148.981944);
    sind(-35.398333)];

stations_ecef.Madrid = (RE + 0.834539) * [...
    cosd(40.427222) * cosd(-355.749444);
    cosd(40.427222) * sind(-355.749444);
    sind(40.427222)];

stations_ecef.Goldstone = (RE + 1.07114904) * [...
    cosd(35.247164) * cosd(243.205);
    cosd(35.247164) * sind(243.205);
    sind(35.247164)];

%% 4. Load Observations
ObsData = importdata('Project2a_Obs.txt');
data = ObsData.data;

obs.times = [];
obs.station = [];
obs.obs = [];

% Loop through each row and extract valid measurements
for i = 1:size(data,1)
    t = data(i,1);

    % Canberra (station 34)
    if ~isnan(data(i,2)) && ~isnan(data(i,5))
        obs.times = [obs.times; t];
        obs.station = [obs.station; 34];
        obs.obs = [obs.obs; [data(i,2), data(i,5)]];
    end

    % Madrid (station 65)
    if ~isnan(data(i,3)) && ~isnan(data(i,6))
        obs.times = [obs.times; t];
        obs.station = [obs.station; 65];
        obs.obs = [obs.obs; [data(i,3), data(i,6)]];
    end

    % Goldstone (station 13)
    if ~isnan(data(i,4)) && ~isnan(data(i,7))
        obs.times = [obs.times; t];
        obs.station = [obs.station; 13];
        obs.obs = [obs.obs; [data(i,4), data(i,7)]];
    end
end

%% 5. Set Up Covariance and Measurement Noise
num_states = 7;

Pbar_init = diag([100^2, 100^2, 100^2, 0.1^2, 0.1^2, 0.1^2, 0.1^2]);

sigma_rho = 0.005;      % km (5 meters)
sigma_rhodot = 5 * 1e-7;  % km/s (0.5 mm/s)
obsW = [1/sigma_rho^2 0; 0 1/sigma_rhodot^2];

%% 6. Options for Filter
options = struct();
options.tol = 1e-12;
options.integ_fcn = @twobodySunSRP;
options.integ_args = {params};
options.obs_fcn = @epoch_range_rangerate;
options.H_fcn = @Htilde_6080;
options.extra_args.stations_ecef = stations_ecef;
options.extra_args.w_earth = w_earth;

options.start_time = 0;
options.end_time = max(obs.times);
options.inc_time = 60; % sec
options.conv_crit = 1e-7;
options.max_iterations = 5;

%% 7. Propagate Truth Trajectory for 200 Days
fprintf('Propagating truth trajectory for 200 days...\n');

X0_truth_propagate = [x_init; reshape(eye(7),49,1)];

dynamics_func = @(t,X) twobodySunSRP(t,X,params);

ode_opts = odeset('RelTol',1e-12,'AbsTol',1e-14);

time_vector = options.start_time:options.inc_time:options.end_time;
[t_truth, X_truth] = ode45(dynamics_func, time_vector, X0_truth_propagate, ode_opts);

%% 8. Run Batch Least Squares Filter
fprintf('Running Batch Filter...\n');

xbar0 = zeros(num_states,1);
[state_out, covar, resids, resids_pf, xhat] = filter_LS(x_init, xbar0, Pbar_init, obs, obsW, options);

%% 9. Plot Residuals
fprintf('Plotting residuals...\n');

sta1Ind = find(obs.station==34);
sta2Ind = find(obs.station==65);
sta3Ind = find(obs.station==13);

figure;
subplot(2,1,1);
plot(obs.times(sta1Ind)/3600, resids(sta1Ind,2),'r.',...
     obs.times(sta2Ind)/3600, resids(sta2Ind,2),'g.',...
     obs.times(sta3Ind)/3600, resids(sta3Ind,2),'b.');
ylabel('Range Pre-fit [km]');
title('Pre-fit Range Residuals');
legend('Canberra','Madrid','Goldstone');

grid on;
subplot(2,1,2);
plot(obs.times(sta1Ind)/3600, resids(sta1Ind,3),'r.',...
     obs.times(sta2Ind)/3600, resids(sta2Ind,3),'g.',...
     obs.times(sta3Ind)/3600, resids(sta3Ind,3),'b.');
ylabel('Range-Rate Pre-fit [km/s]');
xlabel('Time [hours]');
grid on;

figure;
subplot(2,1,1);
plot(obs.times(sta1Ind)/3600, resids_pf(sta1Ind,1),'r.',...
     obs.times(sta2Ind)/3600, resids_pf(sta2Ind,1),'g.',...
     obs.times(sta3Ind)/3600, resids_pf(sta3Ind,1),'b.');
ylabel('Range Post-fit [km]');
title('Post-fit Range Residuals');
legend('Canberra','Madrid','Goldstone');

grid on;
subplot(2,1,2);
plot(obs.times(sta1Ind)/3600, resids_pf(sta1Ind,2),'r.',...
     obs.times(sta2Ind)/3600, resids_pf(sta2Ind,2),'g.',...
     obs.times(sta3Ind)/3600, resids_pf(sta3Ind,2),'b.');
ylabel('Range-Rate Post-fit [km/s]');
xlabel('Time [hours]');
grid on;

%% 10. Propagate Best-Fit Trajectory
fprintf('Propagating best-fit trajectory...\n');

time_vector = options.start_time:options.inc_time:options.end_time;
X0_bestfit = [state_out; reshape(eye(7),49,1)];
[t_bestfit, X_bestfit] = ode45(@(t,x) twobodySunSRP(t,x,params), time_vector, X0_bestfit, ode_opts);

% Propagate covariance
Phist = zeros(length(X_bestfit),7*7);
for ii = 1:length(X_bestfit)
    Phi = reshape(X_bestfit(ii,8:end),7,7);
    Pmapped = Phi*covar*Phi';
    Phist(ii,:) = reshape(Pmapped,1,7*7);
end

% %% 11. Compute State Errors
% X_truth_interp = interp1(t_truth, X_truth(:,1:6), time_vector);
% state_errors = X_bestfit(:,1:6) - X_truth_interp;

%% 12. Propagate State Error and Covariance (for ±3σ bounds)
% Allocate arrays for state errors and covariance
state_errors = zeros(length(X_bestfit), 7);
Phist = zeros(length(X_bestfit), 49);

% Loop through and propagate both state error and covariance
for i = 1:length(X_bestfit)
    % Extract state transition matrix
    Phi = reshape(X_bestfit(i,8:end), 7, 7);
    
    % % Propagate state error using state transition matrix
    % state_errors(i,:) = (Phi * xhat)'; % only for batch

    
    % Propagate covariance
    Pmapped = Phi * covar * Phi';
    Phist(i,:) = reshape(Pmapped, 1, 49);
end

% Calculate 3-sigma bounds
sigma_bounds = zeros(length(time_vector), 7);
for i = 1:length(time_vector)
    P = reshape(Phist(i,:), 7, 7);
    sigma_bounds(i,:) = 3 * sqrt(diag(P));
end

%% 13. Plot Velocity Errors with 3-Sigma Envelopes
figure;
subplot(3,1,1);
plot(time_vector/3600, state_errors(:,4), 'b'); hold on;
plot(time_vector/3600, sigma_bounds(:,4), 'r--');
plot(time_vector/3600, -sigma_bounds(:,4), 'r--');
title('Vx Error with ±3σ'); ylabel('Vx [km/s]'); grid on;

subplot(3,1,2);
plot(time_vector/3600, state_errors(:,5), 'r'); hold on;
plot(time_vector/3600, sigma_bounds(:,5), 'k--');
plot(time_vector/3600, -sigma_bounds(:,5), 'k--');
ylabel('Vy [km/s]'); grid on;

subplot(3,1,3);
plot(time_vector/3600, state_errors(:,6), 'g'); hold on;
plot(time_vector/3600, sigma_bounds(:,6), 'm--');
plot(time_vector/3600, -sigma_bounds(:,6), 'm--');
ylabel('Vz [km/s]'); xlabel('Time [hours]'); grid on;

title('Velocity Errors in ECI Frame with ±3σ Bounds');

% Plot position error with ±3σ bounds
figure;
set(gcf, 'Color', 'w');

subplot(3,1,1);
plot(time_vector/3600, state_errors(:,1), 'b'); hold on;
plot(time_vector/3600, sigma_bounds(:,1), 'r--');
plot(time_vector/3600, -sigma_bounds(:,1), 'r--');
ylabel('X Error [km]');
title('X Position Error with ±3σ Bounds');
grid on;

subplot(3,1,2);
plot(time_vector/3600, state_errors(:,2), 'b'); hold on;
plot(time_vector/3600, sigma_bounds(:,2), 'r--');
plot(time_vector/3600, -sigma_bounds(:,2), 'r--');
ylabel('Y Error [km]');
title('Y Position Error with ±3σ Bounds');
grid on;

subplot(3,1,3);
plot(time_vector/3600, state_errors(:,3), 'b'); hold on;
plot(time_vector/3600, sigma_bounds(:,3), 'r--');
plot(time_vector/3600, -sigma_bounds(:,3), 'r--');
ylabel('Z Error [km]');
xlabel('Time [hours]');
title('Z Position Error with ±3σ Bounds');
grid on;


% %% 14. Plot XY Trajectory Comparison
% figure;
% plot(X_bestfit(:,1), X_bestfit(:,2), 'b');
% hold on;
% plot(X_truth_interp(:,1), X_truth_interp(:,2), 'r--');
% scatter(X_bestfit(1,1), X_bestfit(1,2), 100, 'g', 'filled');
% scatter(X_bestfit(end,1), X_bestfit(end,2), 100, 'k', 'filled');
% xlabel('X [km]');
% ylabel('Y [km]');
% legend('Best-Fit','Truth','Start','End');
% title('XY Plane Trajectory Comparison');
% grid on;
% axis equal;

%% 15. Post-fit Residual RMS Summary
fprintf('\nPost-fit residual RMS analysis:\n');
RMS_rho_pf = sqrt(mean(resids_pf(:,1).^2));
RMS_rhodot_pf = sqrt(mean(resids_pf(:,2).^2));
fprintf('RMS Range Residual = %.6f km\n', RMS_rho_pf);
fprintf('RMS Range-Rate Residual = %.6f km/s\n', RMS_rhodot_pf);

fprintf('--- Converged Project2_exec 2a finished! ---\n');

