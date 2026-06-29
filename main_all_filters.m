tic

%% Project2_exec.m - Full Filter Post-Processing like HW2_exec
% CKF, EKF, Batch with complete post-analysis and visualization

clear; clc; close all;

fprintf("\n--- Starting Project 2 Execution (with full post-analysis) ---\n");

%% 0. Select Filter: 1 = Batch, 2 = CKF, 3 = EKF
useFilter = 3;

% %% 1. Load Truth Trajectory (First State)
% load('Project2_Prob2_truth_traj_50days.mat'); % loads Tt_50 and Xt_50
% 
% X0_truth = Xt_50(1, :)';
% X0_sc = X0_truth(1:6);
% CR0 = X0_truth(7);
 
% %%% disabled for part 3

%% A Priori State

% % Part 2 
% x_init = [
%     -274096790.0;  % X
%     -92859240.0;   % Y
%     -40199490.0;   % Z
%     32.67;         % VX
%     -8.94;         % VY
%     -3.88;         % VZ
%     1.2            % CR
% ];

% % Part 3
x_init = [
    -274096770.76544;
    -92859266.4499061;
    -40199493.6677441;
    32.6704564599943;
    -8.93838913761049;
    -3.87881914050316;
    1
];

%% 2. Setup Parameters
params.mu_earth = 398600.4415;
params.mu_sun = 132712440017.987;
params.AU = 149597870.7;
params.solar_flux = 1357;
params.speed_of_light = 299792.458;
params.area_mass_ratio = 0.01 * 1e-6;
params.epoch = 2456296.25;

RE = 6378.1363;
w_earth = 7.29211585275553e-5;

%% 3. DSN Stations (ECEF at theta=0 deg)
stations_ecef = struct();

stations_ecef.Canberra = (RE + 0.691750) * [...
    cosd(-35.398333) * cosd(148.981944);
    cosd(-35.398333) * sind(148.981944);
    sind(-35.398333)];

stations_ecef.Madrid = (RE + 0.834539) * [...
    cosd(40.427222) * cosd(355.749444);
    cosd(40.427222) * sind(355.749444);
    sind(40.427222)];

stations_ecef.Goldstone = (RE + 1.07114904) * [...
    cosd(35.247164) * cosd(243.205);
    cosd(35.247164) * sind(243.205);
    sind(35.247164)];

%% 4. Load Observations
ObsData = importdata('Project2b_Obs.txt');
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

%% 5. Covariance and Noise
num_states = 7;
Pbar_init = diag([100^2, 100^2, 100^2, 0.1^2, 0.1^2, 0.1^2, 0.1^2]);
sigma_rho = 0.005;
sigma_rhodot = 5e-7;
obsW = [1/sigma_rho^2 0; 0 1/sigma_rhodot^2];

%% 6. Filter Options
options = struct();
options.tol = 1e-12;
options.integ_fcn = @twobodySunSRP;
options.integ_args = {params};
options.obs_fcn = @epoch_range_rangerate;
options.H_fcn = @Htilde_6080;
options.extra_args.stations_ecef = stations_ecef;
options.extra_args.w_earth = w_earth;
% Batch Options
options.start_time     = 0;
options.end_time       = max(obs.times);
options.inc_time       = 60; % The filter updates state estimates at every 20-second interval.
options.conv_crit      = 0.001; % 0.3; %The filter stops iterating when the correction (update to the state) is smaller than 0.001.
options.max_iterations = 3; %The filter will perform a maximum of 10 iterations to refine the solution

% CKF Options
options.num_iters     = 3; % do this many iterations of KF
options.joseph        = 1; % Enables the Joseph stabilized update formula, which improves numerical stability when updating the covariance matrix.

% EKF Options
options.num_obs_switch = 100;
% CKF initializes state before EKF takes over
% run CKF for 100 observations before transitioning to EKF
% If the initial guess of the state is poor, the EKF Jacobians may be inaccurate.
% CKF smoothly refines the estimate, so when EKF starts, the system is already close to the true orbit.


xbar0 = zeros(num_states,1);

%% 7. Run Filter
switch useFilter
    case 1
        [state_out, covar, resids, resids_pf, xhat] = filter_LS(x_init, xbar0, Pbar_init, obs, obsW, options);
    case 2
        [state_out, covar, resids, resids_pf, P_pf, X_pf] = filter_ck(x_init, Pbar_init, obs, obsW, options);
    case 3
        [state_out, covar, resids, resids_pf, P_pf, X_pf] = filter_ekf(x_init, Pbar_init, obs, obsW, options);
end



%% 8. Propagate and Post-Process
ode_options = odeset('RelTol', options.tol, 'AbsTol', 1e-14);

switch useFilter
    case 1
        [~, X_bestfit] = ode45(@(t,x) twobodySunSRP(t,x,params), options.start_time:options.inc_time:options.end_time, [state_out; reshape(eye(num_states),num_states^2,1)], ode_options);
        Phist = zeros(length(X_bestfit), num_states^2);
        for ii = 1:length(X_bestfit)
            Phi = reshape(X_bestfit(ii,(num_states+1):end),num_states,num_states);
            Pmapped = Phi * covar * Phi';
            Phist(ii,:) = reshape(Pmapped,1,num_states^2);
        end
    otherwise
        [~, Xprop] = ode45(@(t,x) twobodySunSRP(t,x,params), fliplr(options.start_time:options.inc_time:options.end_time), [state_out; reshape(eye(num_states),num_states^2,1)], ode_options);
        X_bestfit = flipud(Xprop);
        Phist = zeros(length(X_bestfit), num_states^2);
        for ii = 1:length(X_bestfit)
            Phi = reshape(Xprop(ii,(num_states+1):end),num_states,num_states);
            Pmapped = Phi * covar * Phi';
            Phist(ii,:) = reshape(Pmapped,1,num_states^2);
        end
        Phist = flipud(Phist);
end

%% 7. Propagate Truth Trajectory for 200 Days
fprintf('Propagating truth trajectory for 200 days...\n');

X0_truth_propagate = [x_init; reshape(eye(7),49,1)];

dynamics_func = @(t,X) twobodySunSRP(t,X,params);

ode_opts = odeset('RelTol',1e-12,'AbsTol',1e-14);

time_vector = options.start_time:options.inc_time:options.end_time;
[t_truth, X_truth] = ode45(dynamics_func, time_vector, X0_truth_propagate, ode_opts);

%% Run analysis and plots
% Can adjust starting index to only compute steady state RMS
RMS_zero_ind = 1;

RMS_rho_pf = sqrt(sum(resids_pf(RMS_zero_ind:end,1).*resids_pf(RMS_zero_ind:end,1))/length(resids_pf(RMS_zero_ind:end,1)));
RMS_rhodot_pf = sqrt(sum(resids_pf(RMS_zero_ind:end,2).*resids_pf(RMS_zero_ind:end,2))/length(resids_pf(RMS_zero_ind:end,2)));

% Computes the Root Mean Square (RMS) error of post-fit residuals for range (ρ) and range-rate (ρ̇).
% The post-fit residuals are the differences between the observed and computed measurements after filtering.
% If RMS is small, the filter fits the data well.

fprintf('\nPost-fit residuals\n');
fprintf('Rho    = %g km\n', RMS_rho_pf );
fprintf('RhoDot = %g km/s\n', RMS_rhodot_pf );

fprintf('--------------------------------------------------------\n');


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

%% 10. Additional Final State Plots (Position/Velocity vs Time)
final_state = X_bestfit(end, 1:7)';
final_covariance = reshape(Phist(end, :), 7, 7);

% Optional: plot_covariance_ellipsoid(final_covariance, final_state, 'b', 'r');

time_vector = options.start_time:options.inc_time:options.end_time;
if length(time_vector) ~= size(X_bestfit, 1)
    error('Mismatch between time vector and estimated states.');
end

x_position = X_bestfit(:,1);
y_position = X_bestfit(:,2);
z_position = X_bestfit(:,3);

x_velocity = X_bestfit(:,4);
y_velocity = X_bestfit(:,5);
z_velocity = X_bestfit(:,6);

figure;
subplot(3,1,1);
plot(time_vector, x_position, 'b', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('X Position [km]'); title('Estimated X Position vs Time'); grid on;
subplot(3,1,2);
plot(time_vector, y_position, 'r', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Y Position [km]'); title('Estimated Y Position vs Time'); grid on;
subplot(3,1,3);
plot(time_vector, z_position, 'g', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Z Position [km]'); title('Estimated Z Position vs Time'); grid on;

figure;
subplot(3,1,1);
plot(time_vector, x_velocity, 'b', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('X Velocity [km/s]'); title('Estimated X Velocity vs Time'); grid on;
subplot(3,1,2);
plot(time_vector, y_velocity, 'r', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Y Velocity [km/s]'); title('Estimated Y Velocity vs Time'); grid on;
subplot(3,1,3);
plot(time_vector, z_velocity, 'g', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Z Velocity [km/s]'); title('Estimated Z Velocity vs Time'); grid on;

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

%% 14. Plot XY Trajectory Comparison
figure;
plot(X_bestfit(:,1), X_bestfit(:,2), 'b');
hold on;
plot(X_truth(:,1), X_truth(:,2), 'r--');
scatter(X_bestfit(1,1), X_bestfit(1,2), 100, 'g', 'filled');
scatter(X_bestfit(end,1), X_bestfit(end,2), 100, 'k', 'filled');
xlabel('X [km]');
ylabel('Y [km]');
legend('Best-Fit','Truth','Start','End');
title('XY Plane Trajectory Comparison');
grid on;
axis equal;

%% 14. Filter-based Estimate Error Visualization
if useFilter > 1
    state_errors_est = X_pf(:,1:6) - interp1(t_truth, X_truth(:,1:6), obs.times);

    figure(); set(gcf, 'Color', 'w');
    labels = {'X-ECI [km]', 'Y-ECI [km]', 'Z-ECI [km]'};
    for ii = 1:3
        subplot(3,1,ii);
        plot(obs.times / 3600, state_errors_est(:,ii), '.');
        ylabel(labels{ii});
        if ii == 1, title('KF xhat Position Errors (ECI)'); end
    end
    xlabel('Time [hours]');

    figure(); set(gcf, 'Color', 'w');
    labels = {'Vx-ECI [km/s]', 'Vy-ECI [km/s]', 'Vz-ECI [km/s]'};
    for ii = 4:6
        subplot(3,1,ii-3);
        plot(obs.times / 3600, state_errors_est(:,ii), '.');
        ylabel(labels{ii-3});
        if ii == 4, title('KF xhat Velocity Errors (ECI)'); end
    end
    xlabel('Time [hours]');

    figure(); set(gcf, 'Color', 'w');
    for ii = 1:6
        valid_idx = P_pf(:,(ii-1)*6 + ii) >= 0;
        if any(valid_idx)
            semilogy(obs.times(valid_idx)/3600, sqrt(P_pf(valid_idx,(ii-1)*6 + ii)), '-o');
            hold on;
        end
    end
    ylabel('\sigma [km, km/s]');
    xlabel('Time [hours]');
    legend('x','y','z','vx','vy','vz');
    title('KF State Uncertainty Over Time');

    figure();
    set(gcf, 'Color', 'w');
    diag_inds = (0:5)*6 + 1;
    trace_P = sum(P_pf(:,diag_inds), 2);
    semilogy(obs.times/3600, trace_P, '-o');
    ylabel('Trace(P)');
    xlabel('Time [hours]');
    title('Trace of State Covariance');

    trace_P_position = sum(P_pf(:,1:3), 2);
    trace_P_velocity = sum(P_pf(:,4:6), 2);

    figure('Color','w');
    semilogy(obs.times/3600, trace_P_position, '-o', 'DisplayName', 'Position'); hold on;
    semilogy(obs.times/3600, trace_P_velocity, '-o', 'DisplayName', 'Velocity');
    legend('show');
    xlabel('Time [hours]');
    ylabel('Trace (log scale)');
    title('Position vs Velocity Trace Comparison');
    grid on;
end

toc