%% Project2_exec_Bplane.m - Unified B-Plane Analysis for all Filters

clear; clc; close all;

fprintf("\n--- Starting Project 2 B-Plane Covariance Analysis ---\n");

%% Select Filter: 1 = Batch, 2 = CKF, 3 = EKF
useFilter = 1;

%% 1. Load Truth Trajectory
load('Project2_Prob2_truth_traj_50days.mat');
X0_truth = Xt_50(1,:)';
X0_sc = X0_truth(1:6);
CR0 = X0_truth(7);
 x_init = [
    -274096790.0;  % X
    -92859240.0;   % Y
    -40199490.0;   % Z
    32.67;         % VX
    -8.94;         % VY
    -3.88;         % VZ
    1.2            % CR
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
options.start_time = 0;
options.end_time = max(obs.times);
options.inc_time = 60;
options.conv_crit = 1e-7;
options.max_iterations = 5;
options.num_iters = 3;
options.joseph = 1;
options.num_obs_switch = 100;

xbar0 = zeros(num_states,1);

%% 7 & 8. Run Filter Separately for Each Cutoff Time and Propagate to B-plane

bplane_days = [50, 100, 150, 200];  
BdotR_list = zeros(length(bplane_days),1);
BdotT_list = zeros(length(bplane_days),1);
P_bplane_list = cell(length(bplane_days),1);
labels = {'50d', '100d', '150d', '200d'};

fprintf('\n--- Evaluating Filter and B-plane for Snapshot Epochs ---\n');

for i = 1:length(bplane_days)
    % Cutoff time
    t_end = bplane_days(i) * 86400;
    options.end_time = t_end;

    % Subset observations up to cutoff
    valid_idx = find(obs.times <= t_end);
    obs_cut.times = obs.times(valid_idx);
    obs_cut.station = obs.station(valid_idx);
    obs_cut.obs = obs.obs(valid_idx, :);

    switch useFilter
        case 1  % Batch
            [state_out, covar, ~, ~, ~] = filter_LS(x_init, xbar0, Pbar_init, obs_cut, obsW, options);
            r = state_out(1:3);
            v = state_out(4:6);
            CR = state_out(7);

        otherwise  % CKF or EKF: integrate backward to get initial estimate
            if useFilter == 2
                [state_out, covar, ~, ~, P_pf, X_pf] = filter_ck(x_init, Pbar_init, obs_cut, obsW, options);
            elseif useFilter == 3
                [state_out, covar, ~, ~, P_pf, X_pf] = filter_ekf(x_init, Pbar_init, obs_cut, obsW, options);
            end

            r = X_pf(1, 1:3)';
            v = X_pf(1, 4:6)';
            CR = X_pf(7);
            covar = reshape(P_pf(1, :),7,7);
    end

    % Propagate to B-plane from r, v, covar
    [BdotR_ideal, BdotT_ideal, BdotR, BdotT, r_bplane, v_bplane, STM_bplane, P_bplane] = ...
        propagate_to_bplane(r, v, CR, covar, params);

    % Store results
    BdotR_list(i) = BdotR;
    BdotT_list(i) = BdotT;
    P_bplane_list{i} = P_bplane;

    fprintf('→ %s: BdotR = %.2f km, BdotT = %.2f km\n', labels{i}, BdotR, BdotT);

end


%% Plot 3σ B-plane Ellipses with BdotT on X, BdotR on Y

figure;
set(gcf, 'Color', 'w');
hold on; 
axis equal; 
grid on;

labels = {'50 days', '100 days', '150 days', '200 days'};
colors = {'#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd'};  % blue, orange, green, purple

for i = 1:length(bplane_days)
    % Rotate covariance and center [BdotR, BdotT]
    [ell_x, ell_y] = error_ellipse(P_bplane_list{i}, [BdotR_list(i), BdotT_list(i)], 3);

    % Flip to plot T on x-axis, R on y-axis
    plot(ell_y, ell_x, '--', 'LineWidth', 2, 'Color', colors{i});
    plot(BdotT_list(i), BdotR_list(i), 'o', 'MarkerSize', 6, ...
         'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k');

    % Add text label
    text(BdotT_list(i) + 5, BdotR_list(i) + 5, labels{i}, ...
         'Color', colors{i}, 'FontSize', 10, 'FontWeight', 'bold');
end

% Plot truth marker
target = [14970.824, 9796.737];  % [BdotR, BdotT]
plot(target(2), target(1), 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
text(target(2)+5, target(1)+5, 'Truth', 'FontSize', 10, 'FontWeight', 'bold');

xlabel('BdotT [km]', 'FontSize', 12);
ylabel('BdotR [km]', 'FontSize', 12);
title('\bf\it 3σ B-Plane Covariance Ellipses and Estimates', 'FontSize', 14);

legend({'50d', '100d', '150d', '200d', 'Truth'}, 'Location', 'southeast');



