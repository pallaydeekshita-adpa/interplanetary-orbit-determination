function [state_out, Covariance, residuals, resid_pf, P_pf, X_pf] = filter_ekf(x_init, Pbar_init, obs, obsW, options)
%
% Advanced EKF implementation with multiple process noise strategies
% ---------------------------------------------------------------------
fprintf('Running Extended Kalman Filter with Multiple Process Noise Approaches...\n');

num_states   = length(x_init);
Phi0         = eye(num_states);
Phi0_reshape = reshape(Phi0,num_states*num_states,1);
x0           = x_init;
Pbar0        = Pbar_init;
xbar0        = zeros(num_states,1);
R            = inv(obsW);

abs_tol     = ones(1,num_states+num_states*num_states).*options.tol;
ode_options = odeset('RelTol',options.tol,'AbsTol',abs_tol);
CovarianceI = eye(num_states);

%% Process Noise Configuration
% Base process noise matrices
Q_nominal = zeros(num_states);  % Minimal process noise for well-modeled dynamics
Q_perturbation = diag([5e-3, 5e-3, 5e-3, 2.5e-15, 2.5e-15, 2.5e-15, 0]);  % Increased position noise
% Direction-specific noise (if maneuver was mainly in X direction)
Q_high = diag([0.5e-3, 0.05e-3, 0.05e-3, 2e-13, 5e-14, 5e-14, 0]);

% Define the anomaly time based on observed measurement gap
t_anomaly = 18674940;  % seconds
% Additional time after anomaly when dynamics might still be adapting
t_post_anomaly = t_anomaly + 600;

% Select the approach to use (only enable one at a time)
approach = 1;  % Choose from 1-6:
% 1 = Threshold-Based Process Noise
% 2 = Maneuver Window Detection
% 3 = Explicit Maneuver Modeling
% 4 = Adaptive Process Noise Based on Residuals
% 5 = Fading Memory Filter
% 6 = Schmidt-Kalman Consider Parameters

fprintf('Using approach %d: ', approach);
switch approach
    case 1
        fprintf('Threshold-Based Process Noise\n');
    case 2
        fprintf('Maneuver Window Detection\n');
    case 3
        fprintf('Explicit Maneuver Modeling\n');
    case 4
        fprintf('Adaptive Process Noise Based on Residuals\n');
    case 5
        fprintf('Fading Memory Filter\n');
    
end

%% APPROACH 2: Maneuver Window Detection
% Detect time gaps larger than expected
gap_threshold = 600;  % seconds
gap_indices = find(diff(obs.times) > gap_threshold);
gap_times = obs.times(gap_indices + 1);  % Time points after large gaps

% Build maneuver windows: apply enhanced process noise for steps after each gap
maneuver_windows = [];
window_size = 300;  % Number of steps in maneuver window
for i = 1:length(gap_times)
    t_start = gap_times(i);
    % Create a window of times after the gap
    for j = 0:(window_size-1)
        if i+j <= length(gap_indices)
            window_idx = gap_indices(i) + 1 + j;
            if window_idx <= length(obs.times)
                maneuver_windows = [maneuver_windows; obs.times(window_idx)];
            end
        end
    end
end

% Convert to days for reporting
gap_days = (gap_times - obs.times(1)) / 86400;
if ~isempty(gap_days)
    fprintf('Detected measurement gaps at the following times (days from epoch):\n');
    fprintf('%.2f\n', gap_days);
end

%% APPROACH 3: Explicit Maneuver Modeling
% Parameters for explicit maneuver at the anomaly time
% These can be tuned to minimize residuals
dV_maneuver = [0.0183; -0.0092; 0.0056];  % km/s
fprintf('Maneuver magnitude: %.4f km/s\n', norm(dV_maneuver));

%% APPROACH 4: Adaptive Process Noise Based on Residuals
% Initialize residual monitoring parameters
residual_window_size = 10;  % Number of measurements to consider
residual_threshold = 5.0;   % Multiple of expected measurement noise
residual_history = [];      % Keep track of recent residuals

%% APPROACH 5: Fading Memory Filter
% Set forgetting factor (between 0.95 and 0.99 typically)
alpha = 0.97;  % Smaller values = more aggressive forgetting of past


%% Main Filter Loop
status_chars = ['-' '\' '|' '/'];
fprintf('%c',status_chars(1));

for iterations = 1:1 % only one pass for EKF
    for j = 1:length(obs.times(:))
        fprintf('\b%c',status_chars(mod(j,4)+1));

        if j==1
            x_hat      = xbar0;
            Covariance = Pbar0;
            prev_time  = obs.times(1);
            PhiTotal   = Phi0;
            residuals  = zeros(length(obs.times),2);
            resid_pf   = zeros(length(obs.times),2);
            P_pf       = zeros(length(obs.times),num_states^2);
            X_pf       = zeros(length(obs.times),num_states);
            X_Int(1,:) = [x0; Phi0_reshape];
            Time_Int   = obs.times(1);
        else
            [Time_Int,X_Int] = ode45(@(t, x) options.integ_fcn(t, x, options.integ_args{:}), ...
                [prev_time obs.times(j)], [X_Ref; Phi0_reshape], ode_options);
            prev_time = obs.times(j);
        end

        Phi = reshape(X_Int(end,(num_states+1):end),num_states,num_states);
        PhiTotal = Phi*PhiTotal;

        % Observation model + sensitivity matrix
        state_now = X_Int(end, 1:num_states);
        t_now     = Time_Int(end);

        % Time Update
        xbar = Phi * x_hat;

        % Initialize process noise to nominal
        Q = Q_nominal;
        
        % Initialize Pbar (will be overridden by specific approaches if needed)
        Pbar = [];

        %% Apply selected approach for process noise/dynamics modification
        if approach == 1
            % APPROACH 1: Threshold-Based Process Noise
            if t_now > t_anomaly
                Q = Q_nominal + Q_high;
                if j == find(obs.times > t_anomaly, 1)
                    fprintf('\nApplying higher process noise after t = %.2f days\n', (t_anomaly - obs.times(1))/86400);
                end
            else
                Q = Q_nominal + Q_perturbation;
            end
            
            % Regular covariance propagation
            Pbar = Phi * Covariance * Phi' + Q;
            
        elseif approach == 2
            % APPROACH 2: Maneuver Window Detection
            if any(abs(t_now - maneuver_windows) < 1e-3)
                Q = Q_high;
                fprintf('\nApplying maneuver window process noise at t = %.2f days\n', (t_now - obs.times(1))/86400);
            else
                Q = Q_nominal + Q_perturbation;
            end
            
            % Regular covariance propagation
            Pbar = Phi * Covariance * Phi' + Q;
            
        elseif approach == 3
            % APPROACH 3: Explicit Maneuver Modeling
            % Apply nominal process noise
            Q = Q_nominal + Q_perturbation;

            % If we're crossing the anomaly time in this step, apply the maneuver
            if prev_time <= t_anomaly && t_now >= t_anomaly
                fprintf('\nApplying maneuver at t = %.2f days\n', (t_anomaly - obs.times(1))/86400);
                % Add the velocity change to the predicted state
                xbar(4:6) = xbar(4:6) + dV_maneuver;

                % Could also add uncertainty in the maneuver execution to the covariance
                % through process noise if desired
                Q(4:6, 4:6) = Q(4:6, 4:6) + diag([1e-10, 1e-10, 1e-10]);
            end
            
            % Regular covariance propagation
            Pbar = Phi * Covariance * Phi' + Q;
            
        elseif approach == 4
            % APPROACH 4: Adaptive Process Noise Based on Residuals
            % Add base level of process noise
            Q = Q_nominal + Q_perturbation;

            % After collecting enough measurements, check if residuals indicate model issues
            if j > 1
                % Store normalized residual for this measurement
                norm_resid = norm(residuals(j-1,:) ./ diag(R)');
                residual_history = [residual_history; norm_resid];

                % Keep only the most recent window
                if length(residual_history) > residual_window_size
                    residual_history = residual_history(end-residual_window_size+1:end);
                end

                % If we have enough history, check if residuals are growing
                if length(residual_history) >= 5 && ...
                        mean(residual_history(end-2:end)) > residual_threshold * mean(residual_history(1:3))
                    fprintf('\nDetected growing residuals at t = %.2f days, increasing process noise\n', ...
                        (t_now - obs.times(1))/86400);
                    Q = Q_high;
                end
            end
            
            % Regular covariance propagation
            Pbar = Phi * Covariance * Phi' + Q;
            
        elseif approach == 5
            % APPROACH 5: Fading Memory Filter
            % Base process noise
            Q = Q_nominal + Q_perturbation;

            % Apply fading memory factor to increase covariance, effectively
            % giving more weight to recent measurements
            % This is applied to Pbar directly rather than through Q
            Pbar = (1/alpha) * Phi * Covariance * Phi' + Q;

            % After anomaly, apply more aggressive forgetting
            if t_now > t_anomaly && t_now < t_post_anomaly
                fprintf('\nApplying fading memory at t = %.2f days\n', (t_now - obs.times(1))/86400);
                Pbar = (1/0.9) * Pbar;  % More aggressive forgetting near anomaly
            end
            
        else
            % Default for other approaches - use regular covariance propagation
            Pbar = Phi * Covariance * Phi' + Q;
        end

        % Ensure Pbar is computed
        if isempty(Pbar)
            Pbar = Phi * Covariance * Phi' + Q;
        end
        

        % Standard measurement model and Kalman gain calculation for other approaches
        Htilde = options.H_fcn(state_now, obs, t_now, options);
        KalScratch = Pbar * Htilde';
        KalGain = KalScratch / (Htilde * KalScratch + R);
        

        % Compute measurement and residuals
        C = options.obs_fcn(obs, state_now, t_now, options);
        OminusC = (obs.obs(j,:) - C)';
        residuals(j,:) = OminusC;

        % Measurement Update
        x_hat = xbar + KalGain * (OminusC - Htilde * xbar);

        % Covariance update with Joseph form for better numerical stability
        if options.joseph == 1
            Covariance = (CovarianceI - KalGain * Htilde) * Pbar * (CovarianceI - KalGain * Htilde)' ...
                + KalGain * R * KalGain';
        else
            Covariance = (CovarianceI - KalGain * Htilde) * Pbar;
        end

        % Post-fit residuals
        resid_pf(j,:) = OminusC - Htilde * x_hat;

        % Store state and covariance
        if ((j >= options.num_obs_switch) && (options.num_obs_switch > 0) && (j > 1) && (obs.station(j) == obs.station(j-1)))
            X_Ref = X_Int(end,1:num_states)' + x_hat;
            x_hat = zeros(num_states,1);
            X_pf(j,:) = X_Ref';
        else
            X_Ref = X_Int(end,1:num_states)';
            X_pf(j,:) = (X_Ref + x_hat)';
        end

        P_pf(j,:) = reshape(Covariance,1,num_states^2);
    end

    x_hat0 = PhiTotal \ x_hat;
    x0 = x0 + x_hat0;
    xbar0 = xbar0 - x_hat0;

    clear X_Int
end

residuals = [obs.times residuals];
state_out = X_pf(end,:)';

fprintf('\b');

% RMS Reporting
RMS_rho     = sqrt(mean(residuals(:,2).^2));
RMS_rhodot  = sqrt(mean(residuals(:,3).^2));
RMS_rho_pf  = sqrt(mean(resid_pf(:,1).^2));
RMS_rhodot_pf = sqrt(mean(resid_pf(:,2).^2));

fprintf('EKF Pre-fit residual RMS values:\n');
fprintf('Rho    = %.6f km\n', RMS_rho);
fprintf('RhoDot = %.6f km/s\n', RMS_rhodot);
fprintf('---------------------------------------------\n');
fprintf('EKF Post-fit residual RMS values:\n');
fprintf('Rho    = %.6f km\n', RMS_rho_pf);
fprintf('RhoDot = %.6f km/s\n', RMS_rhodot_pf);
fprintf('---------------------------------------------\n');

end