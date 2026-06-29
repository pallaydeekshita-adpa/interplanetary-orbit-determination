function [state_out, Covariance, residuals, resid_pf, P_pf, X_pf] = filter_ck(x_init, Pbar_init, obs, obsW, options)

fprintf('Running Classical Kalman Filter...\n');

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

status_chars = ['-' '\' '|' '/'];
fprintf('%c',status_chars(1));

for iterations = 1:options.num_iters
    for j = 1:length(obs.times(:))
        fprintf('\b%c',status_chars(mod(j,4)+1));

        if j == 1
            x_hat = xbar0;
            Covariance = Pbar0;
            prev_time = 0;
            PhiTotal = Phi0;
            residuals = zeros(length(obs.times), 2);
            resid_pf = zeros(length(obs.times), 2);
            P_pf = zeros(length(obs.times), num_states^2);
            X_pf = zeros(length(obs.times), num_states);
            X_Int(1,:) = [x0; Phi0_reshape];
            Time_Int = 0;
        else
            [Time_Int, X_Int] = ode45(@(t, x) options.integ_fcn(t, x, options.integ_args{:}), ...
                [prev_time obs.times(j)], [X_Ref; Phi0_reshape], ode_options);
            prev_time = obs.times(j);
        end

        Phi = reshape(X_Int(end,(num_states+1):end),num_states,num_states);
        PhiTotal = Phi * PhiTotal;

        % Time Update
        xbar = Phi * x_hat;
        Pbar = Phi * Covariance * Phi';

        % Observation and Htilde using consistent structure
        state_now = X_Int(end,1:num_states);
        t_now = Time_Int(end);

        C = options.obs_fcn(obs, state_now, t_now, options);
        Htilde = options.H_fcn(state_now, obs, t_now, options);
        OminusC = (obs.obs(j,:) - C)';
        residuals(j,:) = OminusC;

        KalGain = Pbar * Htilde' / (Htilde * Pbar * Htilde' + R);

        % Measurement Update
        x_hat = xbar + KalGain * (OminusC - Htilde * xbar);

        if options.joseph == 1
            Covariance = (CovarianceI - KalGain * Htilde) * Pbar * (CovarianceI - KalGain * Htilde)' + KalGain * R * KalGain';
        else
            Covariance = (CovarianceI - KalGain * Htilde) * Pbar;
        end

        X_Ref = state_now';
        P_pf(j,:) = reshape(Covariance,1,num_states^2);
        resid_pf(j,:) = OminusC - Htilde * x_hat;
        X_pf(j,:) = (X_Ref + x_hat)';
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
RMS_rho     = sqrt(mean(residuals(:,1).^2));
RMS_rhodot  = sqrt(mean(residuals(:,2).^2));
RMS_rho_pf  = sqrt(mean(resid_pf(:,1).^2));
RMS_rhodot_pf = sqrt(mean(resid_pf(:,2).^2));

fprintf('CKF Pre-fit residual RMS values:\n');
fprintf('Rho    = %.6f km\n', RMS_rho);
fprintf('RhoDot = %.6f km/s\n', RMS_rhodot);
fprintf('---------------------------------------------\n');
fprintf('CKF Post-fit residual RMS values:\n');
fprintf('Rho    = %.6f km\n', RMS_rho_pf);
fprintf('RhoDot = %.6f km/s\n', RMS_rhodot_pf);
fprintf('---------------------------------------------\n');

end


