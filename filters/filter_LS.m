function [state_out, Covar, resids, resid_pf, xhat] = filter_LS(x_init, xbar, Pbar_init, obs, obsW, options)
%
% ---------------------------------------------------------------------
% Description: Rewritten Least Squares Filter (Project 2) with structure preserved

fprintf('Running Batch Filter...\n');

num_states    = length(x_init);
Phi0          = eye(num_states);
x0            = [x_init; reshape(Phi0,num_states*num_states,1)];
Pbar0         = Pbar_init;

time          = options.start_time:options.inc_time:options.end_time;
abs_tol       = ones(1,length(x0)).*options.tol;
converge_crit = options.conv_crit;
ode_options   = odeset( 'RelTol', options.tol, 'AbsTol', abs_tol);

status_chars = ['-' '\' '|' '/'];
fprintf('%c',status_chars(1));


j = 0;
xhat_mag = 10*converge_crit;

while xhat_mag > converge_crit
    j = j+1;

    [t,x] = ode45(@(t, x) options.integ_fcn(t, x, options.integ_args{:}), time, x0, ode_options);

    PbarR = chol(Pbar0);
    PbarR_inv = inv(PbarR);
    Pbar_inv = PbarR_inv*PbarR_inv';

    L = Pbar_inv;
    N = Pbar_inv*xbar;
    OminusC = zeros(length(obs.times),2);

    for i = 1:length(obs.times)
        index = find(obs.times(i)== t );
        if length(index) < 1
            continue
        end

        Phi = reshape(x(index,(num_states+1):end),num_states,num_states); 
        state_vector = x(index, 1:num_states);

        C = options.obs_fcn(obs, state_vector, t(index), options);
        OminusC(i,:) = obs.obs(i,:)-C;

        Htilde = options.H_fcn(state_vector, obs, t(index), options);
        H = Htilde*Phi;
        HtW = H'*obsW;

        L = L + HtW*H;
        N = N + HtW*OminusC(i,:)';
    end

    LR = chol(L);
    LR_inv = inv(LR);
    Covar = LR_inv*LR_inv';
    xhat = Covar*N;
    xhat_mag = norm(xhat);
    x0(1:num_states) = x0(1:num_states) + xhat;

    RMS_Rho = sqrt(sum(OminusC(:,1).^2)/length(OminusC(:,1)));
    RMS_RhoDot = sqrt(sum(OminusC(:,2).^2)/length(OminusC(:,2)));

    resid_pf = zeros(length(obs.times),2);
    for i = 1:length(obs.times)
        index = find(obs.times(i)== t );
        if length(index) < 1
            continue
        end
        Phi = reshape(x(index,(num_states+1):end),num_states,num_states);
        Htilde = options.H_fcn(state_vector, obs, t(index), options);
        H = Htilde*Phi;
        resid_pf(i,:) = OminusC(i,:) - (H*xhat)'; 
    end
    RMS_Rho_pf = sqrt(sum(resid_pf(:,1).^2)/length(resid_pf(:,1)));
    RMS_RhoDot_pf = sqrt(sum(resid_pf(:,2).^2)/length(resid_pf(:,2)));

    fprintf('--- Iteration %d info: ---\n', j );
    fprintf('Pre-fit residual RMS values\n');
    fprintf('Rho    = %g meters\n', RMS_Rho );
    fprintf('RhoDot = %g meters\n', RMS_RhoDot );
    fprintf('Post-fit residual RMS values\n');
    fprintf('Rho    = %g meters\n', RMS_Rho_pf );
    fprintf('RhoDot = %g meters\n', RMS_RhoDot_pf );
    fprintf('x_hat magnitude:  %g\n', xhat_mag );
    fprintf('\n---------------------------------------------\n\n');

    resids = [obs.times OminusC];
    sta1Ind = find(obs.station==34);
    sta2Ind = find(obs.station==65);
    sta3Ind = find(obs.station==13);

    figure();
    set(gcf,'Color','w');
    subplot(2,1,1);
    plot(resids(sta1Ind,1)./(3600),resids(sta1Ind,2),'.',resids(sta2Ind,1)./(3600),resids(sta2Ind,2),'.',resids(sta3Ind,1)./(3600),resids(sta3Ind,2),'.');
    title(['Pre-fit Residuals: Iteration ' num2str(j)])
    ylabel('Range [m]');
    subplot(2,1,2);
    plot(resids(sta1Ind,1)./(3600),resids(sta1Ind,3),'.',resids(sta2Ind,1)./(3600),resids(sta2Ind,3),'.',resids(sta3Ind,1)./(3600),resids(sta3Ind,3),'.');
    ylabel('Range-Rate [m/s]')
    xlabel('Time (hours)')

    figure();
    set(gcf,'Color','w');
    subplot(2,1,1);
    plot(resids(sta1Ind,1)./(3600),resid_pf(sta1Ind,1),'.',resids(sta2Ind,1)./(3600),resid_pf(sta2Ind,1),'.',resids(sta3Ind,1)./(3600),resid_pf(sta3Ind,1),'.');
    title(['Post-fit Residuals: Iteration ' num2str(j)])
    ylabel('Range [m]');
    subplot(2,1,2);
    plot(resids(sta1Ind,1)./(3600),resid_pf(sta1Ind,2),'.',resids(sta2Ind,1)./(3600),resid_pf(sta2Ind,2),'.',resids(sta3Ind,1)./(3600),resid_pf(sta3Ind,2),'.');
    ylabel('Range-Rate [m/s]')
    xlabel('Time (hours)')

    xbar = xbar - xhat;
    if j >= options.max_iterations
        break;
    end
end

state_out = x0(1:num_states);
resids = [obs.times OminusC];