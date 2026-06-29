function [X_new, STM] = propagate_state_and_stm(X0, t0, tf, params)
    % This function propagates the state and STM from t0 to tf
    % Inputs:
    %   X0 - Initial state vector [r; v; CR]
    %   t0 - Initial time (seconds from epoch)
    %   tf - Final time (seconds from epoch)
    %   params - Structure containing parameters
    % Outputs:
    %   X_new - Propagated state
    %   STM - State transition matrix from t0 to tf
    
    % State dimension
    n = length(X0);
    
    % Initialize STM
    STM0 = eye(n);
    
    % Augmented state for integration [X; STM(:)]
    X_aug0 = [X0; STM0(:)];
    
    opts = odeset('RelTol',1e-12, 'AbsTol',1e-14);

    % Integrate from t0 to tf
   [t, X_aug] = ode45(@(t, X_aug) augmented_dynamics(t, X_aug, params), [t0 tf], X_aug0, opts);

    
    % Extract final state and STM
    X_new = X_aug(end, 1:n)';
    STM = reshape(X_aug(end, n+1:end), n, n);
end