function [X, Y] = error_ellipse(C, mu, nsig)
% ERROR_ELLIPSE generates the 2D ellipse (X, Y) representing nsig-sigma error boundary
%
% Inputs:
%   C    - 2x2 positive definite covariance matrix
%   mu   - 1x2 or 2x1 vector (center of the ellipse)
%   nsig - Number of standard deviations (e.g., 1, 2, or 3 for 1σ–3σ bounds)
%
% Outputs:
%   X, Y - Coordinates of the ellipse boundary

    % Generate unit circle
    theta = linspace(0, 2*pi, 100);
    circle = [cos(theta); sin(theta)];

    % Eigen-decomposition
    [V, D] = eig(C);
    axes_lengths = nsig * sqrt(diag(D));  % scale eigenvalues for nsig bounds

    % Shape of ellipse
    ellipse = V * diag(axes_lengths) * circle;

    % Offset by mean (mu)
    X = ellipse(1,:) + mu(1);
    Y = ellipse(2,:) + mu(2);
end
