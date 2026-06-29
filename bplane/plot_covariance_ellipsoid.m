function plot_covariance_ellipsoid(Covar, state_out, color, markerColor)
    % Plots the final 3D covariance ellipsoid for position using the covariance matrix.
    %
    % Covar: 18x18 covariance matrix of the estimated state
    % state_out: 18x1 final state estimate vector
    % color: Optional color for the ellipsoid surface

    hold on; % Allow multiple ellipsoids to be plotted on the same figure
    grid on;
    xlabel('X Position [m]'); ylabel('Y Position [m]'); zlabel('Z Position [m]');
    
    % Extract position covariance (assuming first 3 states are position)
    pos_cov = Covar(1:3, 1:3);
    
    % Compute eigenvalues and eigenvectors of the position covariance matrix
    [eigVec, eigVal] = eig(pos_cov);

    % Compute ellipsoid radii
    radii = sqrt(diag(eigVal));

    % Generate an ellipsoid at the origin
    [X, Y, Z] = ellipsoid(0, 0, 0, radii(1), radii(2), radii(3), 20);

    % Rotate the ellipsoid using the eigenvectors
    ellipsoid_points = [X(:), Y(:), Z(:)] * eigVec';

    % Offset the ellipsoid to the estimated position
    x_hat = state_out(1:3);
    X = reshape(ellipsoid_points(:, 1) + x_hat(1), size(X));
    Y = reshape(ellipsoid_points(:, 2) + x_hat(2), size(Y));
    Z = reshape(ellipsoid_points(:, 3) + x_hat(3), size(Z));

    % Plot the ellipsoid
    surf(X, Y, Z, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', color);  % Use the specified color

    % Plot the final state position
    plot3(x_hat(1), x_hat(2), x_hat(3), 'o', 'MarkerFaceColor', markerColor, 'MarkerEdgeColor', markerColor, 'MarkerSize', 8);


    axis equal; % Set equal scaling for all axes
    hold off; % Release the hold
end

