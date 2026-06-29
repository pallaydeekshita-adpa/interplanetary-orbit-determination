function plot_bplane_ellipses(Bdots, P_bplane_list, labels, target)
% PLOT_BPLANE_ELLIPSES
% Plots 3-sigma ellipses and estimate points on the B-plane
%
% Inputs:
%   Bdots         - Nx2 matrix of [BdotR, BdotT] rows
%   P_bplane_list - Cell array of 2x2 projected covariances (one per time)
%   labels        - Cell array of labels (e.g., {'50d','100d','150d','200d'})
%   target        - [BdotR_target, BdotT_target] (optional)

% Default colors
colors = {'r', 'b', 'g', 'm'};

% Start plot
figure; hold on; grid on;
xlabel('B \cdot \hat{R} [km]');
ylabel('B \cdot \hat{T} [km]');
title('B-plane 3σ Covariance Ellipses at Multiple Epochs');
axis equal;

% Plot ellipses and points
for i = 1:size(Bdots,1)
    mu = Bdots(i,:);
    P = P_bplane_list{i};
    error_ellipse(9 * P, mu, 'style', colors{i});  % 3σ ellipse
    plot(mu(1), mu(2), [colors{i} 'x'], 'MarkerSize',10, 'LineWidth',2);
    text(mu(1)+10, mu(2)+10, labels{i}, 'FontSize',9, 'Color',colors{i});
end

% Plot target point if provided
if nargin == 4
    plot(target(1), target(2), 'ko', 'MarkerFaceColor', 'k');
    text(target(1)+10, target(2), 'Target', 'FontSize',10, 'Color','k');
end

legend([labels, 'Target'], 'Location', 'best');
end
