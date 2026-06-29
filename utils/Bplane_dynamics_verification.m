clc;clear;close all;

% --- Load truth trajectory ---
load('Project2_Prob2_truth_traj_50days.mat'); % This loads Tt_50, Xt_50

% --- Constants ---
mu_earth = 398600.4415;                % [km^3/s^2]
r_soi = 925000;                        % [km]
r3soi = 3 * r_soi;    % [km]

params.mu_earth = 398600.4415;           % km^3/s^2
params.mu_sun = 132712440017.987;         % km^3/s^2
params.AU = 149597870.7;                  % km
params.solar_flux = 1357;                 % W/m^2
params.speed_of_light = 299792.458;       % km/s
params.area_mass_ratio = 0.01 * 1e-6;             % km^2/kg
params.epoch = 2456296.25;                % JD at t = 0


N_hat = [0;0;1];                       % Corrected North pole vector

% --- Start from first state (t = 0) ---
r0 = Xt_50(1,1:3)';
v0 = Xt_50(1,4:6)';
CR0 = 1.0;
STM0 = eye(7);
y0 = [r0; v0; CR0; reshape(STM0,49,1)];

% --- Propagate to 3 RSOI ---
options = odeset('RelTol',1e-12,'AbsTol',1e-14, ...
    'Events', @(t,y) reach_3rsoi(t,y,mu_earth,r3soi));

[~,Y,~,Ye] = ode45(@(t,y) twobodySunSRP(t,y,params), [0 1e18], y0, options);

% Extract state at 3 RSOI
r3 = Ye(1:3);
v3 = Ye(4:6);
STM_full = reshape(Ye(8:end),7,7);
STM_3rsoi = STM_full(1:6,1:6);

fprintf('Reached 3 RSOI at r = %.3f km\n', norm(r3));

% --- Orbital elements at 3 RSOI ---
vinf = norm(v3);
e_vec = ((vinf^2 - mu_earth/norm(r3)) * r3 - dot(r3,v3)*v3) / mu_earth;
e = norm(e_vec);
a = -mu_earth / (vinf^2);              % Hyperbolic a < 0
b = abs(a) * sqrt(e^2 - 1);
h_vec = cross(r3,v3);
W_hat = h_vec / norm(h_vec);
P_hat = e_vec / norm(e_vec);

% --- True anomaly and LTOF ---
cos_nu = dot(r3/norm(r3), P_hat);
nu = acos(cos_nu);
f = acosh(1 + (vinf^2 / mu_earth) * (a*(1-e^2)) / (1 + e*cos(nu)));
LTOF = (mu_earth / vinf^3) * (sinh(f) - f);

% B-vector projection using r_inf method
S_hat_3 = v3 / norm(v3);            % v-infinity direction
r_inf = r3;                         % position at 3 RSOI (close to infinity)
B_vec_3 = r_inf - dot(r_inf, S_hat_3) * S_hat_3;  % Projection onto B-plane

% Build B-plane frame
T_hat_3 = cross(S_hat_3, N_hat);
T_hat_3 = T_hat_3 / norm(T_hat_3);
R_hat_3 = cross(S_hat_3, T_hat_3);

% Project B-vector onto TR frame
BdotR_ideal = dot(B_vec_3, R_hat_3);
BdotT_ideal = dot(B_vec_3, T_hat_3);

% --- Propagate from 3 RSOI to B-plane crossing ---
[~,Yb] = ode45(@(t,y) twobodySunSRP(t,y,params), [0 LTOF], Ye');

r_bplane = Yb(end,1:3)';
v_bplane = Yb(end,4:6)';
STM_bplane = reshape(Yb(end,8:end),7,7);
STM_bplane_6x6 = STM_bplane(1:6,1:6);

% --- Actual B-plane projection at crossing ---
[BdotR, BdotT, S_hat, T_hat, R_hat, C_Bplane] = compute_bplane(r_bplane, v_bplane, mu_earth);

fprintf('\n----- Actual B-Plane Crossing Coordinates -----\n');
fprintf('BdotR (real)  = %.3f km\n', BdotR);
fprintf('BdotT (real)  = %.3f km\n', BdotT);

fprintf('\n----- B-Plane Results from Truth -----\n');
fprintf('BdotR (truth) = %.3f km\n', BdotR);
fprintf('BdotT (truth) = %.3f km\n', BdotT);
fprintf('\n----- B-Plane Geometry at 3 RSOI (Ideal Hyperbola) -----\n');
fprintf('BdotR (ideal) = %.3f km\n', BdotR_ideal);
fprintf('BdotT (ideal) = %.3f km\n', BdotT_ideal);

fprintf('Expected Target BdotR = 14970.824 km\n');
fprintf('Expected Target BdotT = 9796.737 km\n');
fprintf('--------------------------------------\n');



%% Helper
function [value,isterminal,direction] = reach_3rsoi(~, y, mu, r_thresh)
    r = norm(y(1:3));
    value = r - r_thresh;
    isterminal = 1;
    direction = 0;
end
