function [BdotR_ideal, BdotT_ideal, BdotR, BdotT, r_bplane, v_bplane, STM_bplane, P_bplane] = ...
    propagate_to_bplane(r0, v0, CR0, covar, params)
% PROPAGATE_TO_BPLANE
% Propagates a state and covariance to B-plane (3*RSOI), computes BdotR/T and rotated covariance

% ---------------- Constants ----------------
mu_earth = params.mu_earth;
r_soi = 925000;                % Earth RSOI [km]
r3soi = 3 * r_soi;             % B-plane crossing distance
N_hat = [0; 0; 1];             % B-plane normal (Z direction)

% ---------------- Initial augmented state (7x7 STM) ----------------
STM0 = eye(7);
y0 = [r0; v0; CR0; reshape(STM0, 49, 1)];

% ---------------- Propagate to 3*RSOI ----------------
options = odeset('RelTol',1e-12,'AbsTol',1e-14, ...
                 'Events', @(t,y) reach_3rsoi(t,y,mu_earth,r3soi));

[~, ~, ~, Ye] = ode45(@(t,y) twobodySunSRP(t,y,params), [0 1e15], y0, options);

if isempty(Ye)
    error('Failed to reach 3 RSOI from given initial state.');
end

% ---------------- State and STM at 3*RSOI ----------------
r3 = Ye(1:3);
v3 = Ye(4:6);
STM_full_3rsoi = reshape(Ye(8:end), 7, 7);
STM_3rsoi = STM_full_3rsoi(1:6, 1:6);

% ---------------- Covariance Propagation ----------------
P_eci = covar(1:6,1:6);
P_rv_3rsoi = STM_3rsoi * P_eci * STM_3rsoi';

% ---------------- Hyperbolic Trajectory Geometry ----------------
vinf = norm(v3);
e_vec = ((vinf^2 - mu_earth/norm(r3)) * r3 - dot(r3,v3)*v3) / mu_earth;
e = norm(e_vec);
a = -mu_earth / vinf^2;
f = acosh(1 + (vinf^2 / mu_earth) * (a*(1 - e^2)) / (1 + e*cos(acos(dot(r3/norm(r3), e_vec/e)))));
LTOF = (mu_earth / vinf^3) * (sinh(f) - f);  % Linearized time-of-flight

% ---------------- BdotR_ideal from r3 and v3 ----------------
S_hat_3 = v3 / norm(v3);
T_hat_3 = cross(S_hat_3, N_hat); T_hat_3 = T_hat_3 / norm(T_hat_3);
R_hat_3 = cross(S_hat_3, T_hat_3);
B_vec_3 = r3 - dot(r3, S_hat_3) * S_hat_3;
BdotR_ideal = dot(B_vec_3, R_hat_3);
BdotT_ideal = dot(B_vec_3, T_hat_3);


% ---------------- Propagate to actual B-plane crossing ----------------
[~, Yb] = ode45(@(t,y) twobodySunSRP(t,y,params), [0 LTOF], Ye');

r_bplane = Yb(end,1:3)';
v_bplane = Yb(end,4:6)';
STM_full_bplane = reshape(Yb(end,8:end), 7, 7);
STM_bplane = STM_full_bplane(1:6, 1:6);

% ---------------- Final BdotR, BdotT at B-plane using same projection ----------------
S_hat_b = v_bplane / norm(v_bplane);
T_hat_b = cross(S_hat_b, N_hat); T_hat_b = T_hat_b / norm(T_hat_b);
R_hat_b = cross(S_hat_b, T_hat_b);
B_vec_b = r_bplane - dot(r_bplane, S_hat_b) * S_hat_b;
BdotR = dot(B_vec_b, R_hat_b);
BdotT = dot(B_vec_b, T_hat_b);
S_hat_b = S_hat_b';
T_hat_b = T_hat_b';
R_hat_b = R_hat_b';

% ---------------- Covariance Rotation into B-plane frame ----------------
P_r_r = P_rv_3rsoi(1:3, 1:3);             % position-position block
C_b = [R_hat_b; T_hat_b; S_hat_b];       % TRS frame as rows
P_r_TR = C_b * P_r_r * C_b';             % rotate into B-plane frame
P_bplane = P_r_TR(1:2, 1:2);             % extract [T,R] 2x2 block

end


%% Event Function: Reach 3 RSOI
function [value, isterminal, direction] = reach_3rsoi(~, y, mu, r_thresh)
    r = norm(y(1:3));
    value = r - r_thresh;
    isterminal = 1;
    direction = 0;
end
