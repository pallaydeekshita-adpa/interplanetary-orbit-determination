function [BdotR, BdotT, S_hat, T_hat, R_hat, C_Bplane] = compute_bplane(r_vec, v_vec, mu_earth)
% COMPUTE_BPLANE
% Computes the B-plane coordinates BdotR and BdotT given a state
%
% Inputs:
%   r_vec      - Position vector [km]
%   v_vec      - Velocity vector [km/s]
%   mu_earth   - Gravitational parameter of Earth [km^3/s^2]
%
% Outputs:
%   BdotR      - B-plane R coordinate [km]
%   BdotT      - B-plane T coordinate [km]
%   S_hat      - Unit vector of incoming v-infinity
%   T_hat      - T vector on B-plane
%   R_hat      - R vector on B-plane
%   C_Bplane   - DCM from ECI frame to B-plane frame

% North Pole direction (corrected)
N_hat = [0; 0; 1];

% Incoming hyperbolic asymptote direction
S_hat = v_vec / norm(v_vec);

% B-plane basis
T_hat = cross(S_hat, N_hat); T_hat = T_hat / norm(T_hat);
R_hat = cross(S_hat, T_hat);

% Build rotation matrix ECI -> B-plane
C_Bplane = [T_hat R_hat S_hat];

% Orbit properties
vinf = norm(v_vec);

e_vec = ((vinf^2 - mu_earth/norm(r_vec)) * r_vec - dot(r_vec,v_vec)*v_vec) / mu_earth;
e = norm(e_vec);
a = -mu_earth / (vinf^2);
b = abs(a) * sqrt(e^2 - 1);

h_vec = cross(r_vec,v_vec);
W_hat = h_vec / norm(h_vec);

% B vector
B_vec = b * cross(S_hat, W_hat);

% Project B onto R and T directions
BdotR = dot(B_vec, R_hat);
BdotT = dot(B_vec, T_hat);

end
