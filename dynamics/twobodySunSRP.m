function dX = twobodySunSRP(t, X, params)
% TWO-BODY + Sun Gravity + Solar Radiation Pressure (SRP) with 7x7 STM
%
% Inputs:
%   t      - Time since epoch (seconds)
%   X      - State vector [56 x 1]: [position (3), velocity (3), CR (1), STM (7x7 flattened)]
%   params - Structure containing constants (mu_earth, mu_sun, AU, etc.)
%
% Output:
%   dX     - Derivative of X (56 x 1)

% ----------------- Extract Parameters -----------------
mu_earth = params.mu_earth;
mu_sun = params.mu_sun;
AU = params.AU;
solar_flux = params.solar_flux;
c_light = params.speed_of_light;
area_mass_ratio = params.area_mass_ratio;
epoch_JD = params.epoch;


% % --- Maneuver Parameters ---
% t_maneuver = 5236 * 3600;  % seconds
% dv_applied = [0.4; 0; 0];    
% maneuver_duration = 10;   % seconds buffer window to apply ΔV once

% ----------------- Unpack State -----------------
x  = X(1);
y  = X(2);
z  = X(3);
vx = X(4);
vy = X(5);
vz = X(6);
CR = X(7);

% persistent dv_applied_once;
% if isempty(dv_applied_once)
%     dv_applied_once = false;
% end
% 
% if ~dv_applied_once && t >= t_maneuver
%     vx = vx + dv_applied(1);
%     vy = vy + dv_applied(2);
%     vz = vz + dv_applied(3);
% 
%     fprintf('🔔 ΔV applied at t = %.2f hr (%.2f s)\n', t/3600, t);
%     dv_applied_once = true;
% end


STM = reshape(X(8:56), 7, 7);  % 7x7 STM

% ----------------- Current Julian Date -----------------
JD = epoch_JD + t/86400;  % seconds to days

% ----------------- Get Planetary Positions -----------------
[r_sun_to_earth, ~, ~] = Ephem(JD, 3, 'EME2000');

% ----------------- Spacecraft Position Vectors -----------------
r_earth_to_sc = [x; y; z];
r_sun_to_sc   = r_sun_to_earth + r_earth_to_sc;

% ----------------- Norms and Unit Vectors -----------------
rE_norm  = norm(r_earth_to_sc);
rS_norm  = norm(r_sun_to_sc);
rSE_norm = norm(r_sun_to_earth);

rhat_earth_to_sc = r_earth_to_sc / rE_norm;
rhat_sun_to_sc   = r_sun_to_sc / rS_norm;

% ----------------- Accelerations -----------------
% Earth gravity
a_earth = -mu_earth * r_earth_to_sc / rE_norm^3;

% Sun third-body gravity
a_sun = -mu_sun * (r_sun_to_sc / rS_norm^3 - r_sun_to_earth / rSE_norm^3);

% Solar Radiation Pressure (SRP)
Phi = solar_flux * (AU / rS_norm)^2;    % Solar flux scaled with distance
P_phi = Phi / c_light;                  % Pressure [N/m²] -> [km²/s²]
a_srp = CR * P_phi * area_mass_ratio * rhat_sun_to_sc;

% Total acceleration
a_total = a_earth + a_sun + a_srp;

% ----------------- Jacobian Partials -----------------
% Earth two-body partials
dax_dx_earth = -mu_earth * (1/rE_norm^3 - 3*x^2/rE_norm^5);
dax_dy_earth = 3 * mu_earth * x * y / rE_norm^5;
dax_dz_earth = 3 * mu_earth * x * z / rE_norm^5;

day_dx_earth = dax_dy_earth;
day_dy_earth = -mu_earth * (1/rE_norm^3 - 3*y^2/rE_norm^5);
day_dz_earth = 3 * mu_earth * y * z / rE_norm^5;

daz_dx_earth = dax_dz_earth;
daz_dy_earth = day_dz_earth;
daz_dz_earth = -mu_earth * (1/rE_norm^3 - 3*z^2/rE_norm^5);

% Sun third-body partials
xs = r_sun_to_sc(1);
ys = r_sun_to_sc(2);
zs = r_sun_to_sc(3);

dax_dx_sun = -mu_sun * (1/rS_norm^3 - 3*xs^2/rS_norm^5);
dax_dy_sun = 3 * mu_sun * xs * ys / rS_norm^5;
dax_dz_sun = 3 * mu_sun * xs * zs / rS_norm^5;

day_dx_sun = dax_dy_sun;
day_dy_sun = -mu_sun * (1/rS_norm^3 - 3*ys^2/rS_norm^5);
day_dz_sun = 3 * mu_sun * ys * zs / rS_norm^5;

daz_dx_sun = dax_dz_sun;
daz_dy_sun = day_dz_sun;
daz_dz_sun = -mu_sun * (1/rS_norm^3 - 3*zs^2/rS_norm^5);

% SRP partials
K_SRP = CR * P_phi * area_mass_ratio;

dax_dx_srp = K_SRP * (1/rS_norm^3 - 3*xs^2/rS_norm^5);
dax_dy_srp = -3 * K_SRP * xs * ys / rS_norm^5;
dax_dz_srp = -3 * K_SRP * xs * zs / rS_norm^5;

day_dx_srp = -3 * K_SRP * ys * xs / rS_norm^5;
day_dy_srp = K_SRP * (1/rS_norm^3 - 3*ys^2/rS_norm^5);
day_dz_srp = -3 * K_SRP * ys * zs / rS_norm^5;

daz_dx_srp = -3 * K_SRP * zs * xs / rS_norm^5;
daz_dy_srp = -3 * K_SRP * zs * ys / rS_norm^5;
daz_dz_srp = K_SRP * (1/rS_norm^3 - 3*zs^2/rS_norm^5);

% ----------------- Total Partial Derivatives -----------------
dax_dx = dax_dx_earth + dax_dx_sun + dax_dx_srp;
dax_dy = dax_dy_earth + dax_dy_sun + dax_dy_srp;
dax_dz = dax_dz_earth + dax_dz_sun + dax_dz_srp;

day_dx = day_dx_earth + day_dx_sun + day_dx_srp;
day_dy = day_dy_earth + day_dy_sun + day_dy_srp;
day_dz = day_dz_earth + day_dz_sun + day_dz_srp;

daz_dx = daz_dx_earth + daz_dx_sun + daz_dx_srp;
daz_dy = daz_dy_earth + daz_dy_sun + daz_dy_srp;
daz_dz = daz_dz_earth + daz_dz_sun + daz_dz_srp;

% Partial w.r.t CR (only from SRP)
dadCR = P_phi * area_mass_ratio * rhat_sun_to_sc;

% ----------------- Build the A Matrix (7x7) -----------------
A = zeros(7,7);

% position derivatives
A(1,4) = 1;
A(2,5) = 1;
A(3,6) = 1;

% acceleration partials
A(4,1) = dax_dx;
A(4,2) = dax_dy;
A(4,3) = dax_dz;
A(4,7) = dadCR(1);

A(5,1) = day_dx;
A(5,2) = day_dy;
A(5,3) = day_dz;
A(5,7) = dadCR(2);

A(6,1) = daz_dx;
A(6,2) = daz_dy;
A(6,3) = daz_dz;
A(6,7) = dadCR(3);

% CR row stays zeros because CR is constant

% ----------------- STM Propagation -----------------
Phi_dot = A * STM;

% ----------------- Assemble Derivative -----------------
dX = zeros(56,1);

% State derivative
dX(1:3) = [vx; vy; vz];
dX(4:6) = a_total;
dX(7)   = 0;              % CR is constant
dX(8:56) = reshape(Phi_dot,49,1); % Flatten 7x7 Phi_dot into column

end
