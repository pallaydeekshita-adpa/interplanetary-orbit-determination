function H_tilde = Htilde_6080(state, obs, t, options)
% Computes the measurement sensitivity matrix Htilde
% Inputs:
%   state - spacecraft state [7x1]: (x,y,z,vx,vy,vz,CR)
%   obs - observations structure
%   t - current time (seconds)
%   options - structure with station info and w_earth
%
% Output:
%   H_tilde - [2x7] partials of measurements w.r.t state

% Find observation index
index = find(obs.times == t);

if isempty(index)
    error('Htilde_6080: Time %f not found in observations.', t);
end

station_id = obs.station(index);

% Extract spacecraft state
x = state(1);
y = state(2);
z = state(3);
vx = state(4);
vy = state(5);
vz = state(6);

% Get station ECEF position
switch station_id
    case 34
        station_ecef = options.extra_args.stations_ecef.Canberra;
    case 65
        station_ecef = options.extra_args.stations_ecef.Madrid;
    case 13
        station_ecef = options.extra_args.stations_ecef.Goldstone;
    otherwise
        error('Htilde_6080: Unknown station ID %d.', station_id);
end

% Rotate station position into ECI
theta = options.extra_args.w_earth * t;
R = [cos(theta) -sin(theta) 0;
     sin(theta)  cos(theta) 0;
     0            0         1];

station_eci = R * station_ecef;
station_vel_eci = cross([0; 0; options.extra_args.w_earth], station_eci);

% Relative position and velocity
rx = x - station_eci(1);
ry = y - station_eci(2);
rz = z - station_eci(3);

vx_rel = vx - station_vel_eci(1);
vy_rel = vy - station_vel_eci(2);
vz_rel = vz - station_vel_eci(3);

rho = sqrt(rx^2 + ry^2 + rz^2);

% ----------------- Range partials -----------------
drho_dx = rx / rho;
drho_dy = ry / rho;
drho_dz = rz / rho;

% No partial w.r.t velocity
drho_dvx = 0;
drho_dvy = 0;
drho_dvz = 0;

% No partial w.r.t CR
drho_dCR = 0;

% ----------------- Range-rate partials -----------------
% First term: (v_rel - (rho_dot / rho) * r_rel)
rho_dot = (rx*vx_rel + ry*vy_rel + rz*vz_rel) / rho;

drhodot_dx = (vx_rel / rho) - (rho_dot * rx) / rho^2;
drhodot_dy = (vy_rel / rho) - (rho_dot * ry) / rho^2;
drhodot_dz = (vz_rel / rho) - (rho_dot * rz) / rho^2;

drhodot_dvx = rx / rho;
drhodot_dvy = ry / rho;
drhodot_dvz = rz / rho;

drhodot_dCR = 0;

% ----------------- Assemble H_tilde -----------------
H_tilde = zeros(2,7);

% Row 1: range partials
H_tilde(1,1) = drho_dx;
H_tilde(1,2) = drho_dy;
H_tilde(1,3) = drho_dz;
H_tilde(1,4) = drho_dvx;
H_tilde(1,5) = drho_dvy;
H_tilde(1,6) = drho_dvz;
H_tilde(1,7) = drho_dCR;

% Row 2: range-rate partials
H_tilde(2,1) = drhodot_dx;
H_tilde(2,2) = drhodot_dy;
H_tilde(2,3) = drhodot_dz;
H_tilde(2,4) = drhodot_dvx;
H_tilde(2,5) = drhodot_dvy;
H_tilde(2,6) = drhodot_dvz;
H_tilde(2,7) = drhodot_dCR;

end
