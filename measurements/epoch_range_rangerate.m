function C = epoch_range_rangerate(obs, state, t, options)
% Computes range and range-rate from spacecraft to ground station
% Inputs:
%   obs - observation structure
%   state - spacecraft state [7x1]: (x,y,z,vx,vy,vz,CR)
%   t - current time (seconds)
%   options - structure with station info and w_earth
%
% Output:
%   C = [range, range_rate]

% Extract spacecraft position and velocity
sat_pos = state(1:3)';  % [km]
sat_vel = state(4:6)';  % [km/s]

% Find index of current time
index = find(abs(obs.times - t) < 1e-6, 1);
if isempty(index)
    error('epoch_range_rangerate: Time %.6f not found in observations.', t);
end


station_id = obs.station(index);

% Extract station ECEF position
switch station_id
    case 34
        station_ecef = options.extra_args.stations_ecef.Canberra;
    case 65
        station_ecef = options.extra_args.stations_ecef.Madrid;
    case 13
        station_ecef = options.extra_args.stations_ecef.Goldstone;
    otherwise
        error('epoch_range_rangerate: Unknown station ID %d.', station_id);
end

% Rotate station position to ECI
theta = options.extra_args.w_earth * t;
R = [cos(theta) -sin(theta) 0;
     sin(theta)  cos(theta) 0;
     0            0         1];

station_eci = R * station_ecef;

% Station velocity in ECI (due to Earth rotation)
w_earth_vec = [0; 0; options.extra_args.w_earth];
station_vel_eci = cross(w_earth_vec, station_eci);

% Compute range
rho_vec = sat_pos - station_eci;
range = norm(rho_vec);

% Compute range-rate
range_rate = dot(rho_vec, sat_vel - station_vel_eci) / range;

% Output
C = [range, range_rate];
end
