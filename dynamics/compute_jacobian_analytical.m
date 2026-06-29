function A = compute_jacobian_analytical(t, X, params)
% COMPUTE_JACOBIAN_ANALYTICAL Computes the analytical Jacobian df/dx
%
% Inputs:
%   t      - Current time (seconds from epoch)
%   X      - State vector [r; v; CR] (7x1)
%   params - Structure containing parameters
%
% Outputs:
%   A      - Jacobian matrix (7x7)

    % Extract state components
    r = X(1:3);  % Position vector (km)
    v = X(4:6);  % Velocity vector (km/s)
    CR = X(7);   % SRP coefficient

    % Extract parameters
    mu_earth = params.mu_earth;
    mu_sun = params.mu_sun;
    AU = params.AU;
    Phi_1AU = params.solar_flux;
    c_light = params.speed_of_light;
    area_mass_ratio = params.area_mass_ratio;

    % Current Julian Date
    JD = params.epoch + t/86400;

    % Get Earth's position relative to Sunnn from Ephem, file
    [r_earth_sun, ~, ~] = Ephem(JD, 3, 'EME2000');

    % Spacecraft position relative to Earth
    r_sc_earth = r;

    % Spacecraft position relative to Sun
    r_sc_sun = r_earth_sun + r_sc_earth;

    % Distances
    R_sc_earth = norm(r_sc_earth);
    R_sc_sun = norm(r_sc_sun);
    R_earth_sun = norm(r_earth_sun);
    delta = r_earth_sun - r_sc_sun;


    % Two-body gravity (Earth) partials
    I3 = eye(3);
    dadr_two_body = -mu_earth * ( (I3 / R_sc_earth^3) - 3*(r_sc_earth*r_sc_earth') / R_sc_earth^5 );

    % Third-body gravity (Sun) partials
    dadr_third_body = mu_sun * ( (I3 / norm(delta)^3) - 3*(delta*delta') / norm(delta)^5 );

    % SRP partials
    Phi = Phi_1AU * (AU / R_sc_sun)^2;
    P_phi = Phi / c_light;
    K_SRP = CR * P_phi * area_mass_ratio;
    dadr_SRP = -K_SRP * ( (I3 / R_sc_sun^3) - 3*(r_sc_sun*r_sc_sun') / R_sc_sun^5 );

    % Total partial wrt position
    dadr_total = dadr_two_body + dadr_third_body + dadr_SRP;

    % Partial wrt CR (only from SRP)
    dadCR = -P_phi * area_mass_ratio * (r_sc_sun / R_sc_sun);

    %  Assemble the Jacobian matrix
    A = zeros(7,7);

    A(1:3,4:6) = eye(3);          % dr/dv
    A(4:6,1:3) = dadr_total;      % dv/dr
    A(4:6,7)   = dadCR;           % dv/dCR

end
