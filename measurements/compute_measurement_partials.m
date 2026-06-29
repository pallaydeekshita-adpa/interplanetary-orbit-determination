function [dp_dX, dp_dY, dp_dZ, dp_dXdot, dp_dYdot, dp_dZdot, ...
          dpdot_dX, dpdot_dY, dpdot_dZ, dpdot_dXdot, dpdot_dYdot, dpdot_dZdot, ...
          dp_dX_s, dp_dY_s, dp_dZ_s, dpdot_dX_s, dpdot_dY_s, dpdot_dZ_s] = compute_measurement_partials(state, station_id,w_earth)

    % Extract spacecraft position and velocity
    X = state(1);  Y = state(2);  Z = state(3);
    Xdot = state(4);  Ydot = state(5);  Zdot = state(6);

    % Use a switch case to assign station positions based on the station ID
    switch station_id
        case 101
            % Station 1 position: Extract position from state vector (elements 10:12)
            station_pos = state(10:12);
            
        case 337
            % Station 337 position: Extract position from state vector (elements 13:15)
            station_pos = state(13:15);
            
        case 394
            % Station 394 position: Extract position from state vector (elements 16:18)
            station_pos = state(16:18);
            
        otherwise
            error('Unknown station ID: %d', station_id);  % Handle unknown station IDs
    end

 
    X_s = station_pos(1);
    Y_s = station_pos(2);
    Z_s = station_pos(3);

    Xdot_s = - w_earth * Y_s;
    Ydot_s =  w_earth * X_s;
    Zdot_s = 0;

    % Compute range
    rho = sqrt((X - X_s)^2 + (Y - Y_s)^2 + (Z - Z_s)^2);

    % Compute range rate
    rho_dot = ((X - X_s) * (Xdot - Xdot_s) + (Y - Y_s) * (Ydot - Ydot_s) + (Z - Z_s) * (Zdot - Zdot_s)) / rho;
    
    % Compute Jacobians
    dp_dX = (X - X_s) / rho;
    dp_dY = (Y - Y_s) / rho;
    dp_dZ = (Z - Z_s) / rho;
    dp_dXdot = 0;
    dp_dYdot = 0;
    dp_dZdot = 0;

    dpdot_dX = (Xdot - Xdot_s) / rho - rho_dot * (X - X_s) / rho^2;
    dpdot_dY = (Ydot - Ydot_s) / rho - rho_dot * (Y - Y_s) / rho^2;
    dpdot_dZ = (Zdot - Zdot_s) / rho - rho_dot * (Z - Z_s) / rho^2;
    dpdot_dXdot = (X - X_s) / rho;
    dpdot_dYdot = (Y - Y_s) / rho;
    dpdot_dZdot = (Z - Z_s) / rho;

    % Partials w.r.t. station positions
    dp_dX_s = -dp_dX;
    dp_dY_s = -dp_dY;
    dp_dZ_s = -dp_dZ;
    
    dpdot_dX_s = ((2*X - 2*X_s)*(Zdot*(Z - Z_s) + (X - X_s)*(Xdot + Y_s*w_earth) + (Y - Y_s)*(Ydot - X_s*w_earth)))/(2*((X - X_s)^2 + (Y - Y_s)^2 + (Z - Z_s)^2)^(3/2)) - (Xdot + Y*w_earth)/((X - X_s)^2 + (Y - Y_s)^2 + (Z - Z_s)^2)^(1/2);
    dpdot_dY_s = ((2*Y - 2*Y_s)*(Zdot*(Z - Z_s) + (X - X_s)*(Xdot + Y_s*w_earth) + (Y - Y_s)*(Ydot - X_s*w_earth)))/(2*((X - X_s)^2 + (Y - Y_s)^2 + (Z - Z_s)^2)^(3/2)) - (Ydot - X*w_earth)/((X - X_s)^2 + (Y - Y_s)^2 + (Z - Z_s)^2)^(1/2);
    dpdot_dZ_s = ((2*Z - 2*Z_s)*(Zdot*(Z - Z_s) + (X - X_s)*(Xdot + Y_s*w_earth) + (Y - Y_s)*(Ydot - X_s*w_earth)))/(2*((X - X_s)^2 + (Y - Y_s)^2 + (Z - Z_s)^2)^(3/2)) - Zdot/((X - X_s)^2 + (Y - Y_s)^2 + (Z - Z_s)^2)^(1/2);

end
