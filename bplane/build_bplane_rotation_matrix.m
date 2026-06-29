function C = build_bplane_rotation_matrix(r, v, mu_earth)
% Returns rotation matrix C: columns are [R_hat, T_hat, S_hat]
% Use in:  P_bplane = C' * P_eci * C

    % Define incoming v-infinity direction
    S_hat = v / norm(v);

    % Reference normal (solar system north)
    N_hat = [0; 0; 1];

    % Construct B-plane frame
    T_hat = cross(S_hat, N_hat);
    T_hat = T_hat / norm(T_hat);

    R_hat = cross(T_hat, S_hat);
    R_hat = R_hat / norm(R_hat);  % might not be unit length

    C = [R_hat, T_hat, S_hat];  % [BdotR, BdotT, S]

    det(C)
end
