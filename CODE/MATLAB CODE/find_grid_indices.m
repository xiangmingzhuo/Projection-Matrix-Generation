function grid_indices = find_grid_indices(midpoint_blh, jdmin, wdmin, gdmin, jdjg, wdjg, gdjg, jdmax, wdmax, gdmax, epsilon,Nx,Ny,Nz)
    indices_x = ceil((midpoint_blh(:, 1) + epsilon - jdmin) / jdjg);
    indices_y = ceil((midpoint_blh(:, 2) + epsilon - wdmin) / wdjg);
    indices_z = ceil((midpoint_blh(:, 3) + epsilon - gdmin) / gdjg);
    indices_x(indices_x == 0) = 1;indices_x(indices_x > Nx) = Nx;
    indices_y(indices_y == 0) = 1;indices_y(indices_y > Ny) = Ny;
    % indices_x = floor((midpoint_blh(:, 1) - epsilon - jdmin) / jdjg)+ 1;
    % indices_y = floor((midpoint_blh(:, 2) - epsilon - wdmin) / wdjg)+ 1;
    % indices_z = floor((midpoint_blh(:, 3) - epsilon - gdmin) / gdjg)+ 1;

    grid_indices = (indices_z - 1) * (Nx * Ny) + (indices_y - 1) * Nx + indices_x;

    % tol = 1e-12;
    % 
    % indices_x = floor((midpoint_blh(:,1) - jdmin) / jdjg) + 1;
    % indices_y = floor((midpoint_blh(:,2) - wdmin) / wdjg) + 1;
    % indices_z = floor((midpoint_blh(:,3) - gdmin) / gdjg) + 1;
    % 
    % valid_mid = indices_x >= 1 & indices_x <= Nx & ...
    %             indices_y >= 1 & indices_y <= Ny & ...
    %             indices_z >= 1 & indices_z <= Nz;
    % 
    % indices_x = indices_x(valid_mid);
    % indices_y = indices_y(valid_mid);
    % indices_z = indices_z(valid_mid);
    % 
    % grid_indices = (indices_z - 1) * (Nx * Ny) + ...
    %                (indices_y - 1) * Nx + indices_x;
end
