function grid_indices = find_grid_indices(midpoint_blh, jdmin, wdmin, gdmin, jdjg, wdjg, gdjg, jdmax, wdmax, gdmax, epsilon,Nx,Ny,Nz)
    indices_x = ceil((midpoint_blh(:, 1) + epsilon - jdmin) / jdjg);
    indices_y = ceil((midpoint_blh(:, 2) + epsilon - wdmin) / wdjg);
    indices_z = ceil((midpoint_blh(:, 3) + epsilon - gdmin) / gdjg);

    grid_indices = (indices_z - 1) * (Nx * Ny) + (indices_y - 1) * Nx + indices_x;
end
