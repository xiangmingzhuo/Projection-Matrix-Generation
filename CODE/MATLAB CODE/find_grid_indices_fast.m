function grid_indices = find_grid_indices_fast(mid_lon, mid_lat, mid_h, ...
    jdmin, wdmin, gdmin, jdjg, wdjg, gdjg, epsilon, Nx, Ny)

    indices_x = ceil((mid_lon(:) + epsilon - jdmin) / jdjg);
    indices_y = ceil((mid_lat(:) + epsilon - wdmin) / wdjg);
    indices_z = ceil((mid_h(:)   + epsilon - gdmin) / gdjg);

    indices_x(indices_x == 0) = 1;
    indices_x(indices_x > Nx) = Nx;

    indices_y(indices_y == 0) = 1;
    indices_y(indices_y > Ny) = Ny;

    grid_indices = (indices_z - 1) * (Nx * Ny) + (indices_y - 1) * Nx + indices_x;
end
