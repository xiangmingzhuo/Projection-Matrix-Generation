function [voxel_indices, grid_indices] = voxel_index_from_inside_points(points, jd_params, wd_params, gd_params)

    jd_min = jd_params(1); jd_max = jd_params(2); jd_spacing = jd_params(3);
    wd_min = wd_params(1); wd_max = wd_params(2); wd_spacing = wd_params(3);
    gd_min = gd_params(1); gd_max = gd_params(2); gd_spacing = gd_params(3);

    Nx = round((jd_max - jd_min) / jd_spacing);
    Ny = round((wd_max - wd_min) / wd_spacing);
    Nz = round((gd_max - gd_min) / gd_spacing);

    i = coord_to_cell(points(1,:), jd_min, jd_spacing, Nx);
    j = coord_to_cell(points(2,:), wd_min, wd_spacing, Ny);
    k = coord_to_cell(points(3,:), gd_min, gd_spacing, Nz);

    voxel_indices = [i; j; k];
    grid_indices = (k - 1) * (Nx * Ny) + (j - 1) * Nx + i;

    bad = i < 1 | i > Nx | j < 1 | j > Ny | k < 1 | k > Nz;
    if any(bad)
        error('SIVT index out of range: check boundary tolerance or intersection filtering.');
    end
end

function idx = coord_to_cell(x, xmin, dx, N)
    tol = max(abs(dx) * 1e-9, 1e-12);
    idx = floor((x - xmin + tol) ./ dx) + 1;
    idx = min(max(idx, 1), N);
end
