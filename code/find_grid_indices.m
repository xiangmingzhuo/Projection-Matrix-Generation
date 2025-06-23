function grid_indices = find_grid_indices(midpoint_blh, jdmin, wdmin, gdmin, jdjg, wdjg, gdjg, jdmax, wdmax, gdmax,epsilon)
    % 计算三个方向的网格数
    Nx = (jdmax - jdmin)/jdjg ;
    Ny = (wdmax - wdmin)/wdjg ;
    Nz = (gdmax - gdmin)/gdjg ;
    % 计算每个点相对于网格起始点的索引
    indices_x = ceil((midpoint_blh(:, 1) + epsilon - jdmin) / jdjg);
    indices_y = ceil((midpoint_blh(:, 2) + epsilon - wdmin) / wdjg);
    indices_z = ceil((midpoint_blh(:, 3) + epsilon - gdmin) / gdjg);
    
    % 确保索引不超出网格范围
    % indices_x(indices_x > Nx) = Nx;
    % indices_y(indices_y > Ny) = Ny;
    % indices_z(indices_z > Nz) = Nz;
    
    % 将索引转换为网格编号
    % 注意：MATLAB中矩阵索引从1开始
    % grid_indices = sub2ind([Nx, Ny, Nz], indices_x, indices_y, indices_z);
    % 计算一维索引编号
    grid_indices = (indices_z -1) * (Nx * Ny) + (indices_y-1 ) * Nx + indices_x;
end
