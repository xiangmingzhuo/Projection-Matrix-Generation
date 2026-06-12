function [voxel_indices, grid_indices] = SIVT_voxel_tracing(ray_intersections, jd_params, wd_params, gd_params)
% SIVT算法中的voxel归属索引确定
% 输入:
%   ray_intersections: 3×n矩阵，[经度; 纬度; 高度]，已按到测站距离排序
%   jd_params: 经度参数 [min, max, spacing]
%   wd_params: 纬度参数 [min, max, spacing]  
%   gd_params: 高度参数 [min, max, spacing]
% 输出:
%   voxel_indices: 3×n矩阵，每个交点的voxel索引 [i; j; k]

    % 提取参数
    jd_min = jd_params(1); jd_max = jd_params(2); jd_spacing = jd_params(3);
    wd_min = wd_params(1); wd_max = wd_params(2); wd_spacing = wd_params(3);
    gd_min = gd_params(1); gd_max = gd_params(2); gd_spacing = gd_params(3);
    
    % 计算每个方向的网格数
    Nx = (jd_max - jd_min) / jd_spacing;
    Ny = (wd_max - wd_min) / wd_spacing;
    Nz = (gd_max - gd_min) / gd_spacing;

    % 生成 splitting surfaces
    jd_surfaces = jd_min:jd_spacing:jd_max;  % 经度分割面
    wd_surfaces = wd_min:wd_spacing:wd_max;  % 纬度分割面  
    gd_surfaces = gd_min:gd_spacing:gd_max;  % 高度分割面
    
    n_intersections = size(ray_intersections, 2);
    voxel_indices = zeros(3, n_intersections-1);
    grid_indices = zeros(1, n_intersections-1);
    
    % 步骤1: 确定第一个交点的voxel索引
    first_point = ray_intersections(:, 1);
    voxel_indices(:, 1) = find_voxel_index(first_point, jd_surfaces, wd_surfaces, gd_surfaces);
    grid_indices(1) = (voxel_indices(3, 1) - 1) * (Nx * Ny) + (voxel_indices(2, 1) - 1) * Nx + voxel_indices(1, 1);
    
    % 如果没有更多交点，直接返回
    if n_intersections <= 2
        return;
    end
    
    % 步骤2: 向量化计算所有后续交点
    % 提取所有前点和当前点
    prev_points = ray_intersections(:, 1:end-2);  % 第1到n-2个点
    current_points = ray_intersections(:, 2:end-1); % 第2到n-1个点
    
    % 向量化计算方向向量
    [I, J, K] = compute_local_orientation_vector_vectorized(prev_points, current_points);
    
    % 向量化确定穿过的分割面类型
    crossed_directions = identify_crossed_surfaces_with_direction_vectorized(...
        prev_points, current_points, jd_surfaces, wd_surfaces, gd_surfaces, I, J, K);
    
    % 使用累积和计算所有后续索引
    cumulative_changes = cumsum(crossed_directions, 2);
    voxel_indices(:, 2:end) = voxel_indices(:, 1) + cumulative_changes;
    
    % 向量化计算grid_indices
    i_indices = voxel_indices(1, :);
    j_indices = voxel_indices(2, :);
    k_indices = voxel_indices(3, :);
    grid_indices = (k_indices - 1) * (Nx * Ny) + (j_indices - 1) * Nx + i_indices;
    
    % fprintf('SIVT voxel tracing完成:\n');
    % fprintf('处理了 %d 个交点\n', n_intersections);
    % fprintf('得到 %d 个线段\n', n_intersections-1);
end

function [I, J, K] = compute_local_orientation_vector_vectorized(prev_points, current_points)
% 向量化计算方向向量
    delta = current_points - prev_points;
    
    % 经度方向 (考虑-180°到180°的跳变)
    lon_diff = delta(1, :);
    I = sign(lon_diff);
    % 处理经度跳变
    jump_mask = abs(lon_diff) > 180;
    I(jump_mask) = sign(-lon_diff(jump_mask));
    
    % 纬度方向
    lat_diff = delta(2, :);
    J = sign(lat_diff);
    
    % 高度方向  
    height_diff = delta(3, :);
    K = sign(height_diff);
    
    % 处理零值情况
    I(I == 0) = 1;
    J(J == 0) = -1;  % 对于纬度，如果差值为零，设为-1（向南）
end

function crossed_directions = identify_crossed_surfaces_with_direction_vectorized(...
    prev_points, current_points, jd_surfaces, wd_surfaces, gd_surfaces, I, J, K)
% 向量化确定射线穿过的分割面类型

    n_segments = size(prev_points, 2);
    crossed_directions = zeros(3, n_segments);
    tolerance = 1e-10;
    
    % 经度方向
    if any(I ~= 0)
        prev_jd = prev_points(1, :);
        curr_jd = current_points(1, :);
        
        % 为每个线段计算索引
        prev_jd_idx = find_voxel_index_for_coordinate_vectorized(prev_jd, jd_surfaces);
        curr_jd_idx = find_voxel_index_for_coordinate_vectorized(curr_jd, jd_surfaces);
        
        % 应用容差调整
        for i = 1:n_segments
            if I(i) == -1
                prev_jd_idx(i) = find_voxel_index_for_coordinate(prev_jd(i) - 0.5*tolerance, jd_surfaces);
                curr_jd_idx(i) = find_voxel_index_for_coordinate(curr_jd(i) - 0.5*tolerance, jd_surfaces);
            else
                prev_jd_idx(i) = find_voxel_index_for_coordinate(prev_jd(i) + 0.5*tolerance, jd_surfaces);
                curr_jd_idx(i) = find_voxel_index_for_coordinate(curr_jd(i) + 0.5*tolerance, jd_surfaces);
            end
        end
        
        % 确定索引变化
        idx_change = curr_jd_idx - prev_jd_idx;
        change_mask = idx_change ~= 0;
        crossed_directions(1, change_mask) = I(change_mask);
        
        % 处理正好在分割面上的情况
        on_surface_mask = ~change_mask & is_point_on_surface_vectorized(curr_jd, jd_surfaces, tolerance);
        crossed_directions(1, on_surface_mask) = I(on_surface_mask);
    end
    
    % 纬度方向
    if any(J ~= 0)
        prev_wd = prev_points(2, :);
        curr_wd = current_points(2, :);
        
        % 为每个线段计算索引
        prev_wd_idx = find_voxel_index_for_coordinate_vectorized(prev_wd, wd_surfaces);
        curr_wd_idx = find_voxel_index_for_coordinate_vectorized(curr_wd, wd_surfaces);
        
        % 应用容差调整
        for i = 1:n_segments
            if J(i) == -1
                prev_wd_idx(i) = find_voxel_index_for_coordinate(prev_wd(i) - 0.5*tolerance, wd_surfaces);
                curr_wd_idx(i) = find_voxel_index_for_coordinate(curr_wd(i) - 0.5*tolerance, wd_surfaces);
            else
                prev_wd_idx(i) = find_voxel_index_for_coordinate(prev_wd(i) + 0.5*tolerance, wd_surfaces);
                curr_wd_idx(i) = find_voxel_index_for_coordinate(curr_wd(i) + 0.5*tolerance, wd_surfaces);
            end
        end
        
        % 确定索引变化
        idx_change = curr_wd_idx - prev_wd_idx;
        change_mask = idx_change ~= 0;
        crossed_directions(2, change_mask) = J(change_mask);
        
        % 处理正好在分割面上的情况
        on_surface_mask = ~change_mask & is_point_on_surface_vectorized(curr_wd, wd_surfaces, tolerance);
        crossed_directions(2, on_surface_mask) = J(on_surface_mask);
    end
    
    % 高度方向
    if any(K ~= 0)
        prev_gd = prev_points(3, :);
        curr_gd = current_points(3, :);
        
        % 为每个线段计算索引
        prev_gd_idx = find_voxel_index_for_coordinate_vectorized(prev_gd, gd_surfaces);
        curr_gd_idx = find_voxel_index_for_coordinate_vectorized(curr_gd, gd_surfaces);
        
        % 确定索引变化
        idx_change = curr_gd_idx - prev_gd_idx;
        change_mask = idx_change ~= 0;
        crossed_directions(3, change_mask) = K(change_mask);
        
        % 处理正好在分割面上的情况
        on_surface_mask = ~change_mask & is_point_on_surface_vectorized(curr_gd, gd_surfaces, tolerance);
        crossed_directions(3, on_surface_mask) = K(on_surface_mask);
    end
end

function indices = find_voxel_index_for_coordinate_vectorized(values, surfaces)
% 向量化查找坐标值在分割面数组中的索引
    n_values = length(values);
    indices = zeros(1, n_values);
    for i = 1:n_values
        indices(i) = find_voxel_index_for_coordinate(values(i), surfaces);
    end
end

function is_on_surface = is_point_on_surface_vectorized(values, surfaces, tolerance)
% 向量化检查点是否正好在分割面上
    n_values = length(values);
    is_on_surface = false(1, n_values);
    for i = 1:n_values
        is_on_surface(i) = is_point_on_surface(values(i), surfaces, tolerance);
    end
end

% 保留原有的标量函数
function [I, J, K] = compute_local_orientation_vector(prev_point, current_point)
% 计算当前线段的局部方向向量 (I, J, K)
    % 经度方向 (考虑-180°到180°的跳变)
    lon_diff = current_point(1) - prev_point(1);
    if abs(lon_diff) > 180
        I = sign(-lon_diff);  % 处理经度跳变
    else
        I = sign(lon_diff);
    end
    
    % 纬度方向
    lat_diff = current_point(2) - prev_point(2);
    J = sign(lat_diff);
    
    % 高度方向  
    height_diff = current_point(3) - prev_point(3);
    K = sign(height_diff);
    
    % 处理零值情况
    if I == 0, I = 1; end
    if J == 0, J = -1; end  % 对于纬度，如果差值为零，设为-1（向南）
end

function voxel_idx = find_voxel_index(point, jd_surfaces, wd_surfaces, gd_surfaces)
% 通过二分查找确定点在voxel模型中的索引
    % 经度索引
    i = find(point(1) >= jd_surfaces, 1, 'last');
    if isempty(i), i = 1; end
    if i == length(jd_surfaces), i = i - 1; end
    
    % 纬度索引  
    j = find(point(2) >= wd_surfaces, 1, 'last');
    if isempty(j), j = 1; end
    if j == length(wd_surfaces), j = j - 1; end
    
    % 高度索引
    k = find(point(3) >= gd_surfaces, 1, 'last'); 
    if isempty(k), k = 1; end
    if k == length(gd_surfaces), k = k - 1; end
    
    voxel_idx = [i; j; k];
end

function idx = find_voxel_index_for_coordinate(value, surfaces)
% 查找坐标值在分割面数组中的索引
    idx = find(value >= surfaces, 1, 'last');
    if isempty(idx), idx = 1; end
    if idx == length(surfaces), idx = idx - 1; end
end

function is_on_surface = is_point_on_surface(value, surfaces, tolerance)
% 检查点是否正好在分割面上
    distances = abs(value - surfaces);
    is_on_surface = any(distances < tolerance);
end