function [A, sortedjd, distance] = Get_A1_trangle_simple1111(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length,jd_grid,wd_grid,h_grid)
% Simplified version - Skips regional checks, directly calculates all rays.

% --- Calculate Grid Dimensions ---
nj = (jdmax - jdmin) / jdjg; % Number of longitude cells
nw = (wdmax - wdmin) / wdjg; % Number of latitude cells
ng = (gdmax - gdmin) / gdjg; % Number of height cells
nx = nj * nw * ng;           % Total number of cells in the 3D grid

% --- Get Number of Rays ---
n_rays = size(station, 1);

% fprintf('Starting to process %d rays...\n', n_rays);

% --- Set Up Parallel Pool (if not already running) ---
% if isempty(gcp('nocreate'))
%     parpool('local');
% end

% --- Use a single parfor loop to avoid complex slicing issues ---
% Pre-allocate cell arrays to store results from each worker
all_rows = cell(n_rays, 1);
all_cols = cell(n_rays, 1);
all_vals = cell(n_rays, 1);
non_zero_counts = zeros(n_rays, 1, 'uint32'); % Count of non-zero elements per ray

sortedjd_cell = cell(n_rays, 1); % To store sorted coordinate data for each ray
distance_cell = cell(n_rays, 1); % To store distance data for each ray

% fprintf('Progress: 0/%d', n_rays);
% wd_grid = [wd_grid;wd_grid];
for i = 1:n_rays
    % Get station and satellite coordinates for the current ray
    sx = station(i, 1);
    sy = station(i, 2);
    sz = station(i, 3);
    
    sat_x = satellite(i, 1);
    sat_y = satellite(i, 2);
    sat_z = satellite(i, 3);
    
    % Directly calculate the projection matrix - with error handling
    try
        [A1, sorted_ddzb, distances] = Get_sta_A1_tradition_optimized3(...
            sx, sy, sz, sat_x, sat_y, sat_z, ...
            jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length,nj,nw,ng,jd_grid,wd_grid,h_grid);
        
        if isempty(A1) 
            % Handle rays with no intersections
            all_rows{i} = [];
            all_cols{i} = [];
            all_vals{i} = [];
            non_zero_counts(i) = 0;
            sortedjd_cell{i} = [];
            distance_cell{i} = [];
        else
            % Extract non-zero elements from the sparse matrix
            [rows, cols, vals] = find(A1);
            all_rows{i} = uint32(rows);
            all_cols{i} = uint32(cols);
            all_vals{i} = single(vals); % Use single precision for memory efficiency
            non_zero_counts(i) = uint32(length(rows));
            sortedjd_cell{i} = sorted_ddzb';
            distance_cell{i} = distances';
        end
    catch ME
        % Error handling: skip problematic rays
        fprintf('Ray %d calculation failed: %s\n', i, ME.message);
        all_rows{i} = [];
        all_cols{i} = [];
        all_vals{i} = [];
        non_zero_counts(i) = 0;
        sortedjd_cell{i} = [];
        distance_cell{i} = [];
    end
    
    % Progress update (using a parallel-safe mechanism)
    % if mod(i, max(1, floor(n_rays/20))) == 0
    %     fprintf('\rProgress: %d/%d', i, n_rays);
    % end
end

fprintf('\rProgress: %d/%d Complete\n', n_rays, n_rays);

% --- Count Valid Rays ---
valid_mask = non_zero_counts > 0;
n_valid = sum(valid_mask);

% if n_valid == 0
%     A = sparse(0, nx);
%     sortedjd = struct();
%     distance = struct();
%     fprintf('Warning: No valid ray data found.\n');
%     return;
% end

% --- Pre-calculate Total Number of Non-Zero Elements ---
total_non_zeros = sum(double(non_zero_counts(valid_mask)));
total_rows = n_valid;

% --- Pre-allocate Arrays for Merging - use efficient data types ---
if total_non_zeros > 2e9
    % For very large arrays, use double to avoid indexing issues
    I = zeros(total_non_zeros, 1, 'double');
    J = zeros(total_non_zeros, 1, 'double');
else
    I = zeros(total_non_zeros, 1, 'uint32');
    J = zeros(total_non_zeros, 1, 'uint32');
end
V = zeros(total_non_zeros, 1, 'single');

% Pre-allocate arrays for struct data
all_sortedjd = cell(1, n_valid);
all_distance = cell(1, n_valid);

element_offset = 0;
valid_idx_counter = 1;
valid_indices = find(valid_mask);

% fprintf('Merging data: 0/%d', n_valid);

for i = 1:length(valid_indices)
    idx = valid_indices(i);
    
    n_current_non_zeros = double(non_zero_counts(idx));
    
    if n_current_non_zeros > 0
        current_rows = all_rows{idx};
        current_cols = all_cols{idx};
        current_vals = all_vals{idx};
        
        % Batch assignment for efficiency
        range = element_offset + 1 : element_offset + n_current_non_zeros;
        I(range) = valid_idx_counter;
        J(range) = current_cols;
        V(range) = current_vals;
        
        element_offset = element_offset + n_current_non_zeros;
    end
    
    % Handle struct data
    all_sortedjd{valid_idx_counter} = sortedjd_cell{idx};
    all_distance{valid_idx_counter} = distance_cell{idx};
    
    valid_idx_counter = valid_idx_counter + 1;
    
    % Release memory from processed cells immediately
    all_rows{idx} = [];
    all_cols{idx} = [];
    all_vals{idx} = [];
    sortedjd_cell{idx} = [];
    distance_cell{idx} = [];
    
    % Progress update
    % if mod(i, max(1, floor(length(valid_indices)/10))) == 0
    %     fprintf('\rMerging data: %d/%d', i, length(valid_indices));
    % end
end

fprintf('\rMerging data: %d/%d Complete\n', length(valid_indices), length(valid_indices));

% --- Create the Final Sparse Matrix ---
A = sparse(double(I), double(J), double(V), total_rows, nx);

% --- Generate Field Names and Create Output Structures ---
field_names = arrayfun(@(x) ['ray_', num2str(x)], 1:n_valid, 'UniformOutput', false);
sortedjd = cell2struct(all_sortedjd, field_names, 2);
distance = cell2struct(all_distance, field_names, 2);

% Final memory cleanup (optional, as function exit will handle it)
% clear all_rows all_cols all_vals sortedjd_cell distance_cell I J V;

end
