function [A, sortedjd, distance] = Get_A1_trangle_simple2(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length)
% Get_A1_trangle_simple2 - Simplified version that skips region checks and directly computes all rays.
% This function is optimized for performance using parallel processing and efficient memory allocation.

% --- Input Parameters ---
% station:     [N x 3] matrix of station coordinates (x, y, z).
% satellite:   [N x 3] matrix of satellite coordinates (x, y, z).
% jdmin, jdmax, jdjg: Longitude min, max, and interval (grid step).
% wdmin, wdmax, wdjg: Latitude min, max, and interval (grid step).
% gdmin, gdmax, gdjg: Depth/Height min, max, and interval (grid step).
%
% --- Output Parameters ---
% A:           Sparse projection matrix (N_valid x N_grid_cells).
% sortedjd:    Structure array where each field 'shexianX' contains the sorted grid intersection points for ray X.
% distance:    Structure array where each field 'shexianX' contains the corresponding distances for ray X.

% Calculate grid dimensions
nj = (jdmax - jdmin) / jdjg;
nw = (wdmax - wdmin) / wdjg;
ng = (gdmax - gdmin) / gdjg;
nx = nj * nw * ng; % Total number of grid cells

% Get the number of rays
n_rays = size(station, 1);

% fprintf('Starting to process %d rays...\n', n_rays);

% Set up parallel pool (if not already started)
% if isempty(gcp('nocreate'))
%     parpool('local');
% end

% Use a single parfor loop to avoid complex slicing issues
all_rows = cell(n_rays, 1);
all_cols = cell(n_rays, 1);
all_vals = cell(n_rays, 1);
non_zero_counts = zeros(n_rays, 1, 'uint32');

sortedjd_cell = cell(n_rays, 1);
distance_cell = cell(n_rays, 1);

% Progress tracking
completed = 0;
% fprintf('Progress: 0/%d', n_rays);

for i = 1:n_rays
    % Get station and satellite coordinates for the current ray
    sx = station(i, 1);
    sy = station(i, 2);
    sz = station(i, 3);
    
    sat_x = satellite(i, 1);
    sat_y = satellite(i, 2);
    sat_z = satellite(i, 3);
    
    % Directly compute the projection matrix - with error handling
    try
        [A1, sorted_ddzb, distances] = Get_sta_A1_optimized(...
            sx, sy, sz, sat_x, sat_y, sat_z, ...
            jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length,nj,nw,ng);
        
        if isempty(A1)
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
            all_vals{i} = single(vals);
            non_zero_counts(i) = uint32(length(rows));
            sortedjd_cell{i} = sorted_ddzb';
            distance_cell{i} = distances';
        end
    catch ME
        % Error handling: skip problematic rays
        fprintf('Ray %d computation failed: %s\n', i, ME.message);
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

% Count valid rays
valid_mask = non_zero_counts > 0;
n_valid = sum(valid_mask);

% if n_valid == 0
%     A = sparse(0, nx);
%     sortedjd = struct();
%     distance = struct();
%     fprintf('Warning: No valid ray data found.\n');
%     return;
% end

% Pre-calculate the total number of non-zero elements
total_non_zeros = sum(double(non_zero_counts(valid_mask)));
total_rows = n_valid;

% Pre-allocate arrays - using more efficient data types
if total_non_zeros > 2e9
    % For very large arrays, use double to avoid indexing issues
    I = zeros(total_non_zeros, 1, 'double');
    J = zeros(total_non_zeros, 1, 'double');
else
    I = zeros(total_non_zeros, 1, 'uint32');
    J = zeros(total_non_zeros, 1, 'uint32');
end
V = zeros(total_non_zeros, 1, 'single');

% Pre-allocate structure data
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
        
        % Batch assignment
        range = element_offset + 1 : element_offset + n_current_non_zeros;
        I(range) = valid_idx_counter;
        J(range) = current_cols;
        V(range) = current_vals;
        
        element_offset = element_offset + n_current_non_zeros;
    end
    
    % Handle structure data
    all_sortedjd{valid_idx_counter} = sortedjd_cell{idx};
    all_distance{valid_idx_counter} = distance_cell{idx};
    
    valid_idx_counter = valid_idx_counter + 1;
    
    % Release memory promptly
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

% Create the final sparse matrix
A = sparse(double(I), double(J), double(V), total_rows, nx);

% Generate field names and create the output structures
field_names = arrayfun(@(x) ['shexian', num2str(x)], 1:n_valid, 'UniformOutput', false);
sortedjd = cell2struct(all_sortedjd, field_names, 2);
distance = cell2struct(all_distance, field_names, 2);

% Final memory cleanup (optional, as function exit will handle it)
% clear all_rows all_cols all_vals sortedjd_cell distance_cell I J V;

end
