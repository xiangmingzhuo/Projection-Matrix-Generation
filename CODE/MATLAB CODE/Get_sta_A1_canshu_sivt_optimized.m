% This function uses the parametric equation + effective splitting plane method
% to generate matrix A for a single ray.
function [one_dimensional_index, sortedjd_ddzb, distances] = Get_sta_A1_canshu_sivt_optimized(XD,YD,ZD,XS,YS,ZS,jd_params,wd_params,gd_params,total_length)
% Initialize output
[one_dimensional_index, sortedjd_ddzb, distances] = deal([]);
jdmin = jd_params(1); jdmax = jd_params(2); jdjg = jd_params(3);
wdmin = wd_params(1); wdmax = wd_params(2); wdjg = wd_params(3);
gdmin = gd_params(1); gdmax = gd_params(2); gdjg = gd_params(3);
% Calculate the first and last intersection points and determine the valid grid range
[is_valid, lon1, lon2, lat1, lat2, h1, h2, orientVector] = get_first_and_last_point3(XD,YD,ZD,XS,YS,ZS,jdjg,wdjg,jdmin,jdmax,wdmin,wdmax,gdmin,gdmax);

% If it's an invalid ray, exit the function immediately
if ~is_valid, return; end

% Pre-calculate common quantities
orient = orientVector(1,:);
R = 6371000; % Earth's radius in meters

%% Vectorized longitude plane calculation
% Generate longitude grid
lon_grid = deg2rad(lon1):deg2rad(jdjg):deg2rad(lon2);
if isempty(lon_grid), lon_grid = deg2rad(lon1); end  % Handle single-point case
% Vectorized calculation
% t_lon = (YD - XD*tan(lon_grid)) ./ ((XS - XD)*tan(lon_grid) - (YS - YD));
t_lon = (YD .* cos(lon_grid) - XD .* sin(lon_grid)) ./ (orientVector(1,1) .* sin(lon_grid) - orientVector(1,2) .* cos(lon_grid));
t_lon = t_lon(:);

%% Vectorized latitude plane calculation
% Generate latitude grid
lat_grid = deg2rad(lat1):deg2rad(wdjg):deg2rad(lat2);
if isempty(lat_grid), lat_grid = deg2rad(lat1); end % Handle single-point case
% Pre-calculate common terms
A_term = orient(3)^2 - (orient(1)^2 + orient(2)^2)*tan(lat_grid).^2;
B_term = 2*(ZD*orient(3) - (XD*orient(1) + YD*orient(2)).*tan(lat_grid).^2);
C_term = ZD^2 - (XD^2 + YD^2).*tan(lat_grid).^2;
% Pre-allocate memory
t_lat1 = zeros(size(lat_grid)) * NaN;
t_lat2 = zeros(size(lat_grid)) * NaN;
% Calculate valid solutions
t_lat1 = (-B_term - sqrt(B_term.^2 - 4*A_term.*C_term)) ./ (2*A_term);
t_lat2 = (-B_term + sqrt(B_term.^2 - 4*A_term.*C_term)) ./ (2*A_term);
% Flatten all solutions
t_lat = [t_lat1(:); t_lat2(:)];
t_lat = t_lat(imag(t_lat) == 0);
%% Vectorized altitude plane calculation
% Generate altitude grid
h_grid = h1:gdjg:h2;
if isempty(h_grid), h_grid = h1; end % Handle single-point case
% Pre-calculate coefficients
A = sum(orient.^2);
B = 2*(XD*orient(1) + YD*orient(2) + ZD*orient(3));
C_term = XD^2 + YD^2 + ZD^2 - (R + h_grid).^2;
% Pre-allocate memory
t_alt1 = zeros(size(h_grid)) * NaN;
t_alt2 = zeros(size(h_grid)) * NaN;
% Calculate valid solutions
t_alt1 = (-B - sqrt(B^2 - 4*A*C_term)) / (2*A);
t_alt2 = (-B + sqrt(B^2 - 4*A*C_term)) / (2*A);
% Organize results
t_alt = [t_alt1(:);t_alt2(:)];

%% Merge and filter logic
% Combine all t values and sort them
% 先过滤 t_all
t_all = [t_lon; t_lat; t_alt];
mask = isfinite(real(t_all)) & abs(imag(t_all)) < 1e-9;
t_all = real(t_all(mask));
t_all = sort(t_all);

% 计算交点
jd_coords_all = bsxfun(@plus, [XD, YD, ZD], bsxfun(@times, t_all, orient));

[lat, lon, alt] = XYZtoBLH_sphere( ...
    jd_coords_all(:,1), jd_coords_all(:,2), jd_coords_all(:,3), R);

jd_point_all = [rad2deg(lon), rad2deg(lat), alt];

tol_jd = max(abs(jdjg) * 1e-9, 1e-10);
tol_wd = max(abs(wdjg) * 1e-9, 1e-10);
tol_gd = max(abs(gdjg) * 1e-9, 1e-4);

valid_row = jd_point_all(:,1) >= jdmin - tol_jd & jd_point_all(:,1) <= jdmax + tol_jd & ...
            jd_point_all(:,2) >= wdmin - tol_wd & jd_point_all(:,2) <= wdmax + tol_wd & ...
            jd_point_all(:,3) >= gdmin - tol_gd & jd_point_all(:,3) <= gdmax + tol_gd;

t_valid = t_all(valid_row);
xyz_valid = jd_coords_all(valid_row, :);
jd_point_valid = jd_point_all(valid_row, :);

% 去掉重复交点，解决顶点/边共点问题
tol_t = max(1e-10, 1e-12 * max(abs(t_valid)));
keep = [true; abs(diff(t_valid)) > tol_t];

t_valid = t_valid(keep);
xyz_valid = xyz_valid(keep, :);
jd_point_valid = jd_point_valid(keep, :);

if length(t_valid) < 2
    return;
end

% segment length
distances = (diff(t_valid) * norm(orient))';

% 用 segment 内部点算所属 voxel
mid_xyz = 0.5 * (xyz_valid(1:end-1, :) + xyz_valid(2:end, :));

[lat_m, lon_m, alt_m] = XYZtoBLH_sphere( ...
    mid_xyz(:,1), mid_xyz(:,2), mid_xyz(:,3), R);

mid_blh = [rad2deg(lon_m), rad2deg(lat_m), alt_m]';

[voxel_indices, grid_indices] = voxel_index_from_inside_points( ...
    mid_blh, jd_params, wd_params, gd_params);

sortedjd_ddzb = round(jd_point_valid', 13);

% one_dimensional_index = sparse(1, grid_indices, distances, 1, total_length);
one_dimensional_index = sparse(grid_indices, 1, distances, total_length, 1);
end
