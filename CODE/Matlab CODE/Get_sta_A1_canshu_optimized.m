% This function uses the parametric equation + effective splitting plane method
% to generate matrix A for a single ray.
function [one_dimensional_index, sortedjd_ddzb, distances] = Get_sta_A1_canshu_optimized(XD,YD,ZD,XS,YS,ZS,jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length,nj,nw,ng)
% Initialize output
[one_dimensional_index, sortedjd_ddzb, distances] = deal([]);

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
t_all = [t_lon; t_lat; t_alt];
t_sorted = sort(t_all);
% Early termination check
% if isempty(t_sorted)
%     return;
% end

%% Batch coordinate conversion
% Batch calculate intersection point coordinates
jd_coords = [XD,YD,ZD] + orient .* t_sorted;
jd_coords = jd_coords(imag(jd_coords(:,1)) == 0,:);
% Batch coordinate conversion
[lat, lon, alt] = XYZtoBLH_sphere(jd_coords(:,1), jd_coords(:,2), jd_coords(:,3),R);
jd_point = [rad2deg(lon), rad2deg(lat), alt];

%% Vectorized validity check
valid_row = (jd_point(:,1) >= jdmin) & (jd_point(:,1) <= jdmax) & ...
            (jd_point(:,2) >= wdmin) & (jd_point(:,2) <= wdmax) & ...
            (jd_point(:,3) >= gdmin - 1e-4) & (jd_point(:,3) <= gdmax + 1e-4);
% Apply filter
t_valid = t_sorted(valid_row);
% if length(t_valid) < 2  % Need at least two points to calculate distance
%     return;
% end

%% Distance calculation optimization
distances = (diff(t_valid) * norm(orient))';
sortedjd = jd_coords(valid_row,:)';
sortedjd_ddzb = jd_point(valid_row,:)';
sortedjd_ddzb(2,:) = round(sortedjd_ddzb(2,:),10);
sortedjd_ddzb(3,:) = round(sortedjd_ddzb(3,:));
%% Midpoint calculation optimization
midpoint = (sortedjd(:,1:end-1) + sortedjd(:,2:end)) / 2;
% Batch convert midpoint coordinates
[lat_mid, lon_mid, alt_mid] = XYZtoBLH_sphere(midpoint(1,:), midpoint(2,:), midpoint(3,:),R);
midpoint_blh = [rad2deg(lon_mid)', rad2deg(lat_mid)', alt_mid'];
midpoint_blh(:,2) = round(midpoint_blh(:,2),10);
%% Index calculation optimization
epsilon = 1e-15;
grid_indices = find_grid_indices(midpoint_blh, jdmin, wdmin, gdmin, jdjg, wdjg, gdjg, jdmax, wdmax, gdmax, epsilon,nj,nw,ng);

%% Sparse matrix optimization
% total_length = round(((wdmax-wdmin)/wdjg) * ((jdmax-jdmin)/jdjg) * ((gdmax-gdmin)/gdjg));
% one_dimensional_index = sparse(1, grid_indices, distances, 1, total_length);
one_dimensional_index = sparse(double(grid_indices), 1, distances, total_length, 1);
end
