% This function uses equation solving combined with the effective splitting plane method
% to generate matrix A for a single ray.
function [one_dimensional_index,sortedjd_ddzb,values] = Get_sta_A1_initial_optimized(XD,YD,ZD,XS,YS,ZS,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,total_length)
%%% Initialize output
[one_dimensional_index, sortedjd_ddzb, values] = deal([]);
% Calculate the first and last intersection points to determine the effective grid range
[lon1,lon2,lat1,lat2,is_valid] = get_first_and_last_point2(jdmin, jdmax, wdmin, wdmax, gdmin, gdmax,XD,YD,ZD,XS,YS,ZS,jdjg,wdjg);
% If the ray is not valid, exit the function immediately
if ~is_valid, return; end 
% Ground coordinates, satellite coordinates, longitude range and interval (in degrees),
% latitude range and interval (in degrees), altitude range and interval, mean Earth radius of the inversion area
R=6371000;
jdmin_ = jdmin;wdmin_ = wdmin;
jdjg_ = jdjg;wdjg_ = wdjg;
jdmax_ = jdmax;wdmax_ = wdmax;
jdmin = deg2rad(jdmin);
jdjg = deg2rad(jdjg);
jdmax = deg2rad(jdmax);
wdmin = deg2rad(wdmin);
wdjg = deg2rad(wdjg);
wdmax = deg2rad(wdmax);
syms X Y Z
f1 = (X-XD)/(XS-XD)-(Y-YD)/(YS-YD);
f2 = (X-XD)/(XS-XD)-(Z-ZD)/(ZS-ZD);
f3 = (Y-YD)/(YS-YD)-(Z-ZD)/(ZS-ZD);
% Initialize result storage arrays
X0_coordinates = [];
Y0_coordinates = [];
Z0_coordinates = [];
% Generate height sequence
height_range = gdmin:gdjg:gdmax;
% Pre-allocate cell arrays
num_heights = length(height_range);
X_cell = cell(1, num_heights);
Y_cell = cell(1, num_heights);
Z_cell = cell(1, num_heights);
% Calculation loop
for idx = 1:num_heights
    current_H = height_range(idx);
    f4 = X^2 + Y^2 + Z^2 - (R+current_H)^2;
    % Solve the system of equations
    solution = solve(f1, f2, f3, f4, X, Y, Z);
    % Store results
    X_cell{idx} = double(solution.X);
    Y_cell{idx} = double(solution.Y);
    Z_cell{idx} = double(solution.Z);
end
% Concatenate all results
for idx = 1:num_heights
    X0_coordinates = [X0_coordinates, X_cell{idx}];
    Y0_coordinates = [Y0_coordinates, Y_cell{idx}];
    Z0_coordinates = [Z0_coordinates, Z_cell{idx}];
end

% Initialize longitude plane intersection result arrays
X1_coordinates = [];
Y1_coordinates = [];
Z1_coordinates = [];
% Generate longitude sequence
longitude_range = deg2rad(lon1):jdjg:deg2rad(lon2);
% Pre-allocate cell arrays
num_longitudes = length(longitude_range);
X_cell_long = cell(1, num_longitudes);
Y_cell_long = cell(1, num_longitudes);
Z_cell_long = cell(1, num_longitudes);
% Calculate longitude plane intersections
for idx = 1:num_longitudes
    current_J = longitude_range(idx);
    f6 = tan(current_J)*X - Y;  % Longitude plane equation 
    % Solve the system of equations
    solution = solve(f1, f2, f3, f6, X, Y, Z);
    % Store results
    X_cell_long{idx} = double(solution.X);
    Y_cell_long{idx} = double(solution.Y);
    Z_cell_long{idx} = double(solution.Z);
end
% Concatenate all longitude plane results
for idx = 1:num_longitudes
    X1_coordinates = [X1_coordinates, X_cell_long{idx}];
    Y1_coordinates = [Y1_coordinates, Y_cell_long{idx}];
    Z1_coordinates = [Z1_coordinates, Z_cell_long{idx}];
end

% Initialize latitude plane intersection result arrays
X2_coordinates = [];
Y2_coordinates = [];
Z2_coordinates = [];
% Generate latitude sequence
latitude_range = deg2rad(lat1):wdjg:deg2rad(lat2);
% Pre-allocate cell arrays
num_latitudes = length(latitude_range);
X_cell_lat = cell(1, num_latitudes);
Y_cell_lat = cell(1, num_latitudes);
Z_cell_lat = cell(1, num_latitudes);
% Calculate latitude plane intersections
for idx = 1:num_latitudes
    current_W = latitude_range(idx);
    f5 = tan(current_W)^2*(X^2 + Y^2) - Z^2;  % Latitude plane equation
    % Solve the system of equations
    solution = solve(f1, f2, f3, f5, X, Y, Z);
    % Handle potential empty solutions
    if ~isempty(solution.X)
        X_cell_lat{idx} = double(solution.X);
        Y_cell_lat{idx} = double(solution.Y); 
        Z_cell_lat{idx} = double(solution.Z);
    else
        X_cell_lat{idx} = [];
        Y_cell_lat{idx} = [];
        Z_cell_lat{idx} = [];
    end
end
% Concatenate all latitude plane results
for idx = 1:num_latitudes
    X2_coordinates = [X2_coordinates, X_cell_lat{idx}];
    Y2_coordinates = [Y2_coordinates, Y_cell_lat{idx}];
    Z2_coordinates = [Z2_coordinates, Z_cell_lat{idx}];
end

JD_ddzb01 = [X0_coordinates(1,:);Y0_coordinates(1,:);Z0_coordinates(1,:)]';
JD_ddzb02 = [X0_coordinates(2,:);Y0_coordinates(2,:);Z0_coordinates(2,:)]';
JD_ddzb1 = [X1_coordinates(1,:);Y1_coordinates(1,:);Z1_coordinates(1,:)]';
JD_ddzb21 = [X2_coordinates(1,:);Y2_coordinates(1,:);Z2_coordinates(1,:)]';
JD_ddzb22 = [X2_coordinates(2,:);Y2_coordinates(2,:);Z2_coordinates(2,:)]';
JD_ddzb = [JD_ddzb01;JD_ddzb02;JD_ddzb1;JD_ddzb21;JD_ddzb22];
real_rows = all(imag(JD_ddzb) == 0, 2);  % Find rows that are entirely real
JD_ddzb = JD_ddzb(real_rows, :);         % Extract these rows
% Batch coordinate conversion
[lat, lon, alt] = XYZtoBLH_sphere(JD_ddzb(:,1), JD_ddzb(:,2), JD_ddzb(:,3),R);
JD_ddzb = [rad2deg(lat),rad2deg(lon),alt]';
epsilon = 1e-6;
% Check if the first row (latitude) is between wdmin and wdmax
valid_wd = JD_ddzb(1,:) >= rad2deg(wdmin) -epsilon & JD_ddzb(1,:) <= rad2deg(wdmax) +epsilon;
% Check if the second row (longitude) is between jdmin and jdmax
valid_jd = JD_ddzb(2,:) >= rad2deg(jdmin) - epsilon & JD_ddzb(2,:) <= rad2deg(jdmax) +epsilon;
% Check if the third row (altitude) is between gdmin and gdmax
valid_gd = JD_ddzb(3,:) >= gdmin -epsilon & JD_ddzb(3,:) <= gdmax +epsilon;
% Combine all conditions to find indices of valid columns
valid_columns = valid_wd & valid_jd & valid_gd;
% If there are no valid intersection points, exit the function
if isempty(valid_columns)
    return;
end
% Get the matrix of valid points based on column indices
valid_points = JD_ddzb(:,valid_columns)';
if isempty(valid_points)
    sortedjd_ddzb = [];
    return;
end
% Convert coordinates
[X,Y,Z] = BLHtoXYZ_sphere(deg2rad(valid_points(:,1)),deg2rad(valid_points(:,2)),valid_points(:,3),R);
sortedjd = [X, Y, Z]';
sortedjd_ddzb = valid_points';

%%% Calculate distances
distances = sqrt((sortedjd(1,:) - XD).^2 + (sortedjd(2,:) - YD).^2 + (sortedjd(3,:) - ZD).^2);

% Sort the distances
[distances, idx] = sort(distances);
% Reorder the sortedjd matrix
sortedjd = sortedjd(:, idx);
sortedjd_ddzb = sortedjd_ddzb(:, idx);
if sortedjd_ddzb(3,1)> gdmin + 1e-10  || sortedjd_ddzb(3,end)< gdmax - 1e-10
    sortedjd_ddzb = [];
    return;
end   
%%% Calculate one-dimensional indices
if length(distances) < 2
    return;
end

% Calculate the midpoint coordinates between each pair of intersection points
midpoints = (sortedjd(:,1:end-1) + sortedjd(:,2:end))/2;
% Batch convert midpoint coordinates
[lat_mid, lon_mid, alt_mid] = XYZtoBLH_sphere(midpoints(1,:), midpoints(2,:), midpoints(3,:),R);
midpoint_blh = [rad2deg(lon_mid)', rad2deg(lat_mid)', alt_mid'];
% Calculate the intercept values between each pair of intersection points
values = distances(2:end) - distances(1:end-1);
%% Index calculation optimization
epsilon = 1e-15;
grid_indices = find_grid_indices(midpoint_blh, jdmin_, wdmin_, gdmin, jdjg_, wdjg_, gdjg, jdmax_, wdmax_, gdmax, epsilon);
%% Sparse matrix optimization
% total_length = round(((rad2deg(wdmax)-rad2deg(wdmin))/rad2deg(wdjg)) * ((rad2deg(jdmax)-rad2deg(jdmin))/rad2deg(jdjg)) * ((gdmax-gdmin)/gdjg));
one_dimensional_index = sparse(1, grid_indices, values, 1, total_length);
