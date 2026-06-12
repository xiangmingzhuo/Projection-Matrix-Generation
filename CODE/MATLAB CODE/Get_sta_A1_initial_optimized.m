% This function uses equation solving combined with the effective splitting plane method
% to generate matrix A for a single ray.
function [one_dimensional_index,sortedjd_ddzb,values] = Get_sta_A1_initial_optimized(XD,YD,ZD,XS,YS,ZS,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,total_length,nj,nw,ng)

%%% Initialize output
[one_dimensional_index, sortedjd_ddzb, values] = deal([]);

% Calculate the first and last intersection points to determine the effective grid range
% [lon1,lon2,lat1,lat2,is_valid] = get_first_and_last_point2(jdmin, jdmax, wdmin, wdmax, gdmin, gdmax,XD,YD,ZD,XS,YS,ZS,jdjg,wdjg);
[is_valid, lon1, lon2, lat1, lat2] = get_first_and_last_point3(XD,YD,ZD,XS,YS,ZS,jdjg,wdjg,jdmin,jdmax,wdmin,wdmax,gdmin,gdmax);

% If the ray is not valid, exit the function immediately
if ~is_valid
    return;
end 

R = 6371000;

jdmin_ = jdmin;
wdmin_ = wdmin;
jdjg_ = jdjg;
wdjg_ = wdjg;

jdmax_ = jdmax;
wdmax_ = wdmax;

jdmin = deg2rad(jdmin);
jdjg = deg2rad(jdjg);
jdmax = deg2rad(jdmax);

wdmin = deg2rad(wdmin);
wdjg = deg2rad(wdjg);
wdmax = deg2rad(wdmax);

syms X Y Z

f1 = (X-XD)/(XS-XD) - (Y-YD)/(YS-YD);
f2 = (X-XD)/(XS-XD) - (Z-ZD)/(ZS-ZD);
f3 = (Y-YD)/(YS-YD) - (Z-ZD)/(ZS-ZD);

%% ============================================================
%  Optimized part 1: height surface intersections
%  Original idea is retained: solve is still used explicitly.
%  Difference: solve only once with symbolic height H.
%% ============================================================

X0_coordinates = [];
Y0_coordinates = [];
Z0_coordinates = [];

height_range = gdmin:gdjg:gdmax;
num_heights = length(height_range);

if num_heights > 0
    syms H
    f4 = X^2 + Y^2 + Z^2 - (R+H)^2;

    solution_H = solve(f1, f2, f3, f4, X, Y, Z);

    if ~isempty(solution_H.X)
        XH_fun = matlabFunction(solution_H.X(:), 'Vars', H);
        YH_fun = matlabFunction(solution_H.Y(:), 'Vars', H);
        ZH_fun = matlabFunction(solution_H.Z(:), 'Vars', H);

        num_H_solutions = numel(solution_H.X);

        X0_coordinates = zeros(num_H_solutions, num_heights);
        Y0_coordinates = zeros(num_H_solutions, num_heights);
        Z0_coordinates = zeros(num_H_solutions, num_heights);

        for idx = 1:num_heights
            current_H = height_range(idx);

            X0_coordinates(:,idx) = XH_fun(current_H);
            Y0_coordinates(:,idx) = YH_fun(current_H);
            Z0_coordinates(:,idx) = ZH_fun(current_H);
        end
    end
end

%% ============================================================
%  Optimized part 2: longitude plane intersections
%  Original idea is retained: solve is still used explicitly.
%  Difference: solve only once with symbolic longitude J.
%% ============================================================

X1_coordinates = [];
Y1_coordinates = [];
Z1_coordinates = [];

longitude_range = deg2rad(lon1):jdjg:deg2rad(lon2);
num_longitudes = length(longitude_range);

if num_longitudes > 0
    syms J
    f6 = tan(J)*X - Y;  % Equation of the longitude plane

    solution_J = solve(f1, f2, f3, f6, X, Y, Z);

    if ~isempty(solution_J.X)
        XJ_fun = matlabFunction(solution_J.X(:), 'Vars', J);
        YJ_fun = matlabFunction(solution_J.Y(:), 'Vars', J);
        ZJ_fun = matlabFunction(solution_J.Z(:), 'Vars', J);

        num_J_solutions = numel(solution_J.X);

        X1_coordinates = zeros(num_J_solutions, num_longitudes);
        Y1_coordinates = zeros(num_J_solutions, num_longitudes);
        Z1_coordinates = zeros(num_J_solutions, num_longitudes);

        for idx = 1:num_longitudes
            current_J = longitude_range(idx);

            X1_coordinates(:,idx) = XJ_fun(current_J);
            Y1_coordinates(:,idx) = YJ_fun(current_J);
            Z1_coordinates(:,idx) = ZJ_fun(current_J);
        end
    end
end

%% ============================================================
%  Optimized part 3: latitude plane intersections
%  Original idea is retained: solve is still used explicitly.
%  Difference: solve only once with symbolic latitude W.
%% ============================================================

X2_coordinates = [];
Y2_coordinates = [];
Z2_coordinates = [];

latitude_range = deg2rad(lat1):wdjg:deg2rad(lat2);
num_latitudes = length(latitude_range);

if num_latitudes > 0
    syms W
    f5 = tan(W)^2*(X^2 + Y^2) - Z^2;  % Equation of the latitude plane

    solution_W = solve(f1, f2, f3, f5, X, Y, Z);

    if ~isempty(solution_W.X)
        XW_fun = matlabFunction(solution_W.X(:), 'Vars', W);
        YW_fun = matlabFunction(solution_W.Y(:), 'Vars', W);
        ZW_fun = matlabFunction(solution_W.Z(:), 'Vars', W);

        num_W_solutions = numel(solution_W.X);

        X2_coordinates = zeros(num_W_solutions, num_latitudes);
        Y2_coordinates = zeros(num_W_solutions, num_latitudes);
        Z2_coordinates = zeros(num_W_solutions, num_latitudes);

        for idx = 1:num_latitudes
            current_W = latitude_range(idx);

            X2_coordinates(:,idx) = XW_fun(current_W);
            Y2_coordinates(:,idx) = YW_fun(current_W);
            Z2_coordinates(:,idx) = ZW_fun(current_W);
        end
    end
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
valid_wd = JD_ddzb(1,:) >= rad2deg(wdmin) - epsilon & JD_ddzb(1,:) <= rad2deg(wdmax) + epsilon;

% Check if the second row (longitude) is between jdmin and jdmax
valid_jd = JD_ddzb(2,:) >= rad2deg(jdmin) - epsilon & JD_ddzb(2,:) <= rad2deg(jdmax) + epsilon;

% Check if the third row (altitude) is between gdmin and gdmax
valid_gd = JD_ddzb(3,:) >= gdmin - epsilon & JD_ddzb(3,:) <= gdmax + epsilon;

% Combine all conditions to find indices of valid columns
valid_columns = valid_wd & valid_jd & valid_gd;

% Exit the function if there are no valid intersection points
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

if sortedjd_ddzb(3,1) > gdmin + 1e-5 || sortedjd_ddzb(3,end) < gdmax - 1e-5
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
grid_indices = find_grid_indices(midpoint_blh, jdmin_, wdmin_, gdmin, jdjg_, wdjg_, gdjg, jdmax_, wdmax_, gdmax, epsilon,nj,nw,ng);

%% Sparse matrix optimization
one_dimensional_index = sparse(grid_indices, 1, values, total_length, 1);

end
