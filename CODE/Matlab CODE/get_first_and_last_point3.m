function [is_valid, lon1, lon2, lat1, lat2, h1, h2, orientVector] = get_first_and_last_point3 (XD, YD, ZD, XS, YS, ZS, jdjg, wdjg, jdmin, jdmax, wdmin, wdmax, gdmin, gdmax)
% Initialize valid range
lon1 = []; lon2 = []; lat1 = []; lat2 = []; h1 = []; h2 = [];

% Calculate ray direction vector
orientVector = [XS,YS,ZS] - [XD,YD,ZD]; % Calculate the direction vector

% Calculate the two intersection points with the top and bottom height planes
R = 6371000; % Earth's radius in meters
A = sum(orientVector(1,:).^2);  % More concise calculation of the sum of squares
B = 2 * dot([XD, YD, ZD], orientVector(1,:));  % Simplify calculation using dot product
Cmin = XD^2 + YD^2 + ZD^2 - (R + gdmin)^2;
Cmax = XD^2 + YD^2 + ZD^2 - (R + gdmax)^2;

% Unify the calculation of tmin and tmax
discriminant_min = B^2 - 4*A*Cmin;
discriminant_max = B^2 - 4*A*Cmax;

% Use vectorized calculation and logical indexing
t_solutions_min = [(-B - sqrt(discriminant_min)); (-B + sqrt(discriminant_min))] / (2*A);
t_solutions_max = [(-B - sqrt(discriminant_max)); (-B + sqrt(discriminant_max))] / (2*A);

% Filter for valid solutions (t must be between 0 and 1)
tmin = t_solutions_min((t_solutions_min > 0) & (t_solutions_min < 1));
tmax = t_solutions_max((t_solutions_max > 0) & (t_solutions_max < 1));
if isempty(tmin) || isempty(tmax)
    is_valid = false;
    return;
end
% Calculate the XYZ coordinates of the entry and exit points
intersection_point = [XD,YD,ZD] + orientVector .* tmin; % Entry point
outersection_point = [XD,YD,ZD] + orientVector .* tmax;  % Exit point

% Calculate the lat, lon, and height values for the entry and exit points
% Note: XYZtoBLH_sphere returns [lat, lon, height]
[enter_point(:,2), enter_point(:,1), enter_point(:,3)] = XYZtoBLH_sphere(intersection_point(:,1), intersection_point(:,2), intersection_point(:,3),R);
enter_point  = [rad2deg(enter_point(:,1)),rad2deg(enter_point(:,2)),enter_point(:,3)];
[exit_point(:,2), exit_point(:,1), exit_point(:,3)] = XYZtoBLH_sphere(outersection_point(:,1), outersection_point(:,2), outersection_point(:,3),R);
exit_point  = [rad2deg(exit_point(:,1)),rad2deg(exit_point(:,2)),exit_point(:,3)];

% Check if the intersection points are within the grid boundaries. If not, the ray is considered invalid.
if enter_point(1,2) >= wdmin && enter_point(1,2) <= wdmax && enter_point(1,1) >= jdmin && enter_point(1,1) <= jdmax
    is_valid = true;
else
    is_valid = false;
    return;
end
if exit_point(1,2) >= wdmin && exit_point(1,2) <= wdmax && exit_point(1,1) >= jdmin && exit_point(1,1) <= jdmax
    is_valid = true;
else
    is_valid = false;
    return;
end

% Calculate the final bounding box for the ray, snapping to the grid and expanding by one grid step
lon1 = max(approximateNumberDown(min(enter_point(:,1),exit_point(:,1)),jdjg)-jdjg,jdmin);
lon2 = min(approximateNumberUp(max(enter_point(:,1),exit_point(:,1)),jdjg)+jdjg,jdmax);
lat1 = max(approximateNumberDown(min(enter_point(:,2),exit_point(:,2)),wdjg)-wdjg,wdmin);
lat2 = min(approximateNumberUp(max(enter_point(:,2),exit_point(:,2)),wdjg)+wdjg,wdmax);
h1 = gdmin;
h2 = gdmax;
end