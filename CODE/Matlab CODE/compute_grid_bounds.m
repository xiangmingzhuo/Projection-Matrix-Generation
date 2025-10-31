function [Bmin, Bmax] = compute_grid_bounds(lon_min, lon_max, lat_min, lat_max, h_min, h_max)
% COMPUTE_GRID_BOUNDS Calculates the Cartesian bounding box for a geographic grid.
%
% This function provides a fast approximation of the Cartesian bounding box
% for a geographic region defined by its minimum/maximum longitude, latitude,
% and height. It operates under the assumption that the minimum Cartesian
% coordinates are derived from the minimum geographic coordinates and the
% maximum Cartesian coordinates from the maximum geographic coordinates.
%
% IMPORTANT CAVEAT: This assumption is NOT generally true for a spherical
% or ellipsoidal Earth model. For example, the point with the maximum
% X-coordinate is often at (lat=0, lon=90) or (lat=0, lon=-90), not at the
% maximum latitude of the grid. This function is therefore only a valid
% and fast approximation if the grid is very small and does not cross
% the equator or the prime meridian. For a robust and accurate calculation,
% one must convert all 8 corners of the grid and find the true min/max.
%
% Inputs:
%   lon_min, lon_max - Minimum and maximum longitude in degrees.
%   lat_min, lat_max - Minimum and maximum latitude in degrees.
%   h_min, h_max     - Minimum and maximum height in meters above the
%                      reference sphere/ellipsoid.
%
% Outputs:
%   Bmin - 3-element vector [X_min, Y_min, Z_min] of the bounding box.
%   Bmax - 3-element vector [X_max, Y_max, Z_max] of the bounding box.

% --- Simplified Calculation ---
% The code below calculates the bounds by only converting the corner
% defined by the minimums and the corner defined by the maximums.
% This is fast but not a true AABB (Axis-Aligned Bounding Box) for a
% spherical grid.

% Radius of the reference sphere in meters.
R = 6371000;

% Convert heights to radial distance from the Earth's center.
h_min = R + h_min;
h_max = R + h_max;

% Calculate the 2D Cartesian coordinates for the min and max corners.
% This assumes a helper function 'latlon_to_cartesian' exists.
% NOTE: This function is not a standard MATLAB function. A typical
% implementation would be:
% x = r * cosd(lat) * cosd(lon);
% y = r * cosd(lat) * sind(lon);
% The height (r) is used here, which is unusual for a 2D conversion.
% We will assume the user's function handles this correctly.
[Bmin(1), Bmin(2)] = latlon_to_cartesian(lat_min, lon_min, h_min);
[Bmax(1), Bmax(2)] = latlon_to_cartesian(lat_max, lon_max, h_max);

% Assign the Z-bounds directly from the radial heights.
Bmin(3) = h_min;
Bmax(3) = h_max;

end
