function [X, Y, Z] = BLHtoXYZ_sphere(B, L, H, R)
%BLHTOXYZ_SPHERE Converts geographic (BLH) coordinates to Cartesian (XYZ) on a sphere.
%
%   [X, Y, Z] = BLHtoXYZ_sphere(B, L, H, R) converts geodetic coordinates
%   (latitude B, longitude L, and height H) to 3D Cartesian coordinates (X, Y, Z)
%   assuming a spherical Earth model.
%
%   INPUTS:
%       B - Latitude in RADIANS. Can be a scalar or an array.
%       L - Longitude in RADIANS. Can be a scalar or an array.
%       H - Height above the sphere's surface (in the same units as R).
%       R - Radius of the sphere.
%
%   OUTPUTS:
%       X - Cartesian X-coordinate.
%       Y - Cartesian Y-coordinate.
%       Z - Cartesian Z-coordinate.

% Calculate the total radial distance from the sphere's center.
% This makes the subsequent calculations cleaner and slightly more efficient.
N = (R + H);

% Calculate the Cartesian X-coordinate
X = N .* cos(B) .* cos(L);

% Calculate the Cartesian Y-coordinate
Y = N .* cos(B) .* sin(L);

% Calculate the Cartesian Z-coordinate
Z = N .* sin(B);

end
