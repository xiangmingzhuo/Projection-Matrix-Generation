function X = generateAndSortGridPoints(lon_all, lat_all, alt_all)
% GENERATEANDSORTGRIDPOINTS Creates a sorted matrix of all 3D grid points.
%
% This function takes 1D vectors for longitude, latitude, and altitude,
% generates the full 3D grid of all possible combinations, and then
% sorts the resulting points in a specific order.
%
% Inputs:
%   lon_all - A 1D vector of longitude values.
%   lat_all - A 1D vector of latitude values.
%   alt_all - A 1D vector of altitude values.
%
% Output:
%   X - An N-by-3 matrix, where each row is a [lon, lat, alt] point.
%     The matrix is sorted primarily by altitude, then latitude, then longitude.

% Create the 3D grid from the 1D vectors.
% meshgrid replicates the vectors to create 3D arrays where each corresponding
% element (i,j,k) represents a single point in the grid.
[lon_grid, lat_grid, alt_grid] = meshgrid(lon_all, lat_all, alt_all);

% Reshape the 3D grid arrays into single column vectors.
% The (:) operator stacks the columns of an array into one long column.
lon_col = lon_grid(:);
lat_col = lat_grid(:);
alt_col = alt_grid(:);

% Combine the three column vectors into an N-by-3 matrix.
% Each row now represents a unique [lon, lat, alt] point from the grid.
X = [lon_col, lat_col, alt_col];

% Sort the rows of the matrix.
% The sort order is specified by the vector [3, 2, 1], which means:
% 1. Sort primarily by the 3rd column (altitude).
% 2. For rows with the same altitude, sort by the 2nd column (latitude).
% 3. For rows with the same altitude and latitude, sort by the 1st column (longitude).
X = sortrows(X, [3, 2, 1]);
end
