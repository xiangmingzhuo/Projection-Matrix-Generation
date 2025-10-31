function grid_indices = find_grid_indices(midpoint_blh, jdmin, wdmin, gdmin, jdjg, wdjg, gdjg, jdmax, wdmax, gdmax, epsilon,Nx,Ny,Nz)
% FIND_GRID_INDICES Converts 3D BLH coordinates to a single linear grid index.
%
%   This function calculates the linear index for a set of points within a
%   3D rectangular grid defined by its minimum/maximum bounds and grid
%   spacing. It uses a manual calculation for maximum performance, avoiding
%   the overhead of MATLAB's built-in sub2ind function.
%
%   The grid is assumed to be axis-aligned in the BLH (Longitude, Latitude,
%   Height) coordinate space.
%
% Inputs:
%   midpoint_blh   - An M-by-3 matrix, where each row is a point [lon, lat, h].
%   jdmin, wdmin, gdmin - The minimum longitude, latitude, and height of the grid.
%   jdjg, wdjg, gdjg   - The grid spacing (interval) for longitude, latitude, and height.
%   jdmax, wdmax, gdmax - The maximum longitude, latitude, and height of the grid.
%   epsilon        - A small value to handle floating-point inaccuracies and
%                    ensure points on the upper boundaries are correctly
%                    indexed into the last cell.
%
% Outputs:
%   grid_indices   - An M-by-1 column vector of linear grid indices.

    % Calculate the total number of cells along each grid dimension.
    % This is effectively the size of the grid.
    % Nx = (jdmax - jdmin) / jdjg;
    % Ny = (wdmax - wdmin) / wdjg;
    % Nz = (gdmax - gdmin) / gdjg;

    % Calculate the integer index for each point along each dimension.
    % The 'ceil' function ensures that points are assigned to the cell
    % they fall into. The 'epsilon' is crucial for correctly handling points
    % that lie exactly on a grid boundary, especially the upper boundary.
    indices_x = ceil((midpoint_blh(:, 1) + epsilon - jdmin) / jdjg);
    indices_y = ceil((midpoint_blh(:, 2) + epsilon - wdmin) / wdjg);
    indices_z = ceil((midpoint_blh(:, 3) + epsilon - gdmin) / gdjg);
    
    % --- Boundary Clamping (Optional but Recommended) ---
    % The following commented-out code can be used to ensure that any
    % points calculated to be outside the grid bounds (e.g., due to
    % floating point errors or points outside the defined volume) are
    % clamped to the nearest valid grid cell. This prevents out-of-bounds
    % indexing errors.
    % indices_x(indices_x > Nx) = Nx;
    % indices_y(indices_y > Ny) = Ny;
    % indices_z(indices_z > Nz) = Nz;
    % indices_x(indices_x < 1) = 1;
    % indices_y(indices_y < 1) = 1;
    % indices_z(indices_z < 1) = 1;

    % --- Convert 3D indices to a single linear index ---
    % This manual calculation is faster than sub2ind.
    % It assumes a column-major ordering, which is standard in MATLAB.
    % The formula is: linear_index = (k-1)*(Nx*Ny) + (j-1)*Nx + i
    % where i, j, k are the indices along x, y, z dimensions respectively.
    grid_indices = (indices_z - 1) * (Nx * Ny) + (indices_y - 1) * Nx + indices_x;
end
