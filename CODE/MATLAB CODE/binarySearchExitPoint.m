function exit_point = binarySearchExitPoint(Bmin, Bmax, r0, d, t_first, t_far, epsilon, max_iters)
% BINARYSEARCHEXITPOINT Uses a binary search to find a ray's exit point from a volume.
%
% This function refines the exit point of a ray from a mesh or grid volume.
% It starts with a known entry point (t_first) and an initial far bound (t_far)
% and uses a binary search to converge on the precise location where the ray
% exits the volume, within a specified tolerance.
%
% INPUTS:
%   Bmin, Bmax - The minimum and maximum Cartesian coordinates defining the
%                bounding box of the volume (e.g., [Xmin, Ymin, Zmin]).
%   r0         - The ray's origin point [3x1].
%   d          - The ray's direction vector (must be normalized) [3x1].
%   t_first    - The ray parameter 't' for the known entry point.
%   t_far      - An initial upper bound for the exit point parameter 't'.
%   epsilon    - The desired tolerance for the search. The algorithm stops
%                when the search interval is smaller than this value.
%   max_iters  - The maximum number of iterations to perform.
%
% OUTPUT:
%   exit_point - The calculated 3D coordinates of the ray's exit point [3x1].

    % Initialize the search interval.
    % The left bound is the entry point. Ensure it's not negative.
    if t_first < 0
        t_left = 0;
    else
        t_left = t_first;
    end
    t_right = t_far;
    
    % Perform the binary search.
    for iter = 1:max_iters
        % Check if the interval is smaller than the tolerance.
        if (t_right - t_left) < epsilon
            break; % Convergence achieved.
        end
        
        % Calculate the midpoint of the interval.
        t_mid = (t_left + t_right) / 2;
        
        % Calculate the 3D position of the midpoint along the ray.
        p = r0 + t_mid * d;
        
        % Determine if the midpoint is inside the volume.
        if isInsideMesh(p, Bmin, Bmax)
            % If inside, the exit point must be further along the ray.
            % Move the left bound to the midpoint.
            t_left = t_mid;
        else
            % If outside, the exit point is between the left bound and the midpoint.
            % Move the right bound to the midpoint.
            t_right = t_mid;
        end
    end
    
    % The exit point is the last point known to be inside the volume.
    % This corresponds to the final left bound of the search interval.
    exit_point = r0 + t_left * d;
end
