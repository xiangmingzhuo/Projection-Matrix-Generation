function intersection_point = binarySearchEntryPoint(Bmin, Bmax, r0, d, t_near, t_far, epsilon, max_iters)
% BINARYSEARCHENTRYPOINT Uses a binary search to find a ray's entry point into a volume.
%
% This function refines the entry point of a ray into a mesh or grid volume.
% It starts with an initial interval [t_near, t_far] where the ray is known
% to be outside at t_near and inside at t_far, and it narrows this interval
% to find the precise boundary.
%
% Inputs:
%   Bmin, Bmax - The minimum and maximum Cartesian coordinates of the bounding box [3x1].
%   r0         - The origin of the ray [3x1].
%   d          - The direction vector of the ray (must be normalized) [3x1].
%   t_near     - The lower bound of the initial search interval (ray parameter t).
%   t_far      - The upper bound of the initial search interval (ray parameter t).
%   epsilon    - The desired tolerance for the interval size. The search stops
%                when (t_right - t_left) < epsilon.
%   max_iters  - The maximum number of iterations to perform.
%
% Outputs:
%   intersection_point - The calculated 3D Cartesian coordinate of the entry point [3x1].

    % Ensure the search interval does not start behind the ray's origin.
    if t_near < 0
        t_left = 0;
    else
        t_left = t_near;
    end
    t_right = t_far;
    
    % Perform the binary search.
    for iter = 1:max_iters
        % Check if the interval is smaller than the tolerance.
        if (t_right - t_left) < epsilon
            break;
        end
        
        % Calculate the midpoint of the current interval.
        t_mid = (t_left + t_right) / 2;
        
        % Calculate the 3D position of the midpoint along the ray.
        p = r0 + t_mid * d;
        
        % Determine if the midpoint is inside the volume.
        if isInsideMesh(p, Bmin, Bmax)
            % If the midpoint is inside, the entry point must be
            % between the left bound and the midpoint.
            % Move the right bound to the midpoint.
            t_right = t_mid;
        else
            % If the midpoint is outside, the entry point must be
            % further along the ray.
            % Move the left bound to the midpoint.
            t_left = t_mid;
        end
    end
    
    % The entry point is the first point known to be inside the volume.
    % This corresponds to the final right bound of the search interval.
    intersection_point = r0 + t_right * d;
end
