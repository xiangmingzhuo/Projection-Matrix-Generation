%% Helper Function: Calculate the intersection points of a ray with a grid's bounding box.
function [entry, exit, is_intersecting] = findGridIntersections(origin, direction, grid_min, grid_max)
    % This function implements the "slab method" for ray-AABB intersection.
    % It calculates the entry and exit points of a ray defined by an origin
    % and a direction vector as it passes through a 3D axis-aligned bounding box.
    %
    % Inputs:
    %   origin    - 1x3 vector, the origin point of the ray in Cartesian coordinates.
    %   direction - 1x3 vector, the direction vector of the ray (does not need to be normalized).
    %   grid_min  - 1x3 vector, the minimum corner of the bounding box.
    %   grid_max  - 1x3 vector, the maximum corner of the bounding box.
    %
    % Outputs:
    %   entry          - 1x3 vector, the coordinates where the ray enters the box.
    %   exit           - 1x3 vector, the coordinates where the ray leaves the box.
    %   is_intersecting- Boolean, true if the ray intersects the box, false otherwise.

    % Initialize the parametric distances (t) for entry and exit for each axis (x, y, z).
    % t_enter stores the near intersection distances for each slab.
    % t_exit stores the far intersection distances for each slab.
    t_enter = -inf(1, 3);
    t_exit = inf(1, 3);
    
    % Iterate over each dimension (x, y, z) to calculate intersection with the "slabs".
    for i = 1:3
        % Check if the ray is parallel to the current pair of slabs (faces).
        if abs(direction(i)) > eps
            % Calculate the parametric distances (t) to the two faces perpendicular to the i-th axis.
            t1 = (grid_min(i) - origin(i)) / direction(i);
            t2 = (grid_max(i) - origin(i)) / direction(i);
            
            % Determine the near and far intersection distances for this axis.
            t_enter(i) = min(t1, t2);
            t_exit(i) = max(t1, t2);
        % If the ray is parallel, it can only intersect if the origin is already within the slab.
        % This case is implicitly handled by the final check on t_entry and t_exit.
        end
    end
    
    % The ray enters the box at the maximum of the near distances (the last slab it enters).
    t_entry_final = max(t_enter);
    % The ray exits the box at the minimum of the far distances (the first slab it exits).
    t_exit_final = min(t_exit);
    
    % A valid intersection occurs only if the final entry point is before the final exit point.
    is_intersecting = (t_entry_final <= t_exit_final);
    
    if is_intersecting
        % Calculate the 3D coordinates of the entry and exit points.
        entry = origin + t_entry_final * direction;
        exit = origin + t_exit_final * direction;
    else
        % If there is no intersection, return NaNs for the coordinates.
        entry = NaN(1, 3);
        exit = NaN(1, 3);
    end
end
