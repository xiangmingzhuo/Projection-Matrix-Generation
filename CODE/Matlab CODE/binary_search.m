% Generic Binary Search Function
function t = binary_search(f, t_low, t_high, epsilon, max_iters, find_entry)
% f: A function handle that returns true if a point at parameter t is inside the volume.
% t_low, t_high: The initial search interval.
% epsilon: The desired precision for the result.
% max_iters: The maximum number of iterations to prevent infinite loops.
% find_entry: A boolean flag. true to find the entry point, false to find the exit point.

% Ensure the lower bound is not negative, as t represents a distance along a ray.
if t_low < 0
    t_low = 0;
end

% Perform the binary search loop.
for i = 1:max_iters
    % Calculate the midpoint of the current interval.
    t_mid = (t_low + t_high) / 2;
    
    % Evaluate the function at the midpoint.
    inside = f(t_mid);

    % Adjust the search interval based on whether we're looking for an entry or exit point.
    if find_entry
        % To find the entry point, we are looking for the boundary from outside to inside.
        % If the midpoint is inside, the boundary must be to its left.
        if inside
            t_high = t_mid;
        else % If the midpoint is outside, the boundary must be to its right.
            t_low = t_mid;
        end
    else
        % To find the exit point, we are looking for the boundary from inside to outside.
        % If the midpoint is inside, the boundary must be to its right.
        if inside
            t_low = t_mid;
        else % If the midpoint is outside, the boundary must be to its left.
            t_high = t_mid;
        end
    end

    % Check if the interval is small enough to meet the precision requirement.
    if (t_high - t_low) < epsilon
        break; % Exit the loop.
    end
end

% Return the appropriate bound based on the search type.
% For entry, t_high converges to the first point that is 'inside'.
% For exit, t_low converges to the last point that is 'inside'.
if find_entry
    t = t_high;
else
    t = t_low;
end
end
