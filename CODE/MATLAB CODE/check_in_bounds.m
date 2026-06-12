function in_bounds = check_in_bounds(value, min_val, max_val)
% CHECK_IN_BOUNDS Determines if a scalar value falls within a specified range.
%
%   in_bounds = check_in_bounds(value, min_val, max_val) returns a logical
%   true if 'value' is greater than or equal to 'min_val' and less than or
%   equal to 'max_val'. Otherwise, it returns false.
%
%   Input Arguments:
%       value   - The scalar value to check.
%       min_val - The lower bound of the range.
%       max_val - The upper bound of the range.
%
%   Output Arguments:
%       in_bounds - A logical scalar (true or false).

    in_bounds = (value >= min_val-1e-6) && (value <= max_val+1e-6);
end
