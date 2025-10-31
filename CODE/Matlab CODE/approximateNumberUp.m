function approx = approximateNumberUp(num, interval)
% APPROXIMATENUMBERUP Rounds a number up to the nearest multiple of an interval.
%
%   approx = approximateNumberUp(num, interval) takes a numerical input 'num'
%   and rounds it up to the largest multiple of 'interval' that is greater
%   than or equal to 'num'.
%
%   Example:
%       approximateNumberUp(23, 10) returns 30
%       approximateNumberUp(20, 10) returns 20
%       approximateNumberUp(3.7, 0.5) returns 4.0

    approx = ceil(num / interval) * interval;
end
