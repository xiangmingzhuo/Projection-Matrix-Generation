% Rounds a number down to the nearest multiple of an interval.
function approx = approximateNumberDown(num, interval)
    approx = floor(num / interval) * interval;
end
