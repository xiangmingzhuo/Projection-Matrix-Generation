% m向上取整到n的最大倍数
function approx = approximateNumberUp(num, interval)
    approx = ceil(num / interval) * interval;
end