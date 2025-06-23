% m向下取整到n的最大倍数
function approx = approximateNumberDown(num, interval)
    approx = floor(num / interval) * interval;
end