% 辅助函数：计算网格边界
function [Bmin, Bmax] = compute_grid_bounds(lon_min, lon_max, lat_min, lat_max, h_min, h_max)
% 只需要计算边界，不需要所有顶点
R = 6371000;
h_min = R + h_min;
h_max = R + h_max;

% 计算经纬度对应的笛卡尔坐标
[Bmin(1), Bmin(2)] = latlon_to_cartesian(lat_min, lon_min, h_min);
[Bmax(1), Bmax(2)] = latlon_to_cartesian(lat_max, lon_max, h_max);

% 处理高度
Bmin(3) = h_min;
Bmax(3) = h_max;
end