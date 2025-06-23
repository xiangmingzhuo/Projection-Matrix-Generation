% 经纬度转笛卡尔坐标 (简化版)
function [x, y] = latlon_to_cartesian(lat, lon, h)
lat_rad = deg2rad(lat);
lon_rad = deg2rad(lon);
x = h * cos(lat_rad) * cos(lon_rad);
y = h * cos(lat_rad) * sin(lon_rad);
end