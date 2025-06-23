function [x, y, z] = blh2xyz(lat, lon, h)
    % WGS84椭球体参数
    a = 6378137.0;       % 长半轴
    f = 1/298.257223563; % 扁率
    e2 = 2*f - f^2;      % 第一偏心率的平方
    
    lat_rad = deg2rad(lat);
    lon_rad = deg2rad(lon);
    
    N = a / sqrt(1 - e2*sin(lat_rad)^2);
    x = (N + h) * cos(lat_rad) * cos(lon_rad);
    y = (N + h) * cos(lat_rad) * sin(lon_rad);
    z = (N*(1-e2) + h) * sin(lat_rad);
end