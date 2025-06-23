function [lat, lon, h] = xyz2blh(x, y, z)
    % WGS84椭球体参数
    a = 6378137.0;
    f = 1/298.257223563;
    b = a*(1-f);
    e2 = 1 - (b/a)^2;
    
    p = sqrt(x.^2 + y.^2);
    lon = atan2(y, x);
    
    % 迭代计算纬度
    lat = atan2(z, p*(1-e2)); % 初始值
    for i = 1:10
        N = a / sqrt(1 - e2*sin(lat).^2);
        h = p./cos(lat) - N;
        lat_new = atan2(z, p*(1 - e2*N./(N + h)));
        if abs(lat_new - lat) < 1e-12
            break;
        end
        lat = lat_new;
    end
    
    lat = rad2deg(lat);
    lon = rad2deg(lon);
end