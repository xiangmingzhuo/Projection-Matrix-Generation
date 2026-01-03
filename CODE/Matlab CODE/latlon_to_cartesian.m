% Convert latitude and longitude to Cartesian coordinates (simplified version)
function [x, y] = latlon_to_cartesian(lat, lon, h)
    % Convert latitude and longitude from degrees to radians
    lat_rad = deg2rad(lat);
    lon_rad = deg2rad(lon);
    
    % Calculate the x and y coordinates
    x = h * cos(lat_rad) * cos(lon_rad);
    y = h * cos(lat_rad) * sin(lon_rad);
end
