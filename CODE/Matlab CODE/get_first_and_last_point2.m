function [lon1, lon2, lat1, lat2, is_valid, B_U, L_U, H_U, A, z, R] = get_first_and_last_point2(jdmin_rad, jdmax_rad, wdmin_rad, wdmax_rad, gdmin, gdmax, XD, YD, ZD, XS, YS, ZS, jdjg_rad,wdjg_rad)
    % Initialize all output variables
    lon1 = []; lon2 = []; lat1 = []; lat2 = []; 
    is_valid = false;
    
    % Constants definition - avoid recalculating in loops
    R = 6371000;  % Earth's radius
    
    % Pre-calculate common values - reduce redundant calculations
    % jdmin_rad = deg2rad(jdmin);
    % jdmax_rad = deg2rad(jdmax);
    % wdmin_rad = deg2rad(wdmin);
    % wdmax_rad = deg2rad(wdmax);
    
    % Calculate receiver position and ray direction
    [B_U, L_U, H_U] = XYZtoBLH_sphere(XD, YD, ZD, R);
    [E, A] = Get_EA(XD, YD, ZD, XS, YS, ZS, R);
    
    % Pre-calculate trigonometric values - avoid redundant calculations
    z = pi/2 - E;
    sin_z = sin(z);
    sin_B_U = sin(B_U);
    cos_B_U = cos(B_U);
    sin_A = sin(A);
    cos_A = cos(A);
    
    % Intersection points with height planes (vectorized)
    heights = [gdmin, gdmax];
    sin_z_gdm = (R + H_U) * sin_z ./ (R + heights);
    
    % Check for valid heights
    if any(abs(sin_z_gdm) > 1)
        return; % Invalid ray
    end
    
    z_gdm = asin(sin_z_gdm);
    alpha = z - z_gdm;
    
    % Pre-calculate trigonometric functions of alpha
    sin_alpha = sin(alpha);
    cos_alpha = cos(alpha);
    
    % Calculate intersection point positions
    term1 = cos_alpha .* sin_B_U;
    term2 = sin_alpha .* cos_B_U .* cos_A;
    B_IPP_GD = asin(term1 + term2);
    
    % Calculate longitude increment
    sin_dL = sin_alpha .* sin_A;
    cos_B_IPP = cos(B_IPP_GD);
    
    % Check validity
    if any(abs(sin_dL) > abs(cos_B_IPP))
        return; % Invalid ray
    end
    
    dL = asin(sin_dL ./ cos_B_IPP);
    L_IPP_GD = L_U + dL;
    
    % Combine results
    gdmjd = [B_IPP_GD; L_IPP_GD; heights];
    
    % Extract first and last points
    gdm_min = gdmjd(:, 1);
    gdm_max = gdmjd(:, 2);
    
    % Consolidated boundary check - reduce repetitive code
    if ~(check_in_bounds(gdm_min(1), wdmin_rad, wdmax_rad) && ...
        check_in_bounds(gdm_min(2), jdmin_rad, jdmax_rad) && ...
        check_in_bounds(gdm_max(1), wdmin_rad, wdmax_rad) && ...
        check_in_bounds(gdm_max(2), jdmin_rad, jdmax_rad))
        return;
    end
    
    % Convert to degrees
    lat11 = gdm_min(1);
    lon11 = gdm_min(2);
    lat22 = gdm_max(1);
    lon22 = gdm_max(2);
    
    % Pre-calculate min/max values
    min_lon = min(lon11, lon22);
    max_lon = max(lon11, lon22);
    min_lat = min(lat11, lat22);
    max_lat = max(lat11, lat22);
    
    % Finalize boundaries
    % lon1 = max(approximateNumberDown(min_lon, jdjg) - jdjg, jdmin);
    % lon2 = min(approximateNumberUp(max_lon, jdjg) + jdjg, jdmax);
    % lat1 = max(approximateNumberDown(min_lat, wdjg) - wdjg, wdmin);
    % lat2 = min(approximateNumberUp(max_lat, wdjg) + wdjg, wdmax);
    lon1 = max(approximateNumberDown(min_lon, jdjg_rad) , jdmin_rad);
    lon2 = min(approximateNumberUp(max_lon, jdjg_rad) , jdmax_rad);
    lat1 = max(approximateNumberDown(min_lat, wdjg_rad) , wdmin_rad);
    lat2 = min(approximateNumberUp(max_lat, wdjg_rad) , wdmax_rad);    
    is_valid = true;
end