function [lon1, lon2, lat1, lat2,is_valid,B_U, L_U, H_U,A,z,R] = get_first_and_last_point2(jdmin, jdmax, wdmin, wdmax, gdmin, gdmax, XD, YD, ZD, XS, YS, ZS, jdjg, wdjg)
    % 初始化所有输出变量
    lon1 = []; lon2 = []; lat1 = []; lat2 = []; 
    is_valid = false;
    
    % 常量定义 - 避免在循环中重复计算
    R = 6371000;  % 地球半径
    
    % 预计算常用值 - 减少重复计算
    jdmin_rad = deg2rad(jdmin);
    jdmax_rad = deg2rad(jdmax);
    wdmin_rad = deg2rad(wdmin);
    wdmax_rad = deg2rad(wdmax);
    
    % 计算接收机位置和射线方向
    [B_U, L_U, H_U] = XYZtoBLH_sphere(XD, YD, ZD,R);
    [E, A] = Get_EA(XD, YD, ZD, XS, YS, ZS,R);
    
    % 预计算三角函数值 - 避免重复计算
    z = pi/2 - E;
    sin_z = sin(z);
    sin_B_U = sin(B_U);
    cos_B_U = cos(B_U);
    sin_A = sin(A);
    cos_A = cos(A);
    
    % 高度面交点计算 (向量化)
    heights = [gdmin, gdmax];
    sin_z_gdm = (R + H_U) * sin_z ./ (R + heights);
    
    % 检查是否有效高度
    if any(abs(sin_z_gdm) > 1)
        return; % 无效射线
    end
    
    z_gdm = asin(sin_z_gdm);
    alpha = z - z_gdm;
    
    % 预计算alpha的三角函数
    sin_alpha = sin(alpha);
    cos_alpha = cos(alpha);
    
    % 计算交点位置
    term1 = cos_alpha .* sin_B_U;
    term2 = sin_alpha .* cos_B_U .* cos_A;
    B_IPP_GD = asin(term1 + term2);
    
    % 计算经度增量
    sin_dL = sin_alpha .* sin_A;
    cos_B_IPP = cos(B_IPP_GD);
    
    % 检查有效性
    if any(abs(sin_dL) > abs(cos_B_IPP))
        return; % 无效射线
    end
    
    dL = asin(sin_dL ./ cos_B_IPP);
    L_IPP_GD = L_U + dL;
    
    % 组合结果
    gdmjd = [B_IPP_GD; L_IPP_GD; heights];
    
    % 提取首尾点
    gdm_min = gdmjd(:, 1);
    gdm_max = gdmjd(:, 2);
    
    % 合并边界检查 - 减少重复代码
    if ~(check_in_bounds(gdm_min(1), wdmin_rad, wdmax_rad) && ...
        check_in_bounds(gdm_min(2), jdmin_rad, jdmax_rad) && ...
        check_in_bounds(gdm_max(1), wdmin_rad, wdmax_rad) && ...
        check_in_bounds(gdm_max(2), jdmin_rad, jdmax_rad))
        return;
    end
    
    % 转换为角度
    lat11 = rad2deg(gdm_min(1));
    lon11 = rad2deg(gdm_min(2));
    lat22 = rad2deg(gdm_max(1));
    lon22 = rad2deg(gdm_max(2));
    
    % 预计算最小最大值
    min_lon = min(lon11, lon22);
    max_lon = max(lon11, lon22);
    min_lat = min(lat11, lat22);
    max_lat = max(lat11, lat22);
    
    % 最终确定边界
    lon1 = max(approximateNumberDown(min_lon, jdjg) - jdjg, jdmin);
    lon2 = min(approximateNumberUp(max_lon, jdjg) + jdjg, jdmax);
    lat1 = max(approximateNumberDown(min_lat, wdjg) - wdjg, wdmin);
    lat2 = min(approximateNumberUp(max_lat, wdjg) + wdjg, wdmax);
    is_valid = true;
end

