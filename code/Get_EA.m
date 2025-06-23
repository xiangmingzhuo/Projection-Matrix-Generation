function [E, A] = Get_EA(sx, sy, sz, x, y, z,R)
    % 计算BLH坐标（只计算需要的B和L）
    [sb, sl] = XYZtoBLH_sphere(sx, sy, sz,R);
    
    % 预计算常用三角函数值
    sin_sb = sin(sb);
    cos_sb = cos(sb);
    sin_sl = sin(sl);
    cos_sl = cos(sl);
    
    % 计算坐标差（向量化）
    deta_xyz = [x - sx; y - sy; z - sz];
    
    % 直接计算NEU分量（避免完整矩阵乘法）
    NEU = zeros(3, 1);
    NEU(1) = -sin_sb*cos_sl*deta_xyz(1) - sin_sb*sin_sl*deta_xyz(2) + cos_sb*deta_xyz(3);
    NEU(2) = -sin_sl*deta_xyz(1) + cos_sl*deta_xyz(2);
    NEU(3) = cos_sb*cos_sl*deta_xyz(1) + cos_sb*sin_sl*deta_xyz(2) + sin_sb*deta_xyz(3);
    
    % 计算高度角E（使用atan2避免除零错误）
    E = atan2(NEU(3), sqrt(NEU(1)^2 + NEU(2)^2));
    
    % 计算方位角A（优化逻辑判断）
    A = atan2(NEU(2), NEU(1));
    if A < 0
        A = A + 2*pi;
    end
end