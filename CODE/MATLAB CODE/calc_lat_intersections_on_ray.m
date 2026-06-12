function [wdmjd, t_wd] = calc_lat_intersections_on_ray( ...
    XD, YD, ZD, XS, YS, ZS, wd_grid, R, ...
    jdmin_rad, jdmax_rad, wdmin_rad, wdmax_rad, gdmin, gdmax)

    % 输出初始化
    wdmjd = [];
    t_wd = [];

    % 射线方向
    vx = XS - XD;
    vy = YS - YD;
    vz = ZS - ZD;

    vv = vx*vx + vy*vy + vz*vz;
    if vv <= 0
        return;
    end

    vnorm = sqrt(vv);

    % 容差设置
    tolB_pre  = 1e-12;   % 纬度预筛选容差，rad
    tolT      = 1e-12;   % t 在线段 [0,1] 附近的容差
    tolLat    = 1e-9;    % 纬度残差检查容差，rad
    tolRange  = 1e-10;   % 经纬度范围容差，rad
    tolH      = 1e-5;    % 高度范围容差，m
    tolDist   = 1e-4;    % 重复交点距离容差，m

    % 纬度网格列向量
    B0 = wd_grid(:);

    % 预筛选：去掉 NaN、Inf、范围外纬度、接近极点的纬度
    % 接近 ±90° 时 tan(B) 会数值爆炸，若你的区域包含极点，需要单独处理
    valid_B = isfinite(B0) & ...
              B0 >= wdmin_rad - tolB_pre & ...
              B0 <= wdmax_rad + tolB_pre & ...
              abs(cos(B0)) > 1e-12;

    B0 = B0(valid_B);

    if isempty(B0)
        return;
    end

    nB = numel(B0);

    % 预计算公共项
    tanB = tan(B0);
    tanB2 = tanB .* tanB;

    xyV2 = vx*vx + vy*vy;
    xyDV = XD*vx + YD*vy;
    xyD2 = XD*XD + YD*YD;

    % 固定纬度面：
    % Z^2 - (X^2 + Y^2) * tan(B0)^2 = 0
    %
    % 射线：
    % P(lambda) = D + lambda * V, lambda in [0,1]
    %
    % 得到：
    % A*lambda^2 + B*lambda + C = 0

    Acoef = vz*vz - xyV2 .* tanB2;
    Bcoef = 2 * (ZD*vz - xyDV .* tanB2);
    Ccoef = ZD*ZD - xyD2 .* tanB2;

    disc = Bcoef.*Bcoef - 4 .* Acoef .* Ccoef;

    % 数值容差
    tolA = 1e-14 * max(vv, 1);
    scale_disc = Bcoef.*Bcoef + 4 * abs(Acoef .* Ccoef) + 1;
    tolDisc = 1e-12 * scale_disc;

    % 直接预分配两个解的位置
    t_all = nan(2*nB, 1);
    B_all = [B0; B0];

    % 二次方程有效解
    valid_quad = abs(Acoef) > tolA & disc >= -tolDisc;

    if any(valid_quad)
        iq = find(valid_quad);

        disc_q = disc(iq);
        disc_q = max(disc_q, 0);

        sqrt_disc = sqrt(disc_q);

        Aq = Acoef(iq);
        Bq = Bcoef(iq);

        t_all(iq) = (-Bq - sqrt_disc) ./ (2 * Aq);
        t_all(nB + iq) = (-Bq + sqrt_disc) ./ (2 * Aq);
    end

    % 退化为一次方程
    valid_linear = abs(Acoef) <= tolA & abs(Bcoef) > tolA;

    if any(valid_linear)
        il = find(valid_linear);
        t_all(il) = -Ccoef(il) ./ Bcoef(il);
        % 第二个解保持 NaN
    end

    % t 范围过滤：只保留线段内部或极小越界的点
    valid_t = isfinite(t_all) & ...
              t_all >= -tolT & ...
              t_all <= 1 + tolT;

    if ~any(valid_t)
        return;
    end

    t_all = t_all(valid_t);
    B_all = B_all(valid_t);

    % 极小越界值压回端点
    t_all(t_all < 0) = 0;
    t_all(t_all > 1) = 1;

    % ECEF 坐标，严格在射线上
    X = XD + t_all * vx;
    Y = YD + t_all * vy;
    Z = ZD + t_all * vz;

    % 球模型 BLH，直接内联计算，避免调用 XYZtoBLH_sphere
    rho = sqrt(X.*X + Y.*Y);
    lat = atan2(Z, rho);
    h = sqrt(rho.*rho + Z.*Z) - R;

    % 先过滤纬度和高度，减少后续 lon 计算量
    valid_basic = ...
        abs(lat - B_all) <= tolLat & ...
        lat >= wdmin_rad - tolRange & ...
        lat <= wdmax_rad + tolRange & ...
        h   >= gdmin     - tolH & ...
        h   <= gdmax     + tolH;

    if ~any(valid_basic)
        return;
    end

    t_all = t_all(valid_basic);
    X = X(valid_basic);
    Y = Y(valid_basic);
    lat = lat(valid_basic);
    h = h(valid_basic);

    % 只有通过前面过滤的点才计算经度
    lon = atan2(Y, X);

    valid_lon = lon >= jdmin_rad - tolRange & ...
                lon <= jdmax_rad + tolRange;

    if ~any(valid_lon)
        return;
    end

    t_all = t_all(valid_lon);
    lat = lat(valid_lon);
    lon = lon(valid_lon);
    h = h(valid_lon);

    % 按射线位置排序
    [t_all, idx] = sort(t_all);
    lat = lat(idx);
    lon = lon(idx);
    h = h(idx);

    % 去重
    if numel(t_all) > 1
        keep = [true; diff(t_all) * vnorm > tolDist];
        t_all = t_all(keep);
        lat = lat(keep);
        lon = lon(keep);
        h = h(keep);
    end

    % 输出
    wdmjd = [lat, lon, h];

    % 第二输出：沿射线的有向距离，单位 m
    % 如果主函数后面愿意改，可以直接用这个，避免再 BLHtoXYZ_sphere 回算 t
    t_wd = t_all * vnorm;
end
