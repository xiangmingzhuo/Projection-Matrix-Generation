function [line_res, lambda, lat_res] = check_lat_point(P, XD,YD,ZD, XS,YS,ZS, B0)
% 检查某个候选点 P 是否：
% 1. 在 XD,YD,ZD -> XS,YS,ZS 这条射线上
% 2. 在目标纬度面 B0 上
%
% 输入：
% P        : 候选交点的 ECEF 坐标，必须是 [X;Y;Z]，单位 m
% XD,YD,ZD : 射线起点 ECEF
% XS,YS,ZS : 射线终点 ECEF
% B0       : 目标纬度，单位 rad
%
% 输出：
% line_res : 点 P 到射线所在直线的垂距，单位 m，越接近 0 越好
% lambda   : P 在 D->S 线段上的比例参数，0~1 表示在线段内部
% lat_res  : P 的纬度与目标纬度 B0 的差，单位 rad，越接近 0 越好

    D = [XD; YD; ZD];
    S = [XS; YS; ZS];
    V = S - D;

    P = P(:);

    % 点到直线的垂距
    line_res = norm(cross(P - D, V)) / norm(V);

    % 线段参数
    % P = D + lambda * V
    lambda = dot(P - D, V) / dot(V, V);

    % 球模型下的地心纬度
    lat = atan2(P(3), hypot(P(1), P(2)));

    % 纬度残差
    lat_res = lat - B0;

    fprintf('line_res = %.12e m\n', line_res);
    fprintf('lambda   = %.12f\n', lambda);
    fprintf('lat_res  = %.12e rad, %.12e deg\n', ...
        lat_res, rad2deg(lat_res));
end
