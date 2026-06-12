%此函数采用二分法+球面三角法来对单独的射线生成A
function [one_dimensional_index,sortedjd_ddzb,values,station,satellite,scope_blh] = Get_sta_A1_optimized(XD,YD,ZD,XS,YS,ZS,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length,nj,nw,ng,h_grid)
%%% 初始化输出
R = 6371000;
station = [XD,YD,ZD];
satellite = [XS,YS,ZS];
[one_dimensional_index, sortedjd_ddzb, values] = deal([]);

%%% 参数预处理
epsilon = 1e-10; max_iters = 100;
[lon1,lon2,lat1,lat2,is_valid] = get_first_and_last_point(jdmin_rad, jdmax_rad, wdmin_rad, wdmax_rad, gdmin, gdmax,epsilon,max_iters,XD,YD,ZD,XS,YS,ZS,jdjg_rad,wdjg_rad);
scope_blh = [lon1,lon2,lat1,lat2,gdmin,gdmax];
if ~is_valid, return; end

%%% 坐标转换预处理
[B_U, L_U, H_U] = XYZtoBLH_sphere(XD,YD,ZD,R);
[E,A] = Get_EA(XD,YD,ZD,XS,YS,ZS,R);
z = pi/2 - E;
sinz = sin(z); cosA = cos(A); sinA = sin(A);
sinBU = sin(B_U); cosBU = cos(B_U);
RH = R+H_U;
%%% 向量化经度面计算
jd_grid = lon1:jdjg_rad:lon2;
if isempty(jd_grid), jd_grid = deg2rad(lon1); end
delta_Ls = jd_grid - L_U;

% 预计算常用项
tan_half_B = tan((pi/2 - B_U)/2);

% 向量化计算
cos_terms = cos((A - delta_Ls)/2) ./ cos((A + delta_Ls)/2);
sin_terms = sin((A - delta_Ls)/2) ./ sin((A + delta_Ls)/2);

% 计算a和b
sum_ab = 2 * atan(cos_terms * tan_half_B);
diff_ab = 2 * atan(sin_terms * tan_half_B);

a = (sum_ab + diff_ab)/2;
b = (sum_ab - diff_ab)/2;

% 计算IPP参数
B_IPP_JD = pi/2 - a;
z_ = z - b;
H_IPP_JD = ((R+H_U)*sin(z)) ./ sin(z_) - R;

% 构建经度面结果
valid_jd = (H_IPP_JD >= gdmin) & (H_IPP_JD <= gdmax);
jdmjd = [B_IPP_JD(valid_jd);jd_grid(valid_jd);  H_IPP_JD(valid_jd)]';

%%% 向量化纬度面计算
wd_grid = lat1:wdjg_rad:lat2;
% % 处理单点的情况
% if isempty(wd_grid)
%     wd_grid = deg2rad(lat1);
% end
% % 使用数值方法
% B_IPP = wd_grid(:)';
% b = pi/2 - B_U;
% a = pi/2 - B_IPP;
% % 计算 B 的两个可能值
% sin_B = sin(b) ./ sin(a) .* sin(A);
% % B1 = real(asin(sin_B));
% sin_B = max(-1, min(1, sin_B));
% B1 = asin(sin_B);
% B2 = pi - B1;  % 第二个可能的 B 值
% % 将两组 B 值合并为一个 2xN 矩阵
% B_all = [B1; B2];
% % 预计算中间变量（扩展维度以同时处理两组解）
% A_plus_B = A + B_all;
% A_minus_B = A - B_all;
% a_plus_b = a + b;
% a_minus_b = a - b;
% % 计算 c 和 C（同时处理两组解）
% tan_c_over_2 = cos(A_plus_B ./ 2) ./ cos(A_minus_B ./ 2) .* tan(a_plus_b ./ 2);
% c = 2 * atan(tan_c_over_2);
% tan_C_over_2 = cos(a_minus_b ./ 2) ./ cos(a_plus_b ./ 2) .* cot(A_plus_B ./ 2);
% C = 2 * atan(tan_C_over_2);
% % 计算 z_ 和 H_（同时处理两组解）
% z_ = z - c;
% H_ = ((H_U + 6371000) .* sin(z)) ./ sin(z_) - 6371000;
% % 筛选有效解（同时对两组解进行筛选）
% valid = (H_ >= gdmin) & (H_ <= gdmax) &  abs(B1)~=pi/2;
% % 提取有效解并合并结果
% % wdmjd = [];
% % for k = 1:2
% %     current_valid = valid(k, :);
% %     if any(current_valid)
% %         jd = [rad2deg(wd_grid(current_valid)); ...
% %                  rad2deg(L_U + C(k, current_valid)); ...
% %                  H_(k, current_valid)]';
% %         wdmjd = [wdmjd; jd];
% %     end
% % end
% if any(valid(:))
%     % 扩展纬度网格以匹配 valid 的维度
%     wd_grid_expanded = repmat(wd_grid, 2, 1);
% 
%     % 一次性计算所有转换
%     valid_lat = wd_grid_expanded(valid);
%     valid_lon = L_U + C(valid);  % 先相加再转换
% 
%     % 一次性构建结果矩阵（减少rad2deg调用）
%     wdmjd = [valid_lat, ...           % 纬度
%              valid_lon, ...           % 经度
%              H_(valid)];                       % 高度
% else
%     wdmjd = [];
% end
% wdmjd = calc_lat_intersections_on_ray( ...
%     XD, YD, ZD, XS, YS, ZS, wd_grid, R, ...
%     jdmin_rad, jdmax_rad, wdmin_rad, wdmax_rad, gdmin, gdmax);
B0 = wd_grid;
tol_rhs   = 1e-13;
tol_alpha = 1e-12;
tol_h     = 1e-6;
tol_box   = 1e-10;
az_tol    = 1e-13;

RHsinz = RH * sinz;



if abs(sinA) < az_tol
    %==============================================================
    % 子午线快速分支：A≈0 或 A≈pi
    % 这是你当前问题最关键的退化情况。
    %
    % A≈0  ：alpha = B0 - B_U, lon = L_U
    % A≈pi ：alpha = B_U - B0, lon = L_U
    %==============================================================

    if cosA >= 0
        alpha = B0 - B_U;      % 正北
    else
        alpha = B_U - B0;      % 正南
    end

    % 极小负数修正
    % alpha(alpha < 0 & alpha > -tol_alpha) = 0;

    den = sin(z - alpha);
    H_  = RHsinz ./ den - R;

    L_ = L_U + zeros(size(B0));

    valid = alpha >= 0 & ...
            alpha < z - tol_alpha & ...
            den > 0 & ...
            H_ >= gdmin - tol_h & H_ <= gdmax + tol_h & ...
            L_ >= jdmin_rad - tol_box & L_ <= jdmax_rad + tol_box;

    if any(valid)
        wdmjd = [B0(valid).', L_(valid).', H_(valid).'];
    else
        wdmjd = [];
    end

else
    %==============================================================
    % 通用全向量化分支
    %
    % sin(B0) = sinBU*cos(alpha) + cosBU*sin(alpha)*cosA
    %
    % p*cos(alpha) + q*sin(alpha) = sin(B0)
    % M*cos(alpha - gamma) = sin(B0)
    %==============================================================

    p = sinBU;
    q = cosBU * cosA;

    M = hypot(p, q);

    % if M < eps
    %     % 极特殊情况：例如赤道附近严格东西向。
    %     % 射线可能沿某个纬度面走，不是普通离散交点。
    %     wdmjd = [];
    % else
        invM = 1 / M;

        rhs = sin(B0) * invM;

        valid_rhs = abs(rhs) <= 1 + tol_rhs;

        % 防止 acos 因舍入误差出现复数或 NaN
        rhs = min(1, max(-1, rhs));

        gamma = atan2(q, p);
        delta = acos(rhs);

        % 两个候选解，拼成 1 × 2N，避免 2 × N 的 B_all / alpha_all / H_ 等矩阵
        alpha = [gamma + delta, gamma - delta];
        Bv    = [B0, B0];

        valid_rhs2 = [valid_rhs, valid_rhs];

        % 极小负数修正
        % alpha(alpha < 0 & alpha > -tol_alpha) = 0;

        sin_alpha = sin(alpha);
        cos_alpha = cos(alpha);

        den = sin(z - alpha);
        H_  = RHsinz ./ den - R;

        C = atan2( ...
            sin_alpha .* sinA, ...
            cosBU .* cos_alpha - sinBU .* sin_alpha .* cosA );

        L_ = L_U + C;

        % 如果你的经度统一为 [-pi, pi]，可打开
        % L_ = mod(L_ + pi, 2*pi) - pi;

        valid = valid_rhs2 & ...
                alpha >= 0 & ...
                alpha < z - tol_alpha & ...
                den > 0 & ...
                H_ >= gdmin - tol_h & H_ <= gdmax + tol_h & ...
                Bv >= wdmin_rad - tol_box & Bv <= wdmax_rad + tol_box & ...
                L_ >= jdmin_rad - tol_box & L_ <= jdmax_rad + tol_box;

        if any(valid)
            wdmjd = [Bv(valid).', L_(valid).', H_(valid).'];
        else
            wdmjd = [];
        end
    % end
end
%%% 向量化高度面计算
% h_grid = gdmin:gdjg:gdmax;
if isempty(h_grid), h_grid = gdmin; end

z_gdm = asin((R+H_U)*sin(z)./(R+h_grid));
alpha = z - z_gdm;
B_IPP_GD = asin(cos(alpha)*sin(B_U) + sin(alpha)*cos(B_U)*cos(A));
L_IPP_GD = L_U + asin(sin(alpha).*sin(A)./cos(B_IPP_GD));

% 构建高度面结果
% valid_gd = (B_IPP_GD >= deg2rad(wdmin)) & (B_IPP_GD <= deg2rad(wdmax)) & ...
%            (L_IPP_GD >= deg2rad(jdmin)) & (L_IPP_GD <= deg2rad(jdmax));
% gdmjd = [rad2deg(B_IPP_GD(valid_gd));rad2deg(L_IPP_GD(valid_gd));  h_grid(valid_gd)]';
gdmjd = [B_IPP_GD(:), L_IPP_GD(:), h_grid(:)];   % n x 3

% 如果 jdmjd 和 wdmjd 也是 n×3 或者是 m×3 等，直接用 vertcat（等价于 [..;..]）
all_points = vertcat(jdmjd, wdmjd, gdmjd);      % 更语义化，也稍微更直观
tol = 1e-6;
valid_cols = (all_points(:,1) >= wdmin_rad - tol) & ...
             (all_points(:,1) <= wdmax_rad + tol) & ...
             (all_points(:,2) >= jdmin_rad - tol) & ...
             (all_points(:,2) <= jdmax_rad + tol) & ...
             (all_points(:,3) >= gdmin - tol) & ...
             (all_points(:,3) <= gdmax + tol);
valid_points = all_points(valid_cols,:);
if isempty(valid_points) || min(valid_points(:,3))> gdmin + tol  || max(valid_points(:,3))< gdmax - tol
    sortedjd_ddzb = [];
    return;
end
% 转换坐标
% [X,Y,Z] = BLHtoXYZ_sphere(valid_points(:,1),valid_points(:,2),valid_points(:,3),R);
% sortedjd = [X, Y, Z]';
% sortedjd_ddzb = valid_points';
[X,Y,Z] = BLHtoXYZ_sphere(valid_points(:,1), valid_points(:,2), valid_points(:,3), R);
sortedjd_ddzb = valid_points';
%%% 计算距离
% [d2, idx] = sort(d2);
vx = XS - XD; vy = YS - YD; vz = ZS - ZD;
vnorm = sqrt(vx*vx + vy*vy + vz*vz);

% 对每个交点 P，t = dot(P-D, v) / |v|   （这是沿射线的有向距离）
% t = ((sortedjd(1,:)-XD)*vx + (sortedjd(2,:)-YD)*vy + (sortedjd(3,:)-ZD)*vz) / vnorm;
t = ((X(:)' - XD) * vx + ...
     (Y(:)' - YD) * vy + ...
     (Z(:)' - ZD) * vz) / vnorm;
[t, idx] = sort(t);
values = diff(t);   % 直接就是每段长度
sortedjd_ddzb = sortedjd_ddzb(:, idx);

% if sortedjd_ddzb(3,1)> gdmin + 1e-10  || sortedjd_ddzb(3,end)< gdmax - 1e-10
%     sortedjd_ddzb = [];
%     return;
% end

%%% 计算一维索引
% if length(distances) < 2
%     return;
% end

% 计算两两交点的中点坐标
% midpoints = (sortedjd(:,1:end-1) + sortedjd(:,2:end))/2;
% % 批量转换中点坐标
% [lat_mid, lon_mid, alt_mid] = XYZtoBLH_sphere(midpoints(1,:), midpoints(2,:), midpoints(3,:),R);
% midpoint_blh = [rad2deg(lon_mid)', rad2deg(lat_mid)', alt_mid'];
% midpoint_blh = ((sortedjd_ddzb(:,1:end-1) + sortedjd_ddzb(:,2:end)) * 0.5);
% % RAD2DEG = 180/pi;
% % midpoint_blh([1,2],:) = midpoint_blh([2,1],:) * RAD2DEG;
% midpoint_blh([1,2],:) = midpoint_blh([2,1],:);
% midpoint_blh(2,:) = round(midpoint_blh(2,:),10);
% midpoint_blh = midpoint_blh';
% mid_lon = round(0.5 * (sortedjd_ddzb(2,1:end-1) + sortedjd_ddzb(2,2:end)), 13);
% mid_lat = round(0.5 * (sortedjd_ddzb(1,1:end-1) + sortedjd_ddzb(1,2:end)), 13);
% mid_h   = 0.5 * (sortedjd_ddzb(3,1:end-1) + sortedjd_ddzb(3,2:end));
% mid_lon = 0.5 * (sortedjd_ddzb(2,1:end-1) + sortedjd_ddzb(2,2:end));
% mid_lat = 0.5 * (sortedjd_ddzb(1,1:end-1) + sortedjd_ddzb(1,2:end));
% mid_h   = 0.5 * (sortedjd_ddzb(3,1:end-1) + sortedjd_ddzb(3,2:end));
% %% 索引计算优化
% epsilon = 1e-13;
% grid_indices = find_grid_indices_fast(mid_lon, mid_lat, mid_h, ...
%     jdmin_rad, wdmin_rad, gdmin, jdjg_rad, wdjg_rad, gdjg, epsilon, nj, nw);
% grid_indices = find_grid_indices(midpoint_blh, jdmin_rad, wdmin_rad, gdmin, jdjg_rad, wdjg_rad, gdjg, jdmax, wdmax, gdmax, epsilon,nj,nw,ng);
t_mid = 0.5 * (t(1:end-1) + t(2:end));

ux = vx / vnorm;
uy = vy / vnorm;
uz = vz / vnorm;

Xmid = XD + ux * t_mid;
Ymid = YD + uy * t_mid;
Zmid = ZD + uz * t_mid;

[lat_mid, lon_mid, mid_h] = XYZtoBLH_sphere(Xmid, Ymid, Zmid, R);

epsilon = 1e-13;
grid_indices = find_grid_indices_fast(lon_mid(:), lat_mid(:), mid_h(:), ...
    jdmin_rad, wdmin_rad, gdmin, jdjg_rad, wdjg_rad, gdjg, ...
    epsilon, nj, nw);
%% 稀疏矩阵优化
% total_length = round(((wdmax-wdmin)/wdjg) * ((jdmax-jdmin)/jdjg) * ((gdmax-gdmin)/gdjg));
% one_dimensional_index = sparse(1, grid_indices, values, 1, total_length);
one_dimensional_index = sparse(grid_indices, 1, values, total_length, 1);

