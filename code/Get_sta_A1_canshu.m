%此函数用参数方程法来对单独的射线生成A
function [one_dimensional_index, sortedjd_ddzb, distances] = Get_sta_A1_canshu(XD,YD,ZD,XS,YS,ZS,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg)
% 初始化输出
[one_dimensional_index, sortedjd_ddzb, distances] = deal([]);
lon1 = jdmin; lon2 = jdmax;
lat1 = wdmin; lat2 = wdmax;
h1 = gdmin; h2 = gdmax;
%计算射线方向向量
orientVector = [XS,YS,ZS] - [XD,YD,ZD]; % 计算方向向量
% 预计算常用量
orient = orientVector(1,:);
tan_i = @(x) tan(x);
R = 6371000;
%% 向量化经度面计算
% 生成经度网格
lon_grid = deg2rad(lon1):deg2rad(jdjg):deg2rad(lon2);
if isempty(lon_grid), lon_grid = deg2rad(lon1); end  % 处理单点情况
% 向量化计算
t_lon = (YD .* cos(lon_grid) - XD .* sin(lon_grid)) ./ (orientVector(1,1) .* sin(lon_grid) - orientVector(1,2) .* cos(lon_grid));
t_lon = t_lon(:);
%% 向量化纬度面计算
% 生成纬度网格
lat_grid = deg2rad(lat1):deg2rad(wdjg):deg2rad(lat2);
if isempty(lat_grid), lat_grid = deg2rad(lat1); end % 处理单点情况
% 预计算公共项
A_term = orient(3)^2 - (orient(1)^2 + orient(2)^2)*tan(lat_grid).^2;
B_term = 2*(ZD*orient(3) - (XD*orient(1) + YD*orient(2)).*tan(lat_grid).^2);
C_term = ZD^2 - (XD^2 + YD^2).*tan(lat_grid).^2;
% 预分配内存
t_lat1 = zeros(size(lat_grid)) * NaN;
t_lat2 = zeros(size(lat_grid)) * NaN;
% 有效解计算
t_lat1 = (-B_term - sqrt(B_term.^2 - 4*A_term.*C_term)) ./ (2*A_term);
t_lat2 = (-B_term + sqrt(B_term.^2 - 4*A_term.*C_term)) ./ (2*A_term);
t_lat = [t_lat1(:); t_lat2(:)];
%% 向量化高度面计算
% 生成高度网格
h_grid = h1:gdjg:h2;
if isempty(h_grid), h_grid = h1; end % 处理单点情况
% 预计算系数
A = sum(orient.^2);
B = 2*(XD*orient(1) + YD*orient(2) + ZD*orient(3));
C_term = XD^2 + YD^2 + ZD^2 - (R + h_grid).^2;
% 预分配内存
t_alt1 = zeros(size(h_grid)) * NaN;
t_alt2 = zeros(size(h_grid)) * NaN;
% 有效解计算
t_alt1 = (-B - sqrt(B^2 - 4*A*C_term)) / (2*A);
t_alt2 = (-B + sqrt(B^2 - 4*A*C_term)) / (2*A);
% 整理结果
t_alt = [t_alt1(:);t_alt2(:)];
%% 合并筛选逻辑
% 合并所有t值并排序
t_all = [t_lon; t_lat; t_alt];
t_sorted = sort(t_all);
% 提前终止判断
if isempty(t_sorted)
    return;
end
%% 批量坐标转换
% 批量计算交点坐标
jd_coords = [XD,YD,ZD] + orient .* t_sorted;
jd_coords = jd_coords(imag(jd_coords(:,1)) == 0,:);
% 批量坐标转换
[lat, lon, alt] = XYZtoBLH_sphere(jd_coords(:,1), jd_coords(:,2), jd_coords(:,3),R);
jd_point = [rad2deg(lon), rad2deg(lat), alt];
%% 向量化有效性判断
valid_row = (jd_point(:,1) >= jdmin) & (jd_point(:,1) <= jdmax) & ...
            (jd_point(:,2) >= wdmin) & (jd_point(:,2) <= wdmax) & ...
            (jd_point(:,3) >= gdmin - 1e-4) & (jd_point(:,3) <= gdmax + 1e-4);
% 应用筛选
t_valid = t_sorted(valid_row);
if length(t_valid) < 2  % 需要至少两个点才能计算距离
    return;
end
%% 距离计算优化
distances = (diff(t_valid) * norm(orient))';
sortedjd = jd_coords(valid_row,:)';
sortedjd_ddzb = jd_point(valid_row,:)';
%% 中点计算优化
midpoint = (sortedjd(:,1:end-1) + sortedjd(:,2:end)) / 2;
% 批量转换中点坐标
[lat_mid, lon_mid, alt_mid] = XYZtoBLH_sphere(midpoint(1,:), midpoint(2,:), midpoint(3,:),R);
midpoint_blh = [rad2deg(lon_mid)', rad2deg(lat_mid)', alt_mid'];
%% 索引计算优化
epsilon = 1e-5;
grid_indices = find_grid_indices(midpoint_blh, jdmin, wdmin, gdmin, jdjg, wdjg, gdjg, jdmax, wdmax, gdmax, epsilon);
%% 稀疏矩阵优化
total_length = round(((wdmax-wdmin)/wdjg) * ((jdmax-jdmin)/jdjg) * ((gdmax-gdmin)/gdjg));
one_dimensional_index = sparse(1, grid_indices, distances, 1, total_length);
end