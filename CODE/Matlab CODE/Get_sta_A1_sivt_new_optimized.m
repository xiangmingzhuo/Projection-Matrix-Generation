%此函数用球面三角+有效分裂面方案二来对单独的射线生成A
function [one_dimensional_index,sortedjd_ddzb,values,jdmjd,wdmjd,gdmjd] =Get_sta_A1_sivt_new_optimized(XD,YD,ZD,XS,YS,ZS,jd_params,wd_params,gd_params,total_length,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad)
% 初始化输出
[one_dimensional_index, sortedjd_ddzb, values] = deal([]);
gdmin = gd_params(1); gdmax = gd_params(2); gdjg = gd_params(3);
%计算首尾交点，根据首尾交点确定有效网格范围
% [lon1,lon2,lat1,lat2,is_valid,B_U, L_U, H_U,A,z,R,gdmjd] = get_first_and_last_point2(jdmin, jdmax, wdmin, wdmax, gdmin, gdmax,XD,YD,ZD,XS,YS,ZS,jdjg,wdjg,gdjg);
[lon1,lon2,lat1,lat2,is_valid,B_U, L_U, H_U,A,z,R] = get_first_and_last_point2(jdmin_rad, jdmax_rad, wdmin_rad, wdmax_rad, gdmin, gdmax,XD,YD,ZD,XS,YS,ZS,jdjg_rad,wdjg_rad);

% 如果不是有效射线，直接跳出函数
if ~is_valid, return; end 
%%% 向量化经度面计算
% jd_grid = deg2rad(approximateNumberDown(lon1,jdjg)):deg2rad(jdjg):deg2rad(approximateNumberUp(lon2,jdjg));
jd_grid = lon1:jdjg_rad:lon2;
% 处理单点的情况
if isempty(jd_grid), jd_grid = lon1; end
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
% jdmjd = [rad2deg(B_IPP_JD(valid_jd)); jd_grid(valid_jd); H_IPP_JD(valid_jd)]';
jdmjd = [B_IPP_JD(valid_jd);jd_grid(valid_jd);  H_IPP_JD(valid_jd)]';

%%% 向量化纬度面计算
% wd_grid = deg2rad(approximateNumberDown(lat1, wdjg)):deg2rad(wdjg):deg2rad(approximateNumberUp(lat2, wdjg));
wd_grid = lat1:wdjg_rad:lat2;
% % 处理单点的情况
% if isempty(wd_grid)
%     wd_grid = deg2rad(lat1);
% end
% 使用数值方法
B_IPP = wd_grid(:)';
b = pi/2 - B_U;
a = pi/2 - B_IPP;
% 计算 B 的两个可能值
sin_B = sin(b) ./ sin(a) .* sin(A);
B1 = asin(sin_B);
B2 = pi - B1;  % 第二个可能的 B 值
% 将两组 B 值合并为一个 2xN 矩阵
B_all = [B1; B2];
% 预计算中间变量（扩展维度以同时处理两组解）
A_plus_B = A + B_all;
A_minus_B = A - B_all;
a_plus_b = a + b;
a_minus_b = a - b;
% 计算 c 和 C（同时处理两组解）
tan_c_over_2 = cos(A_plus_B ./ 2) ./ cos(A_minus_B ./ 2) .* tan(a_plus_b ./ 2);
c = 2 * atan(tan_c_over_2);
tan_C_over_2 = cos(a_minus_b ./ 2) ./ cos(a_plus_b ./ 2) .* cot(A_plus_B ./ 2);
C = 2 * atan(tan_C_over_2);
% 计算 z_ 和 H_（同时处理两组解）
z_ = z - c;
H_ = ((H_U + 6371000) .* sin(z)) ./ sin(z_) - 6371000;
% 筛选有效解（同时对两组解进行筛选）
valid = (H_ >= gdmin) & (H_ <= gdmax) & (imag(H_) == 0);
% 提取有效解并合并结果
% wdmjd = [];
% for k = 1:2
%     current_valid = valid(k, :);
%     if any(current_valid)
%         jd = [rad2deg(wd_grid(current_valid)); ...
%                  rad2deg(L_U + C(k, current_valid)); ...
%                  H_(k, current_valid)]';
%         wdmjd = [wdmjd; jd];
%     end
% end
% 一次性提取所有有效解
% 一次性构建所有有效解
if any(valid(:))
    % 扩展纬度网格以匹配 valid 的维度
    wd_grid_expanded = repmat(wd_grid, 2, 1);
    
    % 一次性计算所有转换
    valid_lat = wd_grid_expanded(valid);
    valid_lon = L_U + C(valid);  % 先相加再转换
    
    % 一次性构建结果矩阵（减少rad2deg调用）
    wdmjd = [valid_lat, ...           % 纬度
             valid_lon, ...           % 经度
             H_(valid)];                       % 高度
else
    wdmjd = [];
end

%%% 向量化高度面计算
h_grid = gdmin:gdjg:gdmax;
% % 处理单点的情况
% if isempty(h_grid), h_grid = gdmin; end
z_gdm = asin((R+H_U)*sin(z)./(R+h_grid));
alpha = z - z_gdm;
B_IPP_GD = asin(cos(alpha)*sin(B_U) + sin(alpha)*cos(B_U)*cos(A));
L_IPP_GD = L_U + asin(sin(alpha).*sin(A)./cos(B_IPP_GD));
% 构建高度面结果
% valid_gd = (B_IPP_GD >= deg2rad(wdmin)) & (B_IPP_GD <= deg2rad(wdmax)) & ...
%            (L_IPP_GD >= deg2rad(jdmin)) & (L_IPP_GD <= deg2rad(jdmax));
% gdmjd = [rad2deg(B_IPP_GD(valid_gd));rad2deg(L_IPP_GD(valid_gd));  h_grid(valid_gd)]';
gdmjd = [B_IPP_GD;L_IPP_GD;  h_grid]';
% %%% 合并所有结果点
all_points = [jdmjd;wdmjd;gdmjd];
% valid_points = [jdmjd;wdmjd;gdmjd];
% 筛选有效交点
valid_cols = (all_points(:,1) >= wdmin_rad) & (all_points(:,1) <= wdmax_rad) & ...
             (all_points(:,2) >= jdmin_rad) & (all_points(:,2) <= jdmax_rad) & ...
             (all_points(:,3) >= gdmin) & (all_points(:,3) <= gdmax);
valid_points = all_points(valid_cols,:);

% 转换坐标
[X,Y,Z] = BLHtoXYZ_sphere(valid_points(:,1),valid_points(:,2),valid_points(:,3),R);
sortedjd = [X, Y, Z]';
sortedjd_ddzb = valid_points';

%%% 计算距离
distances = sqrt((sortedjd(1,:) - XD).^2 + (sortedjd(2,:) - YD).^2 + (sortedjd(3,:) - ZD).^2);
% 对距离进行排序
[distances, idx] = sort(distances);
values = distances(2:end) - distances(1:end-1);
% 重新排列 new_jd 矩阵
sortedjd = sortedjd(:, idx);
sortedjd_ddzb = sortedjd_ddzb(:, idx);
sortedjd_ddzb = round(sortedjd_ddzb,13);
sortedjd_ddzb([1, 2], :) = rad2deg(sortedjd_ddzb([2, 1], :));
%% SIVT caculate grid index
[voxel_indices,grid_indices] = SIVT_voxel_tracing(sortedjd_ddzb, jd_params, wd_params, gd_params);

%% 稀疏矩阵优化
% total_length = round(((wdmax-wdmin)/wdjg) * ((jdmax-jdmin)/jdjg) * ((gdmax-gdmin)/gdjg));
one_dimensional_index = sparse(1, grid_indices, values, 1, total_length);
