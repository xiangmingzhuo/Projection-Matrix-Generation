function [one_dimensional_index,sortedjd_ddzb,values] = Get_sta_A1_new_optimized1(XD,YD,ZD,XS,YS,ZS,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,total_length)
% 初始化输出
[one_dimensional_index, sortedjd_ddzb, values] = deal([]);

%计算首尾交点，根据首尾交点确定有效网格范围
[lon1,lon2,lat1,lat2,is_valid,B_U, L_U, H_U,A,z,R,gdmjd] = get_first_and_last_point2(jdmin, jdmax, wdmin, wdmax, gdmin, gdmax,XD,YD,ZD,XS,YS,ZS,jdjg,wdjg,gdjg);
% 如果不是有效射线，直接跳出函数
if ~is_valid, return; end 

%%% 优化经度面计算 - 保持原有逻辑不变
% 保持原有的 approximateNumberDown/Up 函数调用
jd_grid_deg = approximateNumberDown(lon1,jdjg):jdjg:approximateNumberUp(lon2,jdjg);
% 处理单点的情况
if isempty(jd_grid_deg)
    jd_grid_deg = lon1;
end

jd_grid = deg2rad(jd_grid_deg);
delta_Ls = jd_grid - L_U;

% 预计算常用项 - 但不改变计算顺序
tan_half_B = tan((pi/2 - B_U)/2);

% 向量化计算 - 保持原有公式
cos_terms = cos((A - delta_Ls)/2) ./ cos((A + delta_Ls)/2);
sin_terms = sin((A - delta_Ls)/2) ./ sin((A + delta_Ls)/2);

% 计算a和b - 保持原有计算顺序
sum_ab = 2 * atan(cos_terms * tan_half_B);
diff_ab = 2 * atan(sin_terms * tan_half_B);
a = (sum_ab + diff_ab)/2;
b = (sum_ab - diff_ab)/2;

% 计算IPP参数 - 保持原有公式
B_IPP_JD = pi/2 - a;
z_ = z - b;
H_IPP_JD = ((R+H_U)*sin(z)) ./ sin(z_) - R;

% 构建经度面结果 - 保持原有筛选逻辑
valid_jd = (H_IPP_JD >= gdmin) & (H_IPP_JD <= gdmax);
jdmjd = [rad2deg(B_IPP_JD(valid_jd)); jd_grid_deg(valid_jd); H_IPP_JD(valid_jd)]';

%%% 优化纬度面计算 - 保持原有逻辑
wd_grid_deg = approximateNumberDown(lat1, wdjg):wdjg:approximateNumberUp(lat2, wdjg);
% 处理单点的情况
if isempty(wd_grid_deg)
    wd_grid_deg = lat1;
end

wd_grid = deg2rad(wd_grid_deg);

% 使用数值方法 - 保持原有计算过程
B_IPP = wd_grid(:)';
b_val = pi/2 - B_U;
a_vals = pi/2 - B_IPP;

% 计算 B 的两个可能值 - 保持原有公式
sin_B = sin(b_val) ./ sin(a_vals) .* sin(A);
B1 = asin(sin_B);
B2 = pi - B1;

% 将两组 B 值合并 - 保持原有数据结构
B_all = [B1; B2];

% 预计算中间变量 - 保持原有计算顺序
A_plus_B = A + B_all;
A_minus_B = A - B_all;
a_plus_b = a_vals + b_val;
a_minus_b = a_vals - b_val;

% 计算 c 和 C - 保持原有公式
tan_c_over_2 = cos(A_plus_B ./ 2) ./ cos(A_minus_B ./ 2) .* tan(a_plus_b ./ 2);
c = 2 * atan(tan_c_over_2);
tan_C_over_2 = cos(a_minus_b ./ 2) ./ cos(a_plus_b ./ 2) .* cot(A_plus_B ./ 2);
C = 2 * atan(tan_C_over_2);

% 计算 z_ 和 H_ - 保持原有公式
z_ = z - c;
H_ = ((H_U + 6371000) .* sin(z)) ./ sin(z_) - 6371000;

% 筛选有效解 - 保持原有筛选条件
valid = (H_ >= gdmin) & (H_ <= gdmax) & (imag(H_) == 0);

% 提取有效解并合并结果 - 保持原有循环逻辑
wdmjd = [];
for k = 1:2
    current_valid = valid(k, :);
    if any(current_valid)
        jd = [rad2deg(wd_grid_deg(current_valid)); ...
                 rad2deg(L_U + C(k, current_valid)); ...
                 H_(k, current_valid)]';
        wdmjd = [wdmjd; jd];
    end
end

%%% 高度面计算 - 完全保持原有代码
h_grid = gdmin:gdjg:gdmax;
z_gdm = asin((R+H_U)*sin(z)./(R+h_grid));
alpha = z - z_gdm;
B_IPP_GD = asin(cos(alpha)*sin(B_U) + sin(alpha)*cos(B_U)*cos(A));
L_IPP_GD = L_U + asin(sin(alpha).*sin(A)./cos(B_IPP_GD));
gdmjd = [rad2deg(B_IPP_GD);rad2deg(L_IPP_GD); h_grid]';

%%% 合并所有结果点 - 保持原有逻辑
all_points = [jdmjd; wdmjd; gdmjd];

% 筛选有效交点 - 保持原有条件
valid_cols = (all_points(:,1) >= wdmin) & (all_points(:,1) <= wdmax) & ...
             (all_points(:,2) >= jdmin) & (all_points(:,2) <= jdmax) & ...
             (all_points(:,3) >= gdmin) & (all_points(:,3) <= gdmax);
valid_points = all_points(valid_cols,:);

% 转换坐标 - 保持原有函数调用
[X,Y,Z] = BLHtoXYZ_sphere(deg2rad(valid_points(:,1)),deg2rad(valid_points(:,2)),valid_points(:,3),R);
sortedjd = [X, Y, Z]';
sortedjd_ddzb = valid_points';

%%% 计算距离 - 保持原有公式
distances = sqrt((sortedjd(1,:) - XD).^2 + (sortedjd(2,:) - YD).^2 + (sortedjd(3,:) - ZD).^2);

% 对距离进行排序 - 保持原有逻辑
[distances, idx] = sort(distances);
sortedjd = sortedjd(:, idx);
sortedjd_ddzb = sortedjd_ddzb(:, idx);

%%% 计算一维索引 - 保持原有逻辑
if length(distances) < 2
    return;
end

% 计算两两交点的中点坐标 - 保持原有计算
midpoints = (sortedjd(:,1:end-1) + sortedjd(:,2:end))/2;

% 批量转换中点坐标 - 保持原有函数调用
[lat_mid, lon_mid, alt_mid] = XYZtoBLH_sphere(midpoints(1,:), midpoints(2,:), midpoints(3,:),R);
midpoint_blh = [rad2deg(lon_mid)', rad2deg(lat_mid)', alt_mid'];

% 计算两两交点之间的截距值 - 保持原有公式
values = distances(2:end) - distances(1:end-1);

%% 索引计算 - 保持原有函数调用
epsilon = 1e-15;
grid_indices = find_grid_indices(midpoint_blh, jdmin, wdmin, gdmin, jdjg, wdjg, gdjg, jdmax, wdmax, gdmax, epsilon);

%% 稀疏矩阵 - 保持原有逻辑
one_dimensional_index = sparse(1, grid_indices, values, 1, total_length);