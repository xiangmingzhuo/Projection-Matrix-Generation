%此函数用方程求解的方法来对单独的射线生成A
function [one_dimensional_index,sortedjd_ddzb,values] =Get_sta_A1_initial(XD,YD,ZD,XS,YS,ZS,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg)
%%% 初始化输出
[one_dimensional_index, sortedjd_ddzb, values] = deal([]);
% 预存储网格点
lon1 = jdmin ; lon2 = jdmax; lat1 = wdmin; lat2 = wdmax;
%地面坐标、卫星坐标、经度范围及间隔(以度为单位)、纬度范围及间隔(以度为单位)、高度范围及间隔、反演区域平均地球半径
R=6371000;
jdmin_ = jdmin;wdmin_ = wdmin;
jdjg_ = jdjg;wdjg_ = wdjg;
jdmax_ = jdmax;wdmax_ = wdmax;
jdmin = deg2rad(jdmin);
jdjg = deg2rad(jdjg);
jdmax = deg2rad(jdmax);
wdmin = deg2rad(wdmin);
wdjg = deg2rad(wdjg);
wdmax = deg2rad(wdmax);
syms X Y Z
f1 = (X-XD)/(XS-XD)-(Y-YD)/(YS-YD);
f2 = (X-XD)/(XS-XD)-(Z-ZD)/(ZS-ZD);
f3 = (Y-YD)/(YS-YD)-(Z-ZD)/(ZS-ZD);
% 初始化结果存储数组
X0_coordinates = [];
Y0_coordinates = [];
Z0_coordinates = [];
% 生成高度序列
height_range = gdmin:gdjg:gdmax;
% 预分配单元格数组
num_heights = length(height_range);
X_cell = cell(1, num_heights);
Y_cell = cell(1, num_heights);
Z_cell = cell(1, num_heights);
% 计算循环
for idx = 1:num_heights
    current_H = height_range(idx);
    f4 = X^2 + Y^2 + Z^2 - (R+current_H)^2;
    % 求解方程组
    solution = solve(f1, f2, f3, f4, X, Y, Z);
    % 存储结果
    X_cell{idx} = double(solution.X);
    Y_cell{idx} = double(solution.Y);
    Z_cell{idx} = double(solution.Z);
end
% 合并所有结果
for idx = 1:num_heights
    X0_coordinates = [X0_coordinates, X_cell{idx}];
    Y0_coordinates = [Y0_coordinates, Y_cell{idx}];
    Z0_coordinates = [Z0_coordinates, Z_cell{idx}];
end

% 初始化经度面交点结果数组
X1_coordinates = [];
Y1_coordinates = [];
Z1_coordinates = [];
% 生成经度序列
longitude_range = deg2rad(lon1):jdjg:deg2rad(lon2);
% 预分配单元格数组
num_longitudes = length(longitude_range);
X_cell_long = cell(1, num_longitudes);
Y_cell_long = cell(1, num_longitudes);
Z_cell_long = cell(1, num_longitudes);
% 计算经度面交点
for idx = 1:num_longitudes
    current_J = longitude_range(idx);
    f6 = tan(current_J)*X - Y;  % 经度面方程 
    % 求解方程组
    solution = solve(f1, f2, f3, f6, X, Y, Z);
    % 存储结果
    X_cell_long{idx} = double(solution.X);
    Y_cell_long{idx} = double(solution.Y);
    Z_cell_long{idx} = double(solution.Z);
end
% 合并所有经度面结果
for idx = 1:num_longitudes
    X1_coordinates = [X1_coordinates, X_cell_long{idx}];
    Y1_coordinates = [Y1_coordinates, Y_cell_long{idx}];
    Z1_coordinates = [Z1_coordinates, Z_cell_long{idx}];
end

% 初始化纬度面交点结果数组
X2_coordinates = [];
Y2_coordinates = [];
Z2_coordinates = [];
% 生成纬度序列
latitude_range = deg2rad(lat1):wdjg:deg2rad(lat2);
% 预分配单元格数组
num_latitudes = length(latitude_range);
X_cell_lat = cell(1, num_latitudes);
Y_cell_lat = cell(1, num_latitudes);
Z_cell_lat = cell(1, num_latitudes);
% 计算纬度面交点
for idx = 1:num_latitudes
    current_W = latitude_range(idx);
    f5 = tan(current_W)^2*(X^2 + Y^2) - Z^2;  % 纬度面方程
    % 求解方程组
    solution = solve(f1, f2, f3, f5, X, Y, Z);
    % 处理可能的空解情况
    if ~isempty(solution.X)
        X_cell_lat{idx} = double(solution.X);
        Y_cell_lat{idx} = double(solution.Y); 
        Z_cell_lat{idx} = double(solution.Z);
    else
        X_cell_lat{idx} = [];
        Y_cell_lat{idx} = [];
        Z_cell_lat{idx} = [];
    end
end
% 合并所有纬度面结果
for idx = 1:num_latitudes
    X2_coordinates = [X2_coordinates, X_cell_lat{idx}];
    Y2_coordinates = [Y2_coordinates, Y_cell_lat{idx}];
    Z2_coordinates = [Z2_coordinates, Z_cell_lat{idx}];
end

JD_ddzb01 = [X0_coordinates(1,:);Y0_coordinates(1,:);Z0_coordinates(1,:)]';
JD_ddzb02 = [X0_coordinates(2,:);Y0_coordinates(2,:);Z0_coordinates(2,:)]';
JD_ddzb1 = [X1_coordinates(1,:);Y1_coordinates(1,:);Z1_coordinates(1,:)]';
JD_ddzb21 = [X2_coordinates(1,:);Y2_coordinates(1,:);Z2_coordinates(1,:)]';
JD_ddzb22 = [X2_coordinates(2,:);Y2_coordinates(2,:);Z2_coordinates(2,:)]';
JD_ddzb = [JD_ddzb01;JD_ddzb02;JD_ddzb1;JD_ddzb21;JD_ddzb22];
real_rows = all(imag(JD_ddzb) == 0, 2);  % 找出所有行都是实数的行
JD_ddzb = JD_ddzb(real_rows, :);         % 提取这些行
% 批量坐标转换
[lat, lon, alt] = XYZtoBLH_sphere(JD_ddzb(:,1), JD_ddzb(:,2), JD_ddzb(:,3),R);
JD_ddzb = [rad2deg(lat),rad2deg(lon),alt]';
epsilon = 1e-6;
% 检查第一行是否在 wdmin 和 wdmax 之间
valid_wd = JD_ddzb(1,:) >= rad2deg(wdmin) -epsilon & JD_ddzb(1,:) <= rad2deg(wdmax) +epsilon;
% 检查第二行是否在 jdmin 和 jdmax 之间
valid_jd = JD_ddzb(2,:) >= rad2deg(jdmin) -epsilon & JD_ddzb(2,:) <= rad2deg(jdmax) +epsilon;
% 检查第三行是否在 gdmin 和 gdmax 之间
valid_gd = JD_ddzb(3,:) >= gdmin -epsilon & JD_ddzb(3,:) <= gdmax +epsilon;
% 合并三个条件，找出所有符合条件的列索引
valid_columns = valid_wd & valid_jd & valid_gd;
%若没有有效交点则跳出函数
if isempty(valid_columns)
    return;
end
% 根据列索引获取符合条件的矩阵
valid_points = JD_ddzb(:,valid_columns)';

% 转换坐标
[X,Y,Z] = BLHtoXYZ_sphere(deg2rad(valid_points(:,1)),deg2rad(valid_points(:,2)),valid_points(:,3),R);
sortedjd = [X, Y, Z]';
sortedjd_ddzb = valid_points';

%%% 计算距离
distances = sqrt((sortedjd(1,:) - XD).^2 + (sortedjd(2,:) - YD).^2 + (sortedjd(3,:) - ZD).^2);

% 对距离进行排序
[distances, idx] = sort(distances);
% 重新排列 new_jd 矩阵
sortedjd = sortedjd(:, idx);
sortedjd_ddzb = sortedjd_ddzb(:, idx);

%%% 计算一维索引
if length(distances) < 2
    return;
end

% 计算两两交点的中点坐标
midpoints = (sortedjd(:,1:end-1) + sortedjd(:,2:end))/2;
% 批量转换中点坐标
[lat_mid, lon_mid, alt_mid] = XYZtoBLH_sphere(midpoints(1,:), midpoints(2,:), midpoints(3,:),R);
midpoint_blh = [rad2deg(lon_mid)', rad2deg(lat_mid)', alt_mid'];
% 计算两两交点之间的截距值
values = distances(2:end) - distances(1:end-1);
%% 索引计算优化
epsilon = 1e-5;
grid_indices = find_grid_indices(midpoint_blh, jdmin_, wdmin_, gdmin, jdjg_, wdjg_, gdjg, jdmax_, wdmax_, gdmax, epsilon);
%% 稀疏矩阵优化
total_length = round(((rad2deg(wdmax)-rad2deg(wdmin))/rad2deg(wdjg)) * ((rad2deg(jdmax)-rad2deg(jdmin))/rad2deg(jdjg)) * ((gdmax-gdmin)/gdjg));
one_dimensional_index = sparse(1, grid_indices, values, 1, total_length);

