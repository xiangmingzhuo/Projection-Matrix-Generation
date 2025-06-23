function [lon1, lon2, lat1, lat2,  is_valid] = ...
         get_first_and_last_point(lon_min, lon_max, lat_min, lat_max, h_min, h_max, ...
                                 epsilon, max_iters, XD, YD, ZD, XS, YS, ZS, jdjg, wdjg)
R = 6371000;
% 初始化为有效射线
is_valid = true;
% 转换网格的8个角点到ECEF坐标
[~, Bmin, Bmax] = generate_grid_vertices(lon_min, lon_max, lat_min, lat_max, h_min, h_max);  
% 定义射线起点和终点 (ECEF坐标，单位：米)
P1 = [XD,YD,ZD];       % 起点
P2 = [XS,YS,ZS]; % 终点
% 计算方向向量
d = P2 - P1;
% 计算各轴的相交参数 
t_enter = (Bmin - P1) ./ d;
t_exit = (Bmax - P1) ./ d;
% 处理负方向的情况
% for j = 1:3
%     if d(j) < 0
%         [t_enter(j), t_exit(j)] = deal(t_exit(j), t_enter(j)); % 交换值
%     end
% end
% 处理负方向的情况 - 向量化处理
neg_dir = d < 0;
temp = t_enter(neg_dir);
t_enter(neg_dir) = t_exit(neg_dir);
t_exit(neg_dir) = temp;
% 计算全局相交区间
t_enter_global = max(t_enter);
t_exit_global = min(t_exit);
% 提前检查无效射线
if t_enter_global > 1 || t_exit_global < 0
    is_valid = false;lon1=[];lon2=[];lat1=[];lat2=[];
    return;
end
% 二分法确定首尾点
[intersection_point, outersection_point] = binary_search_ray_points(...
    Bmin, Bmax, P1, d, t_enter_global, t_exit_global, epsilon, max_iters);
% 批量转换两个点的坐标
intersection_point = intersection_point(:)';  % 转为行向量
outersection_point = outersection_point(:)';  % 转为行向量

% 拼接成 2×3 矩阵
points = [intersection_point; outersection_point];  

% 转换坐标
[lat, lon, h] = XYZtoBLH_sphere(points(:,1), points(:,2), points(:,3),R);
BLH = [rad2deg([lon, lat]), h];  % [经度, 纬度, 高度]

% 提取入口点和出口点
enter_point = BLH(1, :);
exit_point  = BLH(2, :);

% 检查射线是否完全在网格外 (优化条件判断)
if (enter_point(1) < lon_min && exit_point(1) < lon_min) || ...
   (enter_point(1) > lon_max && exit_point(1) > lon_max) || ...
   (enter_point(2) < lat_min && exit_point(2) < lat_min) || ...
   (enter_point(2) > lat_max && exit_point(2) > lat_max) || ...
   (enter_point(3) < h_min && exit_point(3) < h_min) || ...
   (enter_point(3) > h_max && exit_point(3) > h_max)
    is_valid = false;lon1=[];lon2=[];lat1=[];lat2=[];
    return;
end

% 检查高度条件 (提前返回)
if enter_point(3) > h_min || exit_point(3) < h_max
    is_valid = false;lon1=[];lon2=[];lat1=[];lat2=[];
    return;
end

% 计算包围盒边界
min_lon = min(enter_point(1), exit_point(1));
max_lon = max(enter_point(1), exit_point(1));
min_lat = min(enter_point(2), exit_point(2));
max_lat = max(enter_point(2), exit_point(2));

% 最终确定首尾点的经纬度高度值
lon1 = max(approximateNumberDown(min_lon, jdjg) - jdjg, lon_min);
lon2 = min(approximateNumberUp(max_lon, jdjg) + jdjg, lon_max);
lat1 = max(approximateNumberDown(min_lat, wdjg) - wdjg, lat_min);
lat2 = min(approximateNumberUp(max_lat, wdjg) + wdjg, lat_max);
end









