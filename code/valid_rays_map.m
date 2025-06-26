% 清除工作区并关闭所有图形
clear; close all; clc;
%导入测站和卫星数据
load station_and_satellite_11.mat;
%% 定义三维网格参数
lon_min = 0;lon_max = 25;grid_res1 = 26;
lat_min = 35;lat_max = 60;grid_res2 = 26;
alt_min = 100;alt_max = 1000;grid_res3 = 10;
grid_min = [lon_min, lat_min, alt_min];   
grid_max = [lon_max, lat_max, alt_max]; 
% R = 6371000;
% [grid_min(1,1),grid_min(1,2),grid_min(1,3)]=BLHtoXYZ_sphere(deg2rad(lat_min), deg2rad(lon_min), alt_min,R);
% [grid_max(1,1),grid_max(1,2),grid_max(1,3)]=BLHtoXYZ_sphere(deg2rad(lat_max), deg2rad(lon_max), alt_max,R);
% 生成网格点
x = linspace(grid_min(1), grid_max(1), grid_res1);
y = linspace(grid_min(2), grid_max(2), grid_res2);
z = linspace(grid_min(3), grid_max(3), grid_res3);
[X, Y, Z] = meshgrid(x, y, z);

%% 创建三维图形
figure('Color', 'white');
hold on; grid on; view(3);
axis equal; 
xlabel('X'); ylabel('Y'); zlabel('Z');
title('卫星-测站射线穿过三维网格示意图');

% 绘制三维网格
plotGrid(X, Y, Z);
daspect([1, 1, 100]);  % Z 轴缩放为 0.1 倍，使其看起来不那么长
% axis equal;
% 绘制测站和卫星
plot3(stations(:,1), stations(:,2), stations(:,3), 'bo', ...
      'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', '测站');
plot3(satellites(:,1), satellites(:,2), satellites(:,3), 'r^', ...
      'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', '卫星');

% 绘制所有射线路径
for i = 1:size(stations, 1)
    for j = 1:size(satellites, 1)
        % 计算射线参数化方程
        ray_start = stations(i, :);
        ray_end = satellites(j, :);
        direction = ray_end - ray_start;
        
        % 找到与网格的交点
        [entry, exit] = findGridIntersections(ray_start, direction, grid_min, grid_max);
        
        % 绘制射线段
        plot3([entry(1), exit(1)], ...
              [entry(2), exit(2)], ...
              [entry(3), exit(3)], ...
              'r-', 'LineWidth', 0.1, 'HandleVisibility', 'off');
    end
end
xlim([min([grid_min(1,1), grid_max(1,1)]),max([grid_min(1,1), grid_max(1,1)])]); % 限制 x 轴显示范围
ylim([min([grid_min(1,2), grid_max(1,2)]),max([grid_min(1,2), grid_max(1,2)])]); % 限制 y 轴显示范围
zlim([grid_min(1,3), grid_max(1,3)]); % 限制 z 轴显示范围
hold off;



