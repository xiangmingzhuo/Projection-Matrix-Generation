clear; close all; clc;
load station_and_satellite_11.mat;

%% 网格参数
lon_min = 0; lon_max = 25; grid_res1 = 26;
lat_min = 35; lat_max = 60; grid_res2 = 26;
alt_min = 100; alt_max = 1000; grid_res3 = 10;
grid_min = [lon_min, lat_min, alt_min];   
grid_max = [lon_max, lat_max, alt_max]; 

%% 创建图形（启用硬件加速）
figure('Color', 'white', 'Renderer', 'opengl', 'Position', [100 100 1200 800]);
hold on; grid on; view(3);
axis equal;
xlabel('Longitude (°)'); ylabel('Latitude (°)'); zlabel('Altitude (km)');
title('Satellite-station ray penetration grid');

%% 增强网格绘制 ---------------------------------------------------
% 生成网格点
x = linspace(grid_min(1), grid_max(1), grid_res1);
y = linspace(grid_min(2), grid_max(2), grid_res2);
z = linspace(grid_min(3), grid_max(3), grid_res3);

% 绘制主要网格平面（增强可见性）
for k = [1:3:grid_res3] % 只绘制底层、中层、顶层
    [X,Y] = meshgrid(x,y);
    surf(X, Y, repmat(z(k), size(X)), ...
        'FaceColor', 'none', ...
        'EdgeColor', [0.3 0.3 0.7], ... % 深蓝色网格线
        'LineWidth', 0.8, ...
        'EdgeAlpha', 0.5);
end

% 绘制垂直网格线（选择性显示）
for i = [1 floor(grid_res1/2) grid_res1]
    for j = [1 floor(grid_res2/2) grid_res2]
        line(repmat(x(i),1,2), repmat(y(j),1,2), [z(1) z(end)], ...
             'Color', [0.7 0.2 0.2], 'LineWidth', 0.5);
    end
end

%% 绘制测站和卫星（增强显示）
scatter3(stations(:,1), stations(:,2), stations(:,3), ...
         120, 'filled', 'MarkerFaceColor', [0 0.5 1], 'MarkerEdgeColor', 'k', ...
         'DisplayName', '地面测站');
scatter3(satellites(:,1), satellites(:,2), satellites(:,3), ...
         120, '^', 'filled', 'MarkerFaceColor', [1 0.3 0.3], 'MarkerEdgeColor', 'k', ...
         'DisplayName', '卫星');

%% 高效绘制所有射线
% 预计算所有射线端点
ray_pairs = zeros(2*size(stations,1)*size(satellites,1), 3);
idx = 1;
for i = 1:size(stations,1)
    for j = 1:size(satellites,1)
        ray_pairs(idx,:) = stations(i,:);
        ray_pairs(idx+1,:) = satellites(j,:);
        idx = idx + 2;
    end
end

% 单次绘制所有射线
h_rays = line(ray_pairs(:,1), ray_pairs(:,2), ray_pairs(:,3), ...
           'Color', [1 0 0 0.15], ... % 半透明红色
           'LineWidth', 0.00001, ...
           'HandleVisibility', 'off');

%% 图形美化
% 设置光照和材质（增强三维感）
material dull;
camlight headlight;
lighting gouraud;

% 坐标轴范围
xlim([min(stations(:,1)) lon_max]);
ylim([min(stations(:,2)) lat_max]);
zlim([min(stations(:,3)) alt_max+200]);

% 添加图例和比例尺
% legend('Location', 'northeast', 'FontSize', 12);
daspect([1 1 100]); % Z轴压缩比例
set(gca, 'FontSize', 12, 'GridAlpha', 0.3);
% colorbar('Ticks', [], 'Title', 'Altitude (km)');
hold off;