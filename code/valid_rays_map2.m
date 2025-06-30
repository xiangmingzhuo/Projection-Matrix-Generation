clear; close all; clc;
load station_and_satellite_11.mat;

%% 网格参数
lon_min = 0; lon_max = 25; grid_res1 = 26;
lat_min = 35; lat_max = 60; grid_res2 = 26;
alt_min = 100; alt_max = 1000; grid_res3 = 10;
grid_min = [lon_min, lat_min, alt_min];   
grid_max = [lon_max, lat_max, alt_max]; 

%% 创建图形
figure('Color', 'white', 'Renderer', 'opengl', 'Position', [100 100 1200 800]);
hold on; grid on; view(3);
axis equal;
xlabel('Longitude (°)'); ylabel('Latitude (°)'); zlabel('Altitude (km)');
title('Valid satellite-station rays');

%% 增强网格绘制
x = linspace(grid_min(1), grid_max(1), grid_res1);
y = linspace(grid_min(2), grid_max(2), grid_res2);
z = linspace(grid_min(3), grid_max(3), grid_res3);

for k = [1:3:grid_res3]
    [X,Y] = meshgrid(x,y);
    surf(X, Y, repmat(z(k), size(X)), ...
        'FaceColor', 'none', ...
        'EdgeColor', [0.3 0.3 0.7], ...
        'LineWidth', 0.8, ...
        'EdgeAlpha', 0.5);
end

%% 随机分配测站颜色组（确保均匀分布）
num_stations = size(stations,1);
station_order = randperm(num_stations); % 随机打乱测站顺序
group_assignment = mod(station_order,3)+1; % 分成3组

% 定义颜色（红、绿、蓝，带透明度）
colors = {[1 0 0 0.15], [0 1 0 0.15], [0 0 1 0.15]};
group_colors = colors(group_assignment);

%% 绘制测站（按分组着色）
for i = 1:num_stations
    scatter3(stations(i,1), stations(i,2), stations(i,3), ...
             120, 'filled', 'MarkerFaceColor', group_colors{i}(1:3), ...
             'MarkerEdgeColor', 'k');
end

%% 绘制卫星（统一颜色）
scatter3(satellites(:,1), satellites(:,2), satellites(:,3), ...
         120, '^', 'filled', 'MarkerFaceColor', [1 0.3 0.3], ...
         'MarkerEdgeColor', 'k', 'DisplayName', '卫星');

%% 绘制射线（按测站分组着色）
for i = 1:num_stations
    current_color = group_colors{i};
    for j = 1:size(satellites,1)
        line([stations(i,1) satellites(j,1)], ...
             [stations(i,2) satellites(j,2)], ...
             [stations(i,3) satellites(j,3)], ...
             'Color', current_color, 'LineWidth', 0.01);
    end
end

%% 创建图例
% h_red = line([NaN NaN], [NaN NaN], [NaN NaN], 'Color', [1 0 0], 'LineWidth', 2);
% h_green = line([NaN NaN], [NaN NaN], [NaN NaN], 'Color', [0 1 0], 'LineWidth', 2);
% h_blue = line([NaN NaN], [NaN NaN], [NaN NaN], 'Color', [0 0 1], 'LineWidth', 2);
% legend([h_red, h_green, h_blue], {'Group 1', 'Group 2', 'Group 3'}, ...
%        'Location', 'northeast');

%% 图形美化
material dull;
camlight headlight;
lighting gouraud;

xlim([min(stations(:,1))-1 lon_max+1]);
ylim([min(stations(:,2))-1 lat_max+1]);
zlim([min(stations(:,3))-50 alt_max+200]);

daspect([1 1 100]);
set(gca, 'FontSize', 12, 'GridAlpha', 0.3);
hold off;