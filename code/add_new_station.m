function [insert_pts , sta_name2 ] = add_new_station(station_xyz)
station = station_xyz; %提取所有测站的坐标
% scatter3(station(:,1),station(:,2),station(:,3));
data = unique(station,'rows');  %筛选重复的测站坐标
x = data(:,1);y = data(:,2);z = data(:,3);  %提取xyz坐标
% 步骤1：计算原始数据的球坐标范围
[theta, phi, r] = cart2sph(x, y, z);
theta_range = [min(theta), max(theta)]; % 经度范围
phi_range = [min(phi), max(phi)];     % 纬度范围

% 步骤2：生成候选点（在原始数据范围内）
n_candidates = round(1/3*length(data(:,1)));
theta_new = theta_range(1) + rand(n_candidates,1) * diff(theta_range); 
phi_new = phi_range(1) + rand(n_candidates,1) * diff(phi_range);
r0 = mode(round(r, 3)); % 根据众数主球面半径

% 转换为笛卡尔坐标
[x_cand, y_cand, z_cand] = sph2cart(theta_new, phi_new, r0);
candidate_pts = [x_cand, y_cand, z_cand];

% 步骤3：密度自适应筛选
[~, dist] = knnsearch([x,y,z], candidate_pts, 'K', 3);
avg_dist = mean(dist, 2);
insert_threshold = prctile(avg_dist, 70); % 取距离最大的30%区域插入

% 筛选插入点
to_insert = avg_dist > insert_threshold;
insert_pts = candidate_pts(to_insert, :);

% 合并结果
final_pts = [x, y, z; insert_pts];

%给添加的虚拟测站取个名字
sta_name2 = [];
for i = 1:length(insert_pts(:,1))
    sta_name2{i,1} = (['XN',num2str(i)]);
end

% 可视化
figure;
subplot(1,2,1);
scatter3(x, y, z, 20, 'filled'); 
title(sprintf('原始数据 (%d点)', length(x)));
axis equal; view(45,30); xlim([min(x)-1 max(x)+1]); ylim([min(y)-1 max(y)+1]);

subplot(1,2,2);
scatter3(final_pts(:,1), final_pts(:,2), final_pts(:,3), 20, 'filled');
title(sprintf('区域插值后 (%d点)', size(final_pts,1)));
axis equal; view(45,30); xlim([min(x)-1 max(x)+1]); ylim([min(y)-1 max(y)+1]);

% % 显示覆盖范围对比
% fprintf('原始数据范围: θ=%.2f~%.2f rad, φ=%.2f~%.2f rad\n',...
%         theta_range(1), theta_range(2), phi_range(1), phi_range(2));