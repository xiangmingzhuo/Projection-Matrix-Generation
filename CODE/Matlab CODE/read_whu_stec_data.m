clear;clc;
% 指定分隔符为逗号
data = readtable('D:\weixin\xwechat_files\wxid_x3qbnf2brxad22_abe1\msg\file\2025-09\EUREF_2015-12-16(1).txt', 'Delimiter', ',');
% 去掉俄罗斯的卫星
data = data(~startsWith(string(data{:,3}), 'R'), :);
% 转换时间列
data.Time = datetime(data{:, 1}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
% 确定层析区域
jdmin=0;jdmax=25;jdjg=1;
wdmin=35;wdmax=60;wdjg=1;
gdmin=100000;gdmax=1000000;gdjg=1000;
jdmin_rad=deg2rad(jdmin);wdmin_rad=deg2rad(wdmin);jdmax_rad=deg2rad(jdmax);wdmax_rad=deg2rad(wdmax);jdjg_rad=deg2rad(jdjg);wdjg_rad=deg2rad(wdjg);
total_length = round(((wdmax-wdmin)/wdjg) * ((jdmax-jdmin)/jdjg) * ((gdmax-gdmin)/gdjg));
epsilon = 1e-10; max_iters = 100;
valid_rays_id = [];
nj=(jdmax-jdmin)/jdjg;
nw=(wdmax-wdmin)/wdjg;
ng=(gdmax-gdmin)/gdjg;
nx=nj*nw*ng;
for t = 1:23
    % 直接定义时间范围
    start_time = datetime(['2015-12-16 ',sprintf('%02d', t-1),':53:00'], 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    end_time = datetime(['2015-12-16 ',sprintf('%02d', t),':08:00'], 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    
    % 筛选数据
    time_mask = data.Time >= start_time & data.Time <= end_time;
    filtered_data = data(time_mask, :);
    % 检查哪些列是数字类型
    numericCols = varfun(@isnumeric, filtered_data, 'OutputFormat', 'uniform');
    filtered_data = filtered_data(:, numericCols);
    filtered_data = table2array(filtered_data);
    A = sparse([], [], [], 0, 0);  % 创建0x0的空稀疏矩阵
    station = filtered_data(:,1:3); satellite = filtered_data(:,6:8)*1000;
    % [is_valid,~,~,~,~,~,~,~] = get_first_and_last_point3 (filtered_data(i,1),filtered_data(i,2),filtered_data(i,3),filtered_data(i,6)*1000,filtered_data(i,7)*1000,filtered_data(i,8)*1000,jdjg,wdjg,jdmin,jdmax,wdmin,wdmax,gdmin,gdmax);
    % if is_valid
    %     valid_rays_id = [valid_rays_id;i];
    % end

    % [A1,sorted_ddzb,distances,station,satellite,scope_blh_one_dimension] = Get_sta_A1_optimized(filtered_data(i,1),filtered_data(i,2),filtered_data(i,3),...
    %     filtered_data(i,6)*1000,filtered_data(i,7)*1000,filtered_data(i,8)*1000,...
    %     jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg);


    [A, sortedjd1, distance1,valid_rays_id] = Get_A1_trangle_simple3(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length);
    % if ~isempty(A1) && max(sorted_ddzb(3,:))<gdmax-100
    %     valid_rays_id = [valid_rays_id;i];
    % end
    % A = [A;A1];
    filtered_data = filtered_data(valid_rays_id,:);

    %保存数据
    fname = sprintf('stec_%d.mat',t);
    savePath = 'C:\Users\Dell\Desktop\GNSS6';
    % 保存数据
    save(fullfile(savePath,fname),'filtered_data','A');
    clear A filtered_data ;
end