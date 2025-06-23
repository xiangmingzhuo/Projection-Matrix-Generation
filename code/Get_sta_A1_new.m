%此函数用来对单独的射线生成A
function [one_dimensional_index,sortedjd_ddzb,distances,jdmjd,wdmjd,gdmjd] =Get_sta_A1_new(XD,YD,ZD,XS,YS,ZS,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg)
%%%
%%%根据反演区域和尺度，基于某测站一条观测路径数据生成层析投影矩阵A
%%%

%预存储网格点
lon_min = jdmin ; lon_max = jdmax; lat_min = wdmin; lat_max = wdmax; h_min = gdmin; h_max = gdmax; 
%计算首尾交点，根据首尾交点确定有效网格范围
[lon1,lon2,lat1,lat2,h1,h2,is_valid] = get_first_and_last_point2(jdmin, jdmax, wdmin, wdmax, gdmin, gdmax,XD,YD,ZD,XS,YS,ZS,jdjg,wdjg);
%如果是有效射线
if is_valid
    %度转弧度
    jdmin = deg2rad(lon1);
    jdjg = deg2rad(jdjg);
    jdmax = deg2rad(lon2);
    wdmin = deg2rad(lat1);
    wdjg = deg2rad(wdjg);
    wdmax = deg2rad(lat2);
    
    [B_U, L_U, H_U] = XYZtoBLH_sphere(XD,YD,ZD,R);
    [E,A]= Get_EA(XD,YD,ZD,XS,YS,ZS);
    z = pi/2 - E;
    sortedjd_ddzb=[];
    distances=[];
    %经度面
    jdmjd = [];
    for i = jdmin:jdjg:jdmax    %遍历经度面的所有分割点
        L_IPP_JD = i;
        delta_L = L_IPP_JD - L_U;
        syms a b
        f1 = tan((a+b)/2)-(cos((A-delta_L)/2))/(cos((A+delta_L)/2))*tan((pi/2-B_U)/2);
        f2 = tan((a-b)/2)-(sin((A-delta_L)/2))/(sin((A+delta_L)/2))*tan((pi/2-B_U)/2);
        s = vpasolve(f1,f2,a,b);
        B_IPP_JD = pi/2 - s.a;
        z_ = z-s.b;
        H_IPP_JD =((6371000+H_U)*sin(z))/ sin(z_)-6371000;
        jd = [B_IPP_JD;L_IPP_JD;H_IPP_JD];
        jdmjd = double([jdmjd,jd]);
    end
    %纬度面
    wdmjd = [];
    for i = wdmin:wdjg:wdmax    %遍历纬度面的所有分割点
        B_IPP = i;
        b=pi/2-B_U;
        a=pi/2-B_IPP;
        B=asin(sin(b)/sin(a)*sin(A));
        syms c C
        f1 = tan(c/2)-(cos((A+B)/2))/(cos((A-B)/2))*tan((a+b)/2);
        f2 = tan(C/2)-(cos((a-b)/2))/(cos((a+b)/2))*cot((A+B)/2);
        s = vpasolve(f1,f2,c,C);
        z_ = z-s.c;
        H_=((H_U+6371000)*sin(z))/sin(z_)-6371000;
        jd = [B_IPP;s.C+L_U;H_];
        wdmjd = double([wdmjd,jd]);
    end
    %高度面
    gdmjd = [];
    for i = gdmin:gdjg:gdmax     %遍历高度面的所有分割点
        z_gdm = asin((6371000+H_U)*sin(z)/(6371000+i));
        alpha = z - z_gdm;
        B_IPP_GD = asin(cos(alpha)*sin(B_U)+sin(alpha)*cos(B_U)*cos(A));
        L_IPP_GD = L_U + asin((sin(alpha)*sin(A))/cos(B_IPP_GD));
        jd = [B_IPP_GD;L_IPP_GD;i];
        gdmjd = double([gdmjd,jd]);
    end
    %求高度面顶点和底点大地坐标,如果高度面顶点或底点经度或纬度不在范围内则跳出函数
    gdm_max = gdmjd(:,end);
    gdm_min = gdmjd(:,1);
    jd_gdm_max_ddzb = [gdm_max(1,1);gdm_max(2,1);gdm_max(3,1)];
    jd_gdm_min_ddzb = [gdm_min(1,1);gdm_min(2,1);gdm_min(3,1)];
    if jd_gdm_min_ddzb(1,1) < wdmin || jd_gdm_min_ddzb(1,1) > wdmax || jd_gdm_min_ddzb(2,1) < jdmin || jd_gdm_min_ddzb(2,1) > jdmax
            one_dimensional_index = [];
            return;
    end
    if jd_gdm_max_ddzb(1,1) < wdmin || jd_gdm_max_ddzb(1,1) > wdmax || jd_gdm_max_ddzb(2,1) < jdmin || jd_gdm_max_ddzb(2,1) > jdmax
            one_dimensional_index = [];
            return;
    end
    
    % 整合目前算的所有交点
    jd = [gdmjd,jdmjd,wdmjd];
    % 检查第一行是否在 wdmin 和 wdmax 之间
    valid_wd = (jd(1,:) >= deg2rad(lat_min)) & (jd(1,:) <= deg2rad(lat_max));
    % 检查第二行是否在 jdmin 和 jdmax 之间
    valid_jd = (jd(2,:) >= deg2rad(lon_min)) & (jd(2,:) <= deg2rad(lon_max));
    % 检查第三行是否在 gdmin 和 gdmax 之间
    valid_gd = (jd(3,:) >= h_min) & (jd(3,:) <= h_max);
    % 合并三个条件，找出所有符合条件的列索引
    valid_columns = valid_wd & valid_jd & valid_gd;
    % 根据列索引获取符合条件的矩阵
    new_jd_ddzb = jd(:,valid_columns);
    % 找到包含虚数的列
    has_complex = any(imag(new_jd_ddzb) ~= 0);
    % 剔除包含虚数的列
    new_jd_ddzb = new_jd_ddzb(:, ~has_complex);
    
    %获取矩阵大小
    [~, n] = size(new_jd_ddzb);
    % 初始化存储转换后坐标的矩阵new_jd
    new_jd = zeros(3, n);
    % 循环处理每一列数据
    for i = 1:n
        % 从new_jd_ddzb矩阵中提取第i列数据
        B = new_jd_ddzb(1, i);
        L = new_jd_ddzb(2, i);
        H = new_jd_ddzb(3, i);
        % 调用BLHtoXYZ_sphere函数将B、L、H转换为X、Y、Z
        [X, Y, Z] = BLHtoXYZ_sphere(B, L, H,R);
        % 将转换后的X、Y、Z存储到new_jd矩阵中
        new_jd(1, i) = X;
        new_jd(2, i) = Y;
        new_jd(3, i) = Z;
    end
    
    % 计算每个交点与测站的距离
    distances = sqrt((new_jd(1,:) - XD).^2 + (new_jd(2,:) - YD).^2 + (new_jd(3,:) - ZD).^2);
    % 对距离进行排序
    [~, idx] = sort(distances);
    % 重新排列 new_jd 矩阵
    sortedjd = new_jd(:, idx);
    sortedjd_ddzb = new_jd_ddzb(:, idx);
    %%%计算截距
    num_cols = size(sortedjd, 2);
    % 初始化一个向量来存储距离
    distances = zeros(1, num_cols - 1);
    % 计算相邻两点之间的距离
    for i = 1:num_cols - 1
        dx = sortedjd(1, i+1) - sortedjd(1, i);
        dy = sortedjd(2, i+1) - sortedjd(2, i);
        dz = sortedjd(3, i+1) - sortedjd(3, i);
        distances(i) = sqrt(dx^2 + dy^2 + dz^2);
    end
    
    %%%计算初始ID
    midpoint = [(sortedjd(1, 1)+sortedjd(1, 2))/2, (sortedjd(2, 1)+sortedjd(2, 2))/2,(sortedjd(3, 1)+sortedjd(3, 2))/2];
    [Bm,Lm,Hm]= XYZtoBLH_sphere(midpoint(1,1),midpoint(1,2),midpoint(1,3),R);
    s_id = [floor((Bm-deg2rad(lat_min))/wdjg)+1;floor((Lm-deg2rad(lon_min))/jdjg)+1;floor((Hm-h_min)/gdjg)+1];
    
    %%%计算标记矩阵 1为高度面交点，2为经度面交点，3为纬度面交点
    [~, n] = size(sortedjd_ddzb);
        bj_matrix = zeros(1, n);
        size_XX0 = size(gdmjd, 2);
        size_XX1 = size(jdmjd, 2);
        size_XX2 = size(wdmjd, 2);
        for i = 1:n
            column = sortedjd_ddzb(:, i);
            match = false;
            for j = 1:size(jd, 2)
                if isequal(column, jd(:, j))
                    if j <= size_XX0
                        bj_matrix(i) = 1;
                    elseif j <= size_XX0 + size_XX1
                        bj_matrix(i) = 2;
                    else
                        bj_matrix(i) = 3;
                    end
                    match = true;
                    break;
                end
            end
            if ~match
                error('无法找到匹配的列');
            end
        end
    bj_matrix = bj_matrix(:, 2:end);
    %%%
    ID = s_id;
    for j = 1:length(bj_matrix)
        if bj_matrix(j) == 1
            % 复制上一列到新列
            ID(:, j + 1) = ID(:, j);
            ID(3, j + 1) = ID(3, j) + 1; % 第三行等于前一列的第三行加1
        elseif bj_matrix(j) == 2 && sortedjd_ddzb(2, 1) > sortedjd_ddzb(2, 2)
            % 复制上一列到新列
            ID(:, j + 1) = ID(:, j);
            ID(2, j + 1) = ID(2, j) - 1; % 第二行-1
        elseif bj_matrix(j) == 2 && sortedjd_ddzb(2, 1) < sortedjd_ddzb(2, 2)
            % 复制上一列到新列
            ID(:, j + 1) = ID(:, j);
            ID(2, j + 1) = ID(2, j) +1; % 第二行+1
        elseif bj_matrix(j) == 3 && sortedjd_ddzb(1, 1) > sortedjd_ddzb(1, 2)
            % 复制上一列到新列
            ID(:, j + 1) = ID(:, j);
            ID(1, j + 1) = ID(1, j) - 1; % 第一行-1
        elseif bj_matrix(j) == 3 && sortedjd_ddzb(1, 1) < sortedjd_ddzb(1, 2)
            % 复制上一列到新列
            ID(:, j + 1) = ID(:, j);
            ID(1, j + 1) = ID(1, j) +1;
        else
            % 其他情况下，复制上一列到新列
            ID(:, j + 1) = ID(:, j);
        end
    end
    ID = ID(:, 1:end-1);
    result = [ID;distances];
    total_length = round(((deg2rad(lat_max)-deg2rad(lat_min))/wdjg)*((deg2rad(lon_max)-deg2rad(lon_min))/jdjg)*((h_max-h_min)/gdjg));
    
    % 初始化结果向量
    one_dimensional_index = sparse(1, total_length);
    
    % 遍历矩阵中的每一列
    n = size(result, 2);
    
    one_dim_idx_matrix = [];
    for i = 1:n
        % 获取当前列的三维编号和对应的值
        wd_num = result(1, i);
        jd_num = result(2, i);
        gd_num = result(3, i);
        value = result(4, i);
        
        % 转换为一维编号，并使用 round 函数确保是整数
        one_dim_idx = (gd_num - 1) * (round((deg2rad(lat_max)-deg2rad(lat_min))/wdjg)) * (round((deg2rad(lon_max)-deg2rad(lon_min))/jdjg)) + ...
                      (wd_num - 1) * (round((deg2rad(lon_max)-deg2rad(lon_min))/jdjg)) + ...
                      (jd_num - 1) + 1; % 加一是因为 MATLAB 索引从 1 开始
        one_dim_idx_matrix = [one_dim_idx_matrix,one_dim_idx];
        % 将第四行的值赋给结果向量中对应位置
        one_dimensional_index(one_dim_idx) = value;
    end
else
    %若为无效射线
    one_dimensional_index=[];
    sortedjd_ddzb=[];
    distances =[];
    jdmjd=[];
    wdmjd=[];
    gdmjd = [];
end