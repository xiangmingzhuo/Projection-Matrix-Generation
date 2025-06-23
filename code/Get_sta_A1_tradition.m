%�˺��������Ե�������������A
function [one_dimensional_index,sortedjd_ddzb,distances,jdmjd,wdmjd,gdmjd] =Get_sta_A1_tradition(XD,YD,ZD,XS,YS,ZS,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg)
%%%
%%%���ݷ�������ͳ߶ȣ�����ĳ��վһ���۲�·���������ɲ���ͶӰ����A
%%%
jdmin = deg2rad(jdmin);
jdjg = deg2rad(jdjg);
jdmax = deg2rad(jdmax);
wdmin = deg2rad(wdmin);
wdjg = deg2rad(wdjg);
wdmax = deg2rad(wdmax);
[B_U, L_U, H_U] = XYZtoBLH_sphere(XD,YD,ZD,R);
[E,A]= Get_EA(XD,YD,ZD,XS,YS,ZS);
z = pi/2 - E;
sortedjd_ddzb=[];
distances=[];
%������
jdmjd = [];
for i = jdmin:jdjg:jdmax
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
%γ����
wdmjd = [];
for i = wdmin:wdjg:wdmax
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
%�߶���
gdmjd = [];
for i = gdmin:gdjg:gdmax
    z_gdm = asin((6371000+H_U)*sin(z)/(6371000+i));
    alpha = z - z_gdm;
    B_IPP_GD = asin(cos(alpha)*sin(B_U)+sin(alpha)*cos(B_U)*cos(A));
    L_IPP_GD = L_U + asin((sin(alpha)*sin(A))/cos(B_IPP_GD));
    jd = [B_IPP_GD;L_IPP_GD;i];
    gdmjd = double([gdmjd,jd]);
end
%��߶��涥��͵׵�������,����߶��涥���׵㾭�Ȼ�γ�Ȳ��ڷ�Χ������������
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

jd = [gdmjd,jdmjd,wdmjd];
% ����һ���Ƿ��� wdmin �� wdmax ֮��
valid_wd = jd(1,:) >= wdmin & jd(1,:) <= wdmax;
% ���ڶ����Ƿ��� jdmin �� jdmax ֮��
valid_jd = jd(2,:) >= jdmin & jd(2,:) <= jdmax;
% ���������Ƿ��� gdmin �� gdmax ֮��
valid_gd = jd(3,:) >= gdmin & jd(3,:) <= gdmax;
% �ϲ������������ҳ����з���������������
valid_columns = valid_wd & valid_jd & valid_gd;
% ������������ȡ���������ľ���
new_jd_ddzb = jd(:,valid_columns);
% �ҵ�������������
has_complex = any(imag(new_jd_ddzb) ~= 0);
% �޳�������������
new_jd_ddzb = new_jd_ddzb(:, ~has_complex);

%��ȡ�����С
[~, n] = size(new_jd_ddzb);
% ��ʼ���洢ת��������ľ���new_jd
new_jd = zeros(3, n);
% ѭ������ÿһ������
for i = 1:n
    % ��new_jd_ddzb��������ȡ��i������
    B = new_jd_ddzb(1, i);
    L = new_jd_ddzb(2, i);
    H = new_jd_ddzb(3, i);
    % ����BLHtoXYZ_sphere������B��L��Hת��ΪX��Y��Z
    [X, Y, Z] = BLHtoXYZ_sphere(B, L, H,R);
    % ��ת�����X��Y��Z�洢��new_jd������
    new_jd(1, i) = X;
    new_jd(2, i) = Y;
    new_jd(3, i) = Z;
end

distances = sqrt((new_jd(1,:) - XD).^2 + (new_jd(2,:) - YD).^2 + (new_jd(3,:) - ZD).^2);
% �Ծ����������
[~, idx] = sort(distances);
% �������� new_jd ����
sortedjd = new_jd(:, idx);
sortedjd_ddzb = new_jd_ddzb(:, idx);
%%%����ؾ�
num_cols = size(sortedjd, 2);
% ��ʼ��һ���������洢����
distances = zeros(1, num_cols - 1);
% ������������֮��ľ���
for i = 1:num_cols - 1
    dx = sortedjd(1, i+1) - sortedjd(1, i);
    dy = sortedjd(2, i+1) - sortedjd(2, i);
    dz = sortedjd(3, i+1) - sortedjd(3, i);
    distances(i) = sqrt(dx^2 + dy^2 + dz^2);
end

%%%�����ʼID
midpoint = [(sortedjd(1, 1)+sortedjd(1, 2))/2, (sortedjd(2, 1)+sortedjd(2, 2))/2,(sortedjd(3, 1)+sortedjd(3, 2))/2];
[Bm,Lm,Hm]= XYZtoBLH_sphere(midpoint(1,1),midpoint(1,2),midpoint(1,3),R);
s_id = [floor((Bm-wdmin)/wdjg)+1;floor((Lm-jdmin)/jdjg)+1;floor((Hm-gdmin)/gdjg)+1];

%%%�����Ǿ��� 1Ϊ�߶��潻�㣬2Ϊ�����潻�㣬3Ϊγ���潻��
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
            error('�޷��ҵ�ƥ�����');
        end
    end
bj_matrix = bj_matrix(:, 2:end);
%%%
ID = s_id;
for j = 1:length(bj_matrix)
    if bj_matrix(j) == 1
        % ������һ�е�����
        ID(:, j + 1) = ID(:, j);
        ID(3, j + 1) = ID(3, j) + 1; % �����е���ǰһ�еĵ����м�1
    elseif bj_matrix(j) == 2 && sortedjd_ddzb(2, 1) > sortedjd_ddzb(2, 2)
        % ������һ�е�����
        ID(:, j + 1) = ID(:, j);
        ID(2, j + 1) = ID(2, j) - 1; % �ڶ���-1
    elseif bj_matrix(j) == 2 && sortedjd_ddzb(2, 1) < sortedjd_ddzb(2, 2)
        % ������һ�е�����
        ID(:, j + 1) = ID(:, j);
        ID(2, j + 1) = ID(2, j) +1; % �ڶ���+1
    elseif bj_matrix(j) == 3 && sortedjd_ddzb(1, 1) > sortedjd_ddzb(1, 2)
        % ������һ�е�����
        ID(:, j + 1) = ID(:, j);
        ID(1, j + 1) = ID(1, j) - 1; % ��һ��-1
    elseif bj_matrix(j) == 3 && sortedjd_ddzb(1, 1) < sortedjd_ddzb(1, 2)
        % ������һ�е�����
        ID(:, j + 1) = ID(:, j);
        ID(1, j + 1) = ID(1, j) +1;
    else
        % ��������£�������һ�е�����
        ID(:, j + 1) = ID(:, j);
    end
end
ID = ID(:, 1:end-1);
result = [ID;distances];
total_length = round(((wdmax-wdmin)/wdjg)*((jdmax-jdmin)/jdjg)*((gdmax-gdmin)/gdjg));

% ��ʼ���������
one_dimensional_index = sparse(1, total_length);

% ���������е�ÿһ��
n = size(result, 2);

one_dim_idx_matrix = [];
for i = 1:n
    % ��ȡ��ǰ�е���ά��źͶ�Ӧ��ֵ
    wd_num = result(1, i);
    jd_num = result(2, i);
    gd_num = result(3, i);
    value = result(4, i);
    
    % ת��Ϊһά��ţ���ʹ�� round ����ȷ��������
    one_dim_idx = (gd_num - 1) * (round((wdmax-wdmin)/wdjg)) * (round((jdmax-jdmin)/jdjg)) + ...
                  (wd_num - 1) * (round((jdmax-jdmin)/jdjg)) + ...
                  (jd_num - 1) + 1; % ��һ����Ϊ MATLAB ������ 1 ��ʼ
    one_dim_idx_matrix = [one_dim_idx_matrix,one_dim_idx];
    % �������е�ֵ������������ж�Ӧλ��
    one_dimensional_index(one_dim_idx) = value;
end

% % ��ʼ�����ȡ�γ�ȡ��߶ȷ�Χ
% jd_range = jdmin:jdjg:jdmax;
% wd_range = wdmin:wdjg:wdmax;
% gd_range = gdmin:gdjg:gdmax;  
% % ���㾭γ�Ⱥ͸߶ȵ����ĵ�
% jd_centers = (jd_range(1:end-1) + jd_range(2:end)) / 2;
% wd_centers = (wd_range(1:end-1) + wd_range(2:end)) / 2;
% gd_centers = (gd_range(1:end-1) + gd_range(2:end)) / 2;
% % ��ȡ result ���������
% num_cols = size(one_dimensional_index, 2);
% % ��ʼ���洢���ĵ㾭γ�Ⱥ͸߶ȵľ���
% grid_centers = zeros(num_cols, 3);
% index = 1;
% % �����߶����ĵ㣬�ٱ���γ�����ĵ㣬�ٱ����������ĵ�
%     for g = 1:length(gd_centers)
%         for w = 1:length(wd_centers)
%             for j = 1:length(jd_centers)
%                 grid_centers(index, :) = [jd_centers(j), wd_centers(w), gd_centers(g)];
%                 index = index + 1;
%             end
%         end
%     end
