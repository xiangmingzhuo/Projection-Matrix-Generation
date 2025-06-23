%�˺��������Ե�������������A
function [one_dimensional_index,sortedjd_ddzb,values,jdmjd,wdmjd,gdmjd] =Get_sta_A1(XD,YD,ZD,XS,YS,ZS,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg)
R = 6371000;
% ��ʼ�����
[one_dimensional_index, sortedjd_ddzb, values, jdmjd, wdmjd, gdmjd] = deal([]);
%������β���㣬������β����ȷ����Ч����Χ
epsilon = 1e-10 ; max_iters = 100;
lon_min = jdmin ; lon_max = jdmax; lat_min = wdmin; lat_max = wdmax; h_min = gdmin; h_max = gdmax; 
[lon1,lon2,lat1,lat2,is_valid] = get_first_and_last_point(lon_min, lon_max, lat_min, lat_max, h_min, h_max,epsilon,max_iters,XD,YD,ZD,XS,YS,ZS,jdjg,wdjg);
% ���������Ч���ߣ�ֱ����������
if ~is_valid, return; end 
%%% ����ת��Ԥ����
[B_U, L_U, H_U] = XYZtoBLH_sphere(XD,YD,ZD,R);
[E,A] = Get_EA(XD,YD,ZD,XS,YS,ZS,R);
z = pi/2 - E;
R = 6371000;

%%% ���������������
jd_grid = deg2rad(approximateNumberDown(lon1,jdjg)):deg2rad(jdjg):deg2rad(approximateNumberUp(lon2,jdjg));
% ����������
if isempty(jd_grid), jd_grid = deg2rad(lon1); end
delta_Ls = jd_grid - L_U;

% Ԥ���㳣����
tan_half_B = tan((pi/2 - B_U)/2);

% ����������
cos_terms = cos((A - delta_Ls)/2) ./ cos((A + delta_Ls)/2);
sin_terms = sin((A - delta_Ls)/2) ./ sin((A + delta_Ls)/2);

% ����a��b
sum_ab = 2 * atan(cos_terms * tan_half_B);
diff_ab = 2 * atan(sin_terms * tan_half_B);

a = (sum_ab + diff_ab)/2;
b = (sum_ab - diff_ab)/2;

% ����IPP����
B_IPP_JD = pi/2 - a;
z_ = z - b;
H_IPP_JD = ((R+H_U)*sin(z)) ./ sin(z_) - R;

% ������������
valid_jd = (H_IPP_JD >= gdmin) & (H_IPP_JD <= gdmax);
jdmjd = [rad2deg(B_IPP_JD(valid_jd));rad2deg(jd_grid(valid_jd));  H_IPP_JD(valid_jd)]';

%%% ������γ�������
wd_grid = deg2rad(approximateNumberDown(lat1, wdjg)):deg2rad(wdjg):deg2rad(approximateNumberUp(lat2, wdjg));
% ����������
if isempty(wd_grid)
    wd_grid = deg2rad(lat1);
end
% ʹ����ֵ����
B_IPP = wd_grid(:)';
b = pi/2 - B_U;
a = pi/2 - B_IPP;
% ���� B ����������ֵ
sin_B = sin(b) ./ sin(a) .* sin(A);
B1 = asin(sin_B);
B2 = pi - B1;  % �ڶ������ܵ� B ֵ
% ������ B ֵ�ϲ�Ϊһ�� 2xN ����
B_all = [B1; B2];
% Ԥ�����м��������չά����ͬʱ��������⣩
A_plus_B = A + B_all;
A_minus_B = A - B_all;
a_plus_b = a + b;
a_minus_b = a - b;
% ���� c �� C��ͬʱ��������⣩
tan_c_over_2 = cos(A_plus_B ./ 2) ./ cos(A_minus_B ./ 2) .* tan(a_plus_b ./ 2);
c = 2 * atan(tan_c_over_2);
tan_C_over_2 = cos(a_minus_b ./ 2) ./ cos(a_plus_b ./ 2) .* cot(A_plus_B ./ 2);
C = 2 * atan(tan_C_over_2);
% ���� z_ �� H_��ͬʱ��������⣩
z_ = z - c;
H_ = ((H_U + 6371000) .* sin(z)) ./ sin(z_) - 6371000;
% ɸѡ��Ч�⣨ͬʱ����������ɸѡ��
valid = (H_ >= gdmin) & (H_ <= gdmax) & (imag(H_) == 0);
% ��ȡ��Ч�Ⲣ�ϲ����
wdmjd = [];
for k = 1:2
    current_valid = valid(k, :);
    if any(current_valid)
        jd = [rad2deg(wd_grid(current_valid)); ...
                 rad2deg(L_U + C(k, current_valid)); ...
                 H_(k, current_valid)]';
        wdmjd = [wdmjd; jd];
    end
end

%%% �������߶������
h_grid = gdmin:gdjg:gdmax;
% ����������
if isempty(h_grid), h_grid = gdmin; end

z_gdm = asin((R+H_U)*sin(z)./(R+h_grid));
alpha = z - z_gdm;
B_IPP_GD = asin(cos(alpha)*sin(B_U) + sin(alpha)*cos(B_U)*cos(A));
L_IPP_GD = L_U + asin(sin(alpha).*sin(A)./cos(B_IPP_GD));

% �����߶�����
valid_gd = (B_IPP_GD >= deg2rad(wdmin)) & (B_IPP_GD <= deg2rad(wdmax)) & ...
           (L_IPP_GD >= deg2rad(jdmin)) & (L_IPP_GD <= deg2rad(jdmax));
gdmjd = [rad2deg(B_IPP_GD(valid_gd));rad2deg(L_IPP_GD(valid_gd));  h_grid(valid_gd)]';

%%% �ϲ����н����
all_points = [jdmjd;wdmjd;gdmjd];
% ɸѡ��Ч����
valid_cols = (all_points(:,1) >= wdmin) & (all_points(:,1) <= wdmax) & ...
             (all_points(:,2) >= jdmin) & (all_points(:,2) <= jdmax) & ...
             (all_points(:,3) >= gdmin) & (all_points(:,3) <= gdmax);
valid_points = all_points(valid_cols,:);

% ת������
[X,Y,Z] = BLHtoXYZ_sphere(deg2rad(valid_points(:,1)),deg2rad(valid_points(:,2)),valid_points(:,3),R);
sortedjd = [X, Y, Z]';
sortedjd_ddzb = valid_points';

%%% �������
distances = sqrt((sortedjd(1,:) - XD).^2 + (sortedjd(2,:) - YD).^2 + (sortedjd(3,:) - ZD).^2);

% �Ծ����������
[distances, idx] = sort(distances);
% �������� new_jd ����
sortedjd = sortedjd(:, idx);
sortedjd_ddzb = sortedjd_ddzb(:, idx);

%%% ����һά����
if length(distances) < 2
    return;
end

% ��������������е�����
midpoints = (sortedjd(:,1:end-1) + sortedjd(:,2:end))/2;
% ����ת���е�����
[lat_mid, lon_mid, alt_mid] = XYZtoBLH_sphere(midpoints(1,:), midpoints(2,:), midpoints(3,:),R);
midpoint_blh = [rad2deg(lon_mid)', rad2deg(lat_mid)', alt_mid'];
% ������������֮��Ľؾ�ֵ
values = distances(2:end) - distances(1:end-1);
%% ���������Ż�
epsilon = 1e-5;
grid_indices = find_grid_indices(midpoint_blh, jdmin, wdmin, gdmin, jdjg, wdjg, gdjg, jdmax, wdmax, gdmax, epsilon);
%% ϡ������Ż�
total_length = round(((wdmax-wdmin)/wdjg) * ((jdmax-jdmin)/jdjg) * ((gdmax-gdmin)/gdjg));
one_dimensional_index = sparse(1, grid_indices, values, 1, total_length);