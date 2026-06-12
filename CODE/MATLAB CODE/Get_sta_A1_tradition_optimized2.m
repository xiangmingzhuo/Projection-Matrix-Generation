% This function uses the direct spherical trigonometry method to generate matrix A for a single ray.
function [one_dimensional_index,sortedjd_ddzb,values,jdmjd,wdmjd,gdmjd] = Get_sta_A1_tradition_optimized2(XD,YD,ZD,XS,YS,ZS,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length,nj,nw,ng)
% Initialize outputs
R = 6371000;
[one_dimensional_index,  values ] = deal([]);
lon1 = jdmin_rad; lon2 = jdmax_rad;
lat1 = wdmin_rad; lat2 = wdmax_rad;
[B_U, L_U, H_U] = XYZtoBLH_sphere(XD,YD,ZD,R);
[E,A] = Get_EA(XD,YD,ZD,XS,YS,ZS,R);
z = pi/2 - E;
sinz = sin(z);  cosA = cos(A);  sinA = sin(A);
sinBU = sin(B_U); cosBU = cos(B_U);
RH = R+H_U;
%%% Vectorized calculation for longitude planes
jd_grid = lon1:jdjg_rad:lon2;
% Handle the case of a single point
if isempty(jd_grid), jd_grid = lon1; end
delta_Ls = jd_grid - L_U;

% Pre-calculate common terms
tan_half_B = tan((pi/2 - B_U)/2);

% Vectorized calculation
cos_terms = cos((A - delta_Ls)/2) ./ cos((A + delta_Ls)/2);
sin_terms = sin((A - delta_Ls)/2) ./ sin((A + delta_Ls)/2);

% Calculate a and b
sum_ab = 2 * atan(cos_terms * tan_half_B);
diff_ab = 2 * atan(sin_terms * tan_half_B);

a = (sum_ab + diff_ab)/2;
b = (sum_ab - diff_ab)/2;

% Calculate IPP parameters
B_IPP_JD = pi/2 - a;
z_ = z - b;
H_IPP_JD = ((RH)*sinz) ./ sin(z_) - R;

% Construct results for longitude planes
valid_jd = (H_IPP_JD >= gdmin) & (H_IPP_JD <= gdmax);
jdmjd = [B_IPP_JD(valid_jd);jd_grid(valid_jd);  H_IPP_JD(valid_jd)]';

%%% Vectorized calculation for latitude planes
wd_grid = lat1:wdjg_rad:lat2;
% Handle the case of a single point
if isempty(wd_grid)
    wd_grid = deg2rad(lat1);
end
% Use numerical methods
B_IPP = wd_grid(:)';
b = pi/2 - B_U;
a = pi/2 - B_IPP;
% Calculate two possible values for B
sin_B = sin(b) ./ sin(a) .* sinA;
B1 = asin(sin_B);
real_rows = all(imag(B1) == 0, 1);
B1 = B1(:, real_rows);
a = a(:,real_rows);
B2 = pi - B1;  % The second possible value for B
% Combine both sets of B values into a 2xN matrix
B_all = [B1; B2];
% Pre-calculate intermediate variables (expand dimensions to handle both sets of solutions simultaneously)
A_plus_B = A + B_all;
A_minus_B = A - B_all;
a_plus_b = a + b;
a_minus_b = a - b;
% Calculate c and C (handle both sets of solutions simultaneously)
tan_c_over_2 = cos(A_plus_B ./ 2) ./ cos(A_minus_B ./ 2) .* tan(a_plus_b ./ 2);
c = 2 * atan(tan_c_over_2);
tan_C_over_2 = (cos(a_minus_b/2) ./ cos(a_plus_b/2)) ./ tan(A_plus_B/2);
C = 2 * atan(tan_C_over_2);
% Calculate z_ and H_ (handle both sets of solutions simultaneously)
z_ = z - c;
H_ = ((H_U + 6371000) .* sinz) ./ sin(z_) - 6371000;
% Filter for valid solutions (filter both sets of solutions simultaneously)
valid = isfinite(H_) & (H_ >= gdmin) & (H_ <= gdmax) & isreal(H_);
% Extract valid solutions and merge results
% wdmjd = [];
% for k = 1:2
%     current_valid = valid(k, :);
%     if any(current_valid)
%         jd = [rad2deg(wd_grid(current_valid)); ...
%                  rad2deg(L_U + C(k, current_valid)); ...
%                  H_(k, current_valid)]';
%         wdmjd = [wdmjd; jd];
%     end
% end
if any(valid(:))
    % 扩展纬度网格以匹配 valid 的维度
    wd_grid_expanded = repmat(wd_grid, 2, 1);
    
    % 一次性计算所有转换
    valid_lat = wd_grid_expanded(valid);
    valid_lon = L_U + C(valid);  % 先相加再转换
    
    % 一次性构建结果矩阵（减少rad2deg调用）
    wdmjd = [valid_lat, ...           % 纬度
             valid_lon, ...           % 经度
             H_(valid)];                       % 高度
else
    wdmjd = [];
end
%%% Vectorized calculation for altitude planes
h_grid = gdmin:gdjg:gdmax;
% Handle the case of a single point
if isempty(h_grid), h_grid = gdmin; end
z_gdm = asin((RH)*sinz./(R+h_grid));
alpha = z - z_gdm;
B_IPP_GD = asin(cos(alpha)*sinBU + sin(alpha)*cosBU*cosA);
L_IPP_GD = L_U + asin(sin(alpha).*sinA./cos(B_IPP_GD));
gdmjd = [B_IPP_GD;L_IPP_GD;  h_grid]';
% Construct results for altitude planes
% valid_gd = (B_IPP_GD >= deg2rad(wdmin)) & (B_IPP_GD <= deg2rad(wdmax)) & ...
%            (L_IPP_GD >= deg2rad(jdmin)) & (L_IPP_GD <= deg2rad(jdmax));
% gdmjd = [rad2deg(B_IPP_GD(valid_gd));rad2deg(L_IPP_GD(valid_gd));  h_grid(valid_gd)]';

%%% Merge all result points
all_points = [jdmjd;wdmjd;gdmjd];
% Filter for valid intersection points
valid_cols = (all_points(:,1) >= wdmin_rad) & (all_points(:,1) <= wdmax_rad) & ...
             (all_points(:,2) >= jdmin_rad) & (all_points(:,2) <= jdmax_rad) & ...
             (all_points(:,3) >= gdmin) & (all_points(:,3) <= gdmax);
valid_points = all_points(valid_cols,:);
if isempty(valid_points)
    sortedjd_ddzb = [];
    return;
end
% Convert coordinates
[X,Y,Z] = BLHtoXYZ_sphere(valid_points(:,1),valid_points(:,2),valid_points(:,3),R);
sortedjd = [X, Y, Z]';
sortedjd_ddzb = valid_points';

%%% Calculate distances

% Reorder the sortedjd matrix
vx = XS - XD; vy = YS - YD; vz = ZS - ZD;
vnorm = sqrt(vx*vx + vy*vy + vz*vz);

% 对每个交点 P，t = dot(P-D, v) / |v|   （这是沿射线的有向距离）
t = ((sortedjd(1,:)-XD)*vx + (sortedjd(2,:)-YD)*vy + (sortedjd(3,:)-ZD)*vz) / vnorm;

[t, idx] = sort(t);
values = diff(t);   % 直接就是每段长度
sortedjd_ddzb = sortedjd_ddzb(:, idx);
if sortedjd_ddzb(3,1)> gdmin + 1e-10  || sortedjd_ddzb(3,end)< gdmax - 1e-10
    sortedjd_ddzb = [];
    return;
end    
%%% Calculate one-dimensional indices
% if length(distances) < 2
%     return;
% end

% Calculate the midpoint coordinates between each pair of intersection points
% midpoints = (sortedjd(:,1:end-1) + sortedjd(:,2:end))/2;
% % Batch convert midpoint coordinates
% [lat_mid, lon_mid, alt_mid] = XYZtoBLH_sphere(midpoints(1,:), midpoints(2,:), midpoints(3,:),R);
% midpoint_blh = [rad2deg(lon_mid)', rad2deg(lat_mid)', alt_mid'];
midpoint_blh =((sortedjd_ddzb(:,1:end-1) + sortedjd_ddzb(:,2:end))/2)';
RAD2DEG = 180/pi;
midpoint_blh(:,[1,2]) = midpoint_blh(:,[2,1]) * RAD2DEG;
% Calculate the intercept values between each pair of intersection points

%% Index calculation optimization
epsilon = 1e-15;
grid_indices = find_grid_indices(midpoint_blh, jdmin, wdmin, gdmin, jdjg, wdjg, gdjg, jdmax, wdmax, gdmax, epsilon,nj,nw,ng);
%% Sparse matrix optimization
% total_length = round(((wdmax-wdmin)/wdjg) * ((jdmax-jdmin)/jdjg) * ((gdmax-gdmin)/gdjg));
one_dimensional_index = sparse(1, grid_indices, values, 1, total_length);
