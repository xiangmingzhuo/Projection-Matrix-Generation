function [is_valid,lon1,lon2,lat1,lat2,h1,h2,orientVector] = get_first_and_last_point3 (XD,YD,ZD,XS,YS,ZS,jdjg,wdjg,jdmin,jdmax,wdmin,wdmax,gdmin,gdmax)
%初始化有效范围
lon1=[];lon2=[];lat1=[];lat2=[];h1=[];h2=[];
%计算射线方向向量
orientVector = [XS,YS,ZS] - [XD,YD,ZD]; % 计算方向向量
%计算首尾高度面的两个交点
R = 6371000;
A = sum(orientVector(1,:).^2);  % 更简洁的平方和计算
B = 2 * dot([XD, YD, ZD], orientVector(1,:));  % 使用点积简化计算
Cmin = XD^2 + YD^2 + ZD^2 - (R + gdmin)^2;
Cmax = XD^2 + YD^2 + ZD^2 - (R + gdmax)^2;
% 统一计算 tmin 和 tmax
discriminant_min = B^2 - 4*A*Cmin;
discriminant_max = B^2 - 4*A*Cmax;
% 使用向量化计算和逻辑索引
t_solutions_min = [(-B - sqrt(discriminant_min)); (-B + sqrt(discriminant_min))] / (2*A);
t_solutions_max = [(-B - sqrt(discriminant_max)); (-B + sqrt(discriminant_max))] / (2*A);
% 筛选符合条件的解
tmin = t_solutions_min((t_solutions_min > 0) & (t_solutions_min < 1));
tmax = t_solutions_max((t_solutions_max > 0) & (t_solutions_max < 1));
%计算首尾点的xyz坐标
intersection_point = [XD,YD,ZD] + orientVector * tmin;
outersection_point = [XD,YD,ZD] + orientVector * tmax;
%计算首尾点的经纬度高度值
[enter_point(:,2), enter_point(:,1), enter_point(:,3)] = XYZtoBLH_sphere(intersection_point(:,1), intersection_point(:,2), intersection_point(:,3),R);
enter_point  = [rad2deg(enter_point(:,1)),rad2deg(enter_point(:,2)),enter_point(:,3)];
[exit_point(:,2), exit_point(:,1), exit_point(:,3)] = XYZtoBLH_sphere(outersection_point(:,1), outersection_point(:,2), outersection_point(:,3),R);
exit_point  = [rad2deg(exit_point(:,1)),rad2deg(exit_point(:,2)),exit_point(:,3)];
%判断交点是否在网格内，否则视为无效射线
if enter_point(1,2) >= wdmin && enter_point(1,2) <= wdmax && enter_point(1,1) >= jdmin && enter_point(1,1) <= jdmax
    is_valid = true;
else
    is_valid =false;
    return;
end
if exit_point(1,2) >= wdmin && exit_point(1,2) <= wdmax && exit_point(1,1) >= jdmin && exit_point(1,1) <= jdmax
    is_valid = true;
else
    is_valid =false;
    return;
end
lon1 = max(approximateNumberDown(min(enter_point(:,1),exit_point(:,1)),jdjg)-jdjg,jdmin);
lon2 = min(approximateNumberUp(max(enter_point(:,1),exit_point(:,1)),jdjg)+jdjg,jdmax);
lat1 = max(approximateNumberDown(min(enter_point(:,2),exit_point(:,2)),wdjg)-wdjg,wdmin);
lat2 = min(approximateNumberUp(max(enter_point(:,2),exit_point(:,2)),wdjg)+wdjg,wdmax);
h1 = gdmin;
h2 = gdmax;

