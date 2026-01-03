tic
[A1, sortedjd1, distance1] = Get_A1_canshu_simple1(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
toc
tic
[A1, sortedjd1, distance1] = Get_A1_canshu_simple2(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
toc
% Parametric Equation Method（Effective cleavage plane）
tic
[A2, sortedjd2, distance2] = Get_A1_trangle_simple11(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length,jd_grid,wd_grid,h_grid);
toc
tic
[A2, sortedjd1, distance1] = Get_A1_trangle_simple31(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length,h_grid);
toc
tic
[A2, sortedjd, distance] = Get_A1_trangle_simple2(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length,h_grid);
toc
% Spherical Trigonometry
tic
[A3, sortedjd3, distance3] = Get_A1_trangle_simple1(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length);
toc



% % 
% % % 方法1：获取所有字段名，然后取第3个
field_names = fieldnames(sortedjd);
third_field = field_names{20240}
% 
a = A2-A5;
% mean(a(:))
% a
max(a(:))