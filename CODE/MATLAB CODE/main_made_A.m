clear;clc;
%Inversion region
jdmin=0;jdmax=25;jdjg=5;
wdmin=35;wdmax=60;wdjg=5;
gdmin=100000;gdmax=1000000;gdjg=100000;
jdmin_rad=deg2rad(jdmin);wdmin_rad=deg2rad(wdmin);jdmax_rad=deg2rad(jdmax);wdmax_rad=deg2rad(wdmax);jdjg_rad=deg2rad(jdjg);wdjg_rad=deg2rad(wdjg);
total_length = round(((wdmax-wdmin)/wdjg) * ((jdmax-jdmin)/jdjg) * ((gdmax-gdmin)/gdjg));
Nx = (jdmax - jdmin) / jdjg;
Ny = (wdmax - wdmin) / wdjg;
Nz = (gdmax - gdmin) / gdjg;
jd_grid = jdmin_rad:jdjg_rad:jdmax_rad;
wd_grid = wdmin_rad:wdjg_rad:wdmax_rad;
h_grid = gdmin:gdjg:gdmax;
load data_5.mat;
% load valid_station_satellite_site_11.mat;
% station = stations ; satellite = satellites;
% load satellite.txt; load station.txt;
% satellite=[6827771.25089481,18077022.78276832,18537032.50751910];
% station=  [3275753.91200000,00321110.86510000,05445041.88290000];
% satellite=satellite(27659,:);
% station=station(27659,:);
% sx = station(1,1);sy=station(1,2);sz=station(1,3);
% sat_x = satellite(1,1);sat_y=satellite(1,2);sat_z=satellite(1,3);
% XD = station(1,1);YD=station(1,2);ZD=station(1,3);
% XS = satellite(1,1);YS=satellite(1,2);ZS=satellite(1,3);
% nj = (jdmax - jdmin) / jdjg; % Number of longitude cells
% nw = (wdmax - wdmin) / wdjg; % Number of latitude cells
% ng = (gdmax - gdmin) / gdjg; % Number of height cells
% nx = nj * nw * ng;           % Total number of cells in the 3D grid

% Retrieve the total number of CPU cores
% numCores = 5;  
% desiredCores = max(1, floor(1 * numCores));  
% % If the parallel pool is not started, it will be started based on the number of cores.
% if isempty(gcp('nocreate'))
%     parpool('local', desiredCores);
% end

%% 
% Traditional methods
% tic
% [A01, sortedjd01, distance01] = Get_A1_initial_simple1(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
% toc
% % Traditional methods（Effective cleavage plane）
% tic
% [A02, sortedjd02, distance02] = Get_A1_initial_simple2(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
% toc
% % Parametric Equation Method
tic
[A1, sortedjd1, distance1] = Get_A1_canshu_simple1(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
toc
% Parametric Equation Method（Effective cleavage plane）
tic
[A2, sortedjd2, distance2] = Get_A1_canshu_simple2(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
toc
% Spherical Trigonometry
tic
[A3, sortedjd3, distance3] = Get_A1_trangle_simple11(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length,jd_grid,wd_grid,h_grid);
toc
% Spherical Trigonometry（Effective cleavage planeⅠ）
tic
[A4, sortedjd4, distance4] = Get_A1_trangle_simple2(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length,h_grid);
toc
% Spherical Trigonometry（Effective cleavage planeⅡ）
tic
[A5, sortedjd5, distance5] = Get_A1_trangle_simple31(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length,h_grid);
toc
tic
[A6, sortedjd6, distance6] = Get_A1_canshu_sivt_simple2(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
toc
tic
[A7, sortedjd7, distance7, valid_indices] = Get_A1_trangle_simple_full(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length,jd_grid,wd_grid,h_grid);
toc
