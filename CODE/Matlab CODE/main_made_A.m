
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
load data_11.mat;
% load satellite.txt; load station.txt;
% Retrieve the total number of CPU cores
numCores = 10;  
desiredCores = max(1, floor(1 * numCores));  
% If the parallel pool is not started, it will be started based on the number of cores.
if isempty(gcp('nocreate'))
    parpool('local', desiredCores);
end

%% part1 ： Use midpoint coordinates to obtain indices
% Traditional methods
tic
[A1, sortedjd, distance] = Get_A1_initial_simple1(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
toc
% Traditional methods（Effective cleavage plane）
tic
[A1, sortedjd, distance] = Get_A1_initial_simple2(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
toc
% Parametric Equation Method
tic
[A1, sortedjd, distance] = Get_A1_canshu_simple1(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
toc
% Parametric Equation Method（Effective cleavage plane）
tic
[A, sortedjd, distance] = Get_A1_canshu_simple2(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
toc
% Spherical Trigonometry
tic
[A1, sortedjd, distance] = Get_A1_trangle_simple1(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length);
toc
% Spherical Trigonometry（Effective cleavage planeⅠ）
tic
[A, sortedjd, distance] = Get_A1_trangle_simple2(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length);
toc
% Spherical Trigonometry（Effective cleavage planeⅡ）
tic
[A1, sortedjd1, distance1] = Get_A1_trangle_simple3(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length);
toc

%% %% part2 ： Use sivt Algorithm to obtain indices
% Traditional methods
tic
[A, sortedjd, distance] = Get_A1_initial_sivt_simple1(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
toc
% Traditional methods（Effective cleavage plane）
tic
[A, sortedjd, distance] = Get_A1_initial_sivt_simple2(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
toc
% Parametric Equation Method
tic
[A, sortedjd, distance] = Get_A1_canshu_sivt_simple1(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
toc
% Parametric Equation Method（Effective cleavage plane）
tic
[A1, sortedjd, distance] = Get_A1_canshu_sivt_simple2(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
toc
% Spherical Trigonometry
tic
[A, sortedjd, distance] = Get_A1_trangle_sivt_simple1(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length);
toc
% Spherical Trigonometry（Effective cleavage planeⅠ）
tic
[A, sortedjd, distance] = Get_A1_trangle_sivt_simple2(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length);
toc
% Spherical Trigonometry（Effective cleavage planeⅡ）
tic
[A, sortedjd, distance] = Get_A1_trangle_sivt_simple3(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length);
toc

%% get valid station and satellite sites 
[A, sortedjd, distance,station,satellite] = Get_A1_trangle_simple33(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length);
% 将测站卫星数据转为矩阵形式
fields = fieldnames(station);
stations = zeros(length(station.(fields{1})), length(fields))';
for i = 1:length(fields)
    stations(i,:) = station.(fields{i});
end
fields = fieldnames(satellite);
satellites = zeros(length(satellite.(fields{1})), length(fields))';
for i = 1:length(fields)
    satellites(i,:) = satellite.(fields{i});
end