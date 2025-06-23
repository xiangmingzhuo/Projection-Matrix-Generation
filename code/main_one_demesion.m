%输入卫星和测站xyz坐标，确定经纬度高度网格范围和大小
clear;clc;
satellite=[6827771.25089481,18077022.78276832,18537032.50751910];
station=  [3275753.91200000,00321110.86510000,05445041.88290000];
load station_and_satellite.mat;
satellite=ray_satellite.shexian6;
station=ray_station.shexian6;
XD=station(1,1);YD=station(1,2);ZD=station(1,3);
XS=satellite(1,1);YS=satellite(1,2);ZS=satellite(1,3);
jdmin=0;jdmax=25;jdjg=1;
wdmin=35;wdmax=60;wdjg=1;
gdmin=100000;gdmax=1000000;gdjg=100000;
%计算生成投影矩阵
tic
%二分法确定首尾点+球面三角
[one_dimensional_index1,sortedjd_ddzb1,distances1] =Get_sta_A1(XD,YD,ZD,XS,YS,ZS,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg);
toc
tic
%直接球面三角
[one_dimensional_index2,sortedjd_ddzb2,distances2] =Get_sta_A1_tradition_optimized(XD,YD,ZD,XS,YS,ZS,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg);
toc
tic
%球面三角确定首尾点+球面三角
[one_dimensional_index3,sortedjd_ddzb3,distances3] =Get_sta_A1_new_optimized(XD,YD,ZD,XS,YS,ZS,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg);
toc
tic
%参数方程法(有效分裂面筛选)
[one_dimensional_index4,sortedjd_ddzb4,distances4] =Get_sta_A1_canshu_optimized(XD,YD,ZD,XS,YS,ZS,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg);
toc
tic
%参数方程法(有效分裂面筛选隐式求解)
[one_dimensional_index5,sortedjd_ddzb5,distances5] =Get_sta_A1_initial_optimized(XD,YD,ZD,XS,YS,ZS,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg);
toc
tic
%参数方程
[one_dimensional_index6,sortedjd_ddzb6,distances6]=Get_sta_A1_canshu(XD,YD,ZD,XS,YS,ZS,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg);
toc

%变换一下单位
A=1e+3*one_dimensional_index;

%生成多条射线的投影矩阵和STEC
[A,L]=Get_A(doy,Sites_Info,sate,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,t);
