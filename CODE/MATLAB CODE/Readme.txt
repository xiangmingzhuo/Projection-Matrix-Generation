This file contains some algorithms for generating projection matrices：
Midpoint Coordinate Index of a Line Segment:
1.parametric equation method；(Get_A1_canshu_simple1.m)
2.Parametric equations + effective splitting surface scheme 2；	(Get_A1_canshu_simple2.m)
3.spherical trigonometry；(Get_A1_trangle_simple11.m)
4.Spherical trigonometry + effective splitting surface scheme 1；	(Get_A1_trangle_simple2.m)
5.Spherical trigonometry + effective splitting surface scheme 2；	(Get_A1_trangle_simple31.m)

SIVT Algorithm Index:
6.Parametric equations + effective splitting surface scheme 2；	(Get_A1_canshu_sivt_simple2.m)

main_get_A.m : the main function
time_test.m : Record the computational time consumed by the algorithm
time_cost_map.m : Algorithmic computation time consuming visualization

input : 
	station.txt    : station coordinates
	satellite.txt : satellite coordinates
	or data_5.mat,data_11.mat : station and satellite coordinates
	jdmin、jdmax、jdjg、wdmin、wdmax、wdjg、gdmin、gdmax、gdjg 		 (Angular System)
	jdmin_rad、jdmax_rad、jdjg_rad、wdmin_rad、wdmax_rad、wdjg_rad 	 (Radian System)
	total_length: Total number of voxels
	t ：hour
output ： 	
	A : projection matrix
	sortedjd : The intersection of each ray with the grid
	distance : Intersection intercept of each ray with the grid

For any question, please contact Mingzhuo Xiang (1035952535@qq.com)