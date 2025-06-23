# Projection-Matrix-Generation

This file contains seven algorithms for generating projection matrices：
1.Traditional methods；	(Get_sta_A1_initial.m)
2.Traditional Methods + Effective Split Surface Program 2；	(Get_sta_A1_initial_optimized.m)
3.parametric equation method；(Get_sta_A1_canshu.m)
4.Parametric equations + effective splitting surface scheme 2；	Get_sta_A1_canshu_optimized.m
5.spherical trigonometry；(Get_sta_A1.m)
6.Spherical trigonometry + effective splitting surface scheme 1；	(Get_sta_A1_tradition_optimized.m)
7.Spherical trigonometry + effective splitting surface scheme 2；	(Get_sta_A1_new_optimized.m)

main_get_A.m : the main function
time_test.m : Record the computational time consumed by the algorithm
time_cost_map.m : Algorithmic computation time consuming visualization

input : 
	doy : the time
	Sites_Info15.mat : Station data
	RDCB_ION.mat : Satellite differential code deviation
	SDCB_REF.mat : GPS Receiver Differential Code Offset Differential
    	P4 : Pseudo-distance observations for phase-smoothing codes
	sp3.mat : satellite coordinates
	jdmin、jdmax、jdjg、wdmin、wdmax、wdjg、gdmin、gdmax、gdjg
	t ：hour
output ： 	
	A : projection matrix
	L : STEC
	sortedjd : The intersection of each ray with the grid
	distance : Intersection intercept of each ray with the grid

For any question, please contact mingzhuo Xiang (1035952535@qq.com)
