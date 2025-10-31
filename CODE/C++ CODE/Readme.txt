This file contains one algorithm for generating projection matrices：
Midpoint Coordinate Index of a Line Segment:
1.Spherical trigonometry + effective splitting surface scheme 2；	(main_functions_optimized.cpp)
input : 
	station.txt : station coordinates
	satellite.txt : satellite coordinates
	jdmin、jdmax、jdjg、wdmin、wdmax、wdjg、gdmin、gdmax、gdjg 		 (Angular System)
	jdmin_rad、jdmax_rad、jdjg_rad、wdmin_rad、wdmax_rad、wdjg_rad 	 (Radian System)
	total_length: Total number of voxels
	omp_set_num_threads : Set the number of OpenMP threads (set to 1 for single-threaded execution)
output ： 	
	project_matrix : projection matrix
	time.txt : the time cost

For any question, please contact Mingzhuo Xiang (1035952535@qq.com)