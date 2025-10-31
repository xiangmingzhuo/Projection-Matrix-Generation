#include "main_functions.h"
#include "config.h"
#include "ray_tracing_optimized.h"
#include "sparse_matrix.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <omp.h>

void Get_A1_trangle_simple3(const std::vector<std::vector<double>>& station,
                           const std::vector<std::vector<double>>& satellite,
                           double jdmin, double jdmax, double jdjg,
                           double wdmin, double wdmax, double wdjg,
                           double gdmin, double gdmax, double gdjg,
                           double jdmin_rad, double wdmin_rad,
                           double jdmax_rad, double wdmax_rad,
                           double jdjg_rad, double wdjg_rad,
                           int total_length, SparseMatrix& A) {
    
    int n_rays = station.size();
    int nj = static_cast<int>((jdmax - jdmin) / jdjg);
    int nw = static_cast<int>((wdmax - wdmin) / wdjg);
    int ng = static_cast<int>((gdmax - gdmin) / gdjg);
    
    // std::cout << "Starting to process " << n_rays << " rays..." << std::endl;
    
    // 使用线程安全的存储结构
    struct ThreadResult {
        int ray_index;
        std::vector<int> grid_indices;
        std::vector<double> values;
        bool valid;
    };
    
    std::vector<ThreadResult> thread_results(n_rays);
    
    // 初始化结果
    for (int i = 0; i < n_rays; i++) {
        thread_results[i].ray_index = i;
        thread_results[i].valid = false;
    }
    
    // OpenMP并行计算 - 保持顺序
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n_rays; i++) {
        Eigen::Vector3d sta_vec(station[i][0], station[i][1], station[i][2]);
        Eigen::Vector3d sat_vec(satellite[i][0], satellite[i][1], satellite[i][2]);
        
        try {
            RayResult ray_result = RayTracer::trace_ray_optimized(
                sta_vec, sat_vec, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg,
                gdmin, gdmax, gdjg, jdmin_rad, wdmin_rad, jdmax_rad, wdmax_rad,
                jdjg_rad, wdjg_rad, total_length, nj, nw, ng);
            
            if (!ray_result.grid_indices.empty()) {
                thread_results[i].grid_indices = ray_result.grid_indices;
                thread_results[i].values = ray_result.values;
                thread_results[i].valid = true;
            }
            
        } catch (const std::exception& e) {
            #pragma omp critical
            {
                std::cerr << "Ray " << i << " calculation failed: " << e.what() << std::endl;
            }
        }
        
        // // 进度显示
        // if ((i + 1) % std::max(1, n_rays / 20) == 0) {
        //     #pragma omp critical
        //     {
        //         std::cout << "Progress: " << (i + 1) << "/" << n_rays << std::endl;
        //     }
        // }
    }
    
    // 按顺序构建稀疏矩阵
    int valid_ray_count = 0;
    for (int i = 0; i < n_rays; i++) {
        if (thread_results[i].valid) {
            for (size_t j = 0; j < thread_results[i].grid_indices.size(); j++) {
                int col_index = thread_results[i].grid_indices[j] - 1; // MATLAB索引从1开始，C++从0开始
                double value = thread_results[i].values[j];
                
                if (col_index >= 0 && col_index < total_length) {
                    A.add_value(valid_ray_count, col_index, value);
                }
            }
            valid_ray_count++;
        }
    }
    
    // 更新稀疏矩阵的行数
    A.rows = valid_ray_count;
    
    // std::cout << "Processing completed. Valid rays: " << valid_ray_count << "/" << n_rays << std::endl;
    // std::cout << "OpenMP threads used: " << omp_get_max_threads() << std::endl;
}