#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <omp.h>
#include "config.h"
#include "main_functions.h"
#include "sparse_matrix.h"

std::vector<std::vector<double>> readMatrixFromFile(const std::string& filename) {
    std::vector<std::vector<double>> matrix;
    std::ifstream file(filename);
    
    // if (!file.is_open()) {
    //     std::cerr << "Error opening file: " << filename << std::endl;
    //     return matrix;
    // }
    
    double x, y, z;
    while (file >> x >> y >> z) {
        matrix.push_back({x, y, z});
    }
    
    file.close();
    return matrix;
}

int main() {
    // 设置OpenMP线程数
    omp_set_num_threads(1);
    // 反演区域参数
    double jdmin = -10, jdmax = 50, jdjg = 2;
    double wdmin = 10, wdmax = 70, wdjg = 2;
    double gdmin = 100000, gdmax = 2100000, gdjg = 100000;
    
    double jdmin_rad = deg2rad(jdmin);
    double wdmin_rad = deg2rad(wdmin);
    double jdmax_rad = deg2rad(jdmax);
    double wdmax_rad = deg2rad(wdmax);
    double jdjg_rad = deg2rad(jdjg);
    double wdjg_rad = deg2rad(wdjg);
    
    int total_length = static_cast<int>(((wdmax - wdmin) / wdjg) * 
                                       ((jdmax - jdmin) / jdjg) * 
                                       ((gdmax - gdmin) / gdjg));
    
    // 读取数据
    std::vector<std::vector<double>> station = readMatrixFromFile("station.txt");
    std::vector<std::vector<double>> satellite = readMatrixFromFile("satellite.txt");
    
    // if (station.empty() || satellite.empty()) {
    //     std::cerr << "Error reading input files!" << std::endl;
    //     return 1;
    // }
    
    // 创建稀疏矩阵
    SparseMatrix A(station.size(), total_length);
    
    // 开始计时
    auto start = std::chrono::high_resolution_clock::now();
    
    // 主要计算
    Get_A1_trangle_simple3(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg,
                          gdmin, gdmax, gdjg, jdmin_rad, wdmin_rad, jdmax_rad, wdmax_rad,
                          jdjg_rad, wdjg_rad, total_length, A);
    
    // 结束计时
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // 保存时间
    std::ofstream time_file("time.txt");
    if (time_file.is_open()) {
        time_file << duration.count() << " ms" << std::endl;
        time_file.close();
    }
    
    // 保存投影矩阵
    A.save_to_file("project_matrix.txt");
    
    std::cout << "Calculation completed in " << duration.count() << " ms" << std::endl;
    std::cout << "Results saved to time.txt and project_matrix.txt" << std::endl;
    
    return 0;
}