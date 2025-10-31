#ifndef RAY_TRACING_OPTIMIZED_H
#define RAY_TRACING_OPTIMIZED_H

#include "ray_tracing_optimized.h"
#include "config.h"
#include "math_utils.h"
#include "coordinate_conversion.h"
#include "grid_utils.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include <numeric>  // 添加这个头文件用于std::iota

struct RayResult {
    std::vector<int> grid_indices;
    std::vector<double> values;
    
    RayResult() = default;
};

class RayTracer {
public:
    static void get_first_and_last_point_optimized(
        double jdmin_rad, double jdmax_rad, double wdmin_rad, double wdmax_rad,
        double gdmin, double gdmax, const Eigen::Vector3d& station, const Eigen::Vector3d& satellite,
        double jdjg_rad, double wdjg_rad, Eigen::Vector2d& lon_bounds, Eigen::Vector2d& lat_bounds, 
        bool& is_valid, Eigen::Vector3d& blh_station, double& A, double& z);

    static RayResult trace_ray_optimized(
        const Eigen::Vector3d& station, const Eigen::Vector3d& satellite,
        double jdmin, double jdmax, double jdjg, double wdmin, double wdmax, double wdjg,
        double gdmin, double gdmax, double gdjg, double jdmin_rad, double wdmin_rad,
        double jdmax_rad, double wdmax_rad, double jdjg_rad, double wdjg_rad,
        int total_length, int nj, int nw, int ng);

private:
    static void calculate_longitude_intersections(
        const Eigen::Vector3d& blh_station, double A, double z, double H_U,
        const Eigen::Vector2d& lon_bounds, double jdjg_rad, double gdmin, double gdmax,
        std::vector<Eigen::Vector3d>& intersections);

    static void calculate_latitude_intersections(
        const Eigen::Vector3d& blh_station, double A, double z, double H_U,
        const Eigen::Vector2d& lat_bounds, double wdjg_rad, double gdmin, double gdmax,
        std::vector<Eigen::Vector3d>& intersections);

    static void calculate_height_intersections(
        const Eigen::Vector3d& blh_station, double A, double z, double H_U,
        double gdmin, double gdmax, double gdjg,
        std::vector<Eigen::Vector3d>& intersections);

    static constexpr size_t INITIAL_INTERSECTION_CAPACITY = 1000;
};

#endif