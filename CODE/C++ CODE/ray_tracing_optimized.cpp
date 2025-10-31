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


void RayTracer::get_first_and_last_point_optimized(
    double jdmin_rad, double jdmax_rad, double wdmin_rad, double wdmax_rad,
    double gdmin, double gdmax, const Eigen::Vector3d& station, const Eigen::Vector3d& satellite,
    double jdjg_rad, double wdjg_rad, Eigen::Vector2d& lon_bounds, Eigen::Vector2d& lat_bounds, 
    bool& is_valid, Eigen::Vector3d& blh_station, double& A, double& z) {
    
    is_valid = false;
    double R_val = 6371000.0;
    
    // 坐标转换
    double B_U, L_U, H_U;
    XYZtoBLH_sphere(station[0], station[1], station[2], R_val, B_U, L_U, H_U);
    blh_station = Eigen::Vector3d(B_U, L_U, H_U);
    
    double E_temp;
    Get_EA(station[0], station[1], station[2], satellite[0], satellite[1], satellite[2], 
           R_val, E_temp, A);
    z = PI/2 - E_temp;
    
    // 数值稳定性检查
    if (std::isnan(z) || std::isnan(A)) {
        return;
    }
    
    // 预计算三角函数
    double sin_z = std::sin(z);
    double sin_B_U = std::sin(B_U);
    double cos_B_U = std::cos(B_U);
    double sin_A = std::sin(A);
    double cos_A = std::cos(A);
    
    // 检查数值有效性
    if (std::abs(sin_z) > 1.0 || std::isnan(sin_z)) {
        return;
    }
    
    // 向量化计算两个高度面的交点
    Eigen::Vector2d heights(gdmin, gdmax);
    Eigen::Vector2d sin_z_gdm = (R_val + H_U) * sin_z / (R_val + heights.array());
    
    // 数值边界检查
    if ((sin_z_gdm.array().abs() > 1.0).any() || sin_z_gdm.hasNaN()) {
        return;
    }
    
    Eigen::Vector2d z_gdm = sin_z_gdm.array().asin();
    Eigen::Vector2d alpha = z - z_gdm.array();
    
    // 向量化计算BLH坐标
    Eigen::Vector2d sin_alpha = alpha.array().sin();
    Eigen::Vector2d cos_alpha = alpha.array().cos();
    
    Eigen::Vector2d B_IPP_GD = (cos_alpha * sin_B_U + sin_alpha * cos_B_U * cos_A).array().asin();
    Eigen::Vector2d sin_dL = sin_alpha * sin_A;
    Eigen::Vector2d cos_B_IPP = B_IPP_GD.array().cos();
    
    // 避免除零和无效值
    if ((sin_dL.array().abs() > cos_B_IPP.array().abs()).any() || 
        B_IPP_GD.hasNaN() || cos_B_IPP.hasNaN()) {
        return;
    }
    
    Eigen::Vector2d dL = (sin_dL.array() / cos_B_IPP.array()).asin();
    Eigen::Vector2d L_IPP_GD = L_U + dL.array();
    
    // 边界检查
    if (!(check_in_bounds(B_IPP_GD[0], wdmin_rad, wdmax_rad) &&
          check_in_bounds(L_IPP_GD[0], jdmin_rad, jdmax_rad) &&
          check_in_bounds(B_IPP_GD[1], wdmin_rad, wdmax_rad) &&
          check_in_bounds(L_IPP_GD[1], jdmin_rad, jdmax_rad))) {
        return;
    }
    
    // 计算边界 
    double min_lon = std::min(L_IPP_GD[0], L_IPP_GD[1]);
    double max_lon = std::max(L_IPP_GD[0], L_IPP_GD[1]);
    double min_lat = std::min(B_IPP_GD[0], B_IPP_GD[1]);
    double max_lat = std::max(B_IPP_GD[0], B_IPP_GD[1]);
    
    // 近似函数
    lon_bounds[0] = std::max(approximateNumberDown(min_lon, jdjg_rad), jdmin_rad);
    lon_bounds[1] = std::min(approximateNumberUp(max_lon, jdjg_rad), jdmax_rad);
    lat_bounds[0] = std::max(approximateNumberDown(min_lat, wdjg_rad), wdmin_rad);
    lat_bounds[1] = std::min(approximateNumberUp(max_lat, wdjg_rad), wdmax_rad);
    
    // 确保边界有效
    if (lon_bounds[0] > lon_bounds[1] || lat_bounds[0] > lat_bounds[1]) {
        return;
    }
    
    is_valid = true;
}

void RayTracer::calculate_longitude_intersections(
    const Eigen::Vector3d& blh_station, double A, double z, double H_U,
    const Eigen::Vector2d& lon_bounds, double jdjg_rad, double gdmin, double gdmax,
    std::vector<Eigen::Vector3d>& intersections) {
    
    // 预分配内存
    intersections.reserve(intersections.size() + INITIAL_INTERSECTION_CAPACITY);
    
    double B_U = blh_station[0];
    double L_U = blh_station[1];
    
    // 生成经度网格
    int num_lons = std::max(1, static_cast<int>((lon_bounds[1] - lon_bounds[0]) / jdjg_rad) + 1);
    Eigen::VectorXd jd_grid = Eigen::VectorXd::LinSpaced(num_lons, lon_bounds[0], lon_bounds[1]);
    
    // 预计算常用值
    const double sin_z_val = std::sin(z);
    const double R_plus_H = R + H_U;
    
    // 向量化计算
    Eigen::VectorXd delta_Ls = jd_grid.array() - L_U;
    double tan_half_B = std::tan((PI/2 - B_U)/2);
    
    // 使用Eigen的array操作，避免临时对象
    auto A_minus_delta_Ls = (A - delta_Ls.array()) / 2;
    auto A_plus_delta_Ls = (A + delta_Ls.array()) / 2;
    
    Eigen::VectorXd cos_terms = A_minus_delta_Ls.cos() / A_plus_delta_Ls.cos();
    Eigen::VectorXd sin_terms = A_minus_delta_Ls.sin() / A_plus_delta_Ls.sin();
    
    // 处理数值问题
    cos_terms = cos_terms.unaryExpr([](double x) { return std::isfinite(x) ? x : 0.0; });
    sin_terms = sin_terms.unaryExpr([](double x) { return std::isfinite(x) ? x : 0.0; });
    
    Eigen::VectorXd sum_ab = 2 * (cos_terms * tan_half_B).array().atan();
    Eigen::VectorXd diff_ab = 2 * (sin_terms * tan_half_B).array().atan();
    
    Eigen::VectorXd a = (sum_ab + diff_ab) / 2;
    Eigen::VectorXd b = (sum_ab - diff_ab) / 2;
    
    Eigen::VectorXd B_IPP_JD = PI/2 - a.array();
    Eigen::VectorXd z_current = z - b.array();
    
    // 避免逐元素循环，使用向量化操作
    Eigen::VectorXd sin_z_current = z_current.array().sin();
    sin_z_current = sin_z_current.unaryExpr([](double x) { 
        return (std::abs(x) < 1e-15 || !std::isfinite(x)) ? 1e-15 : x; 
    });
    
    Eigen::VectorXd H_IPP_JD = (R_plus_H * sin_z_val) / sin_z_current.array() - R;
    
    // 批量添加有效点
    for (int i = 0; i < H_IPP_JD.size(); ++i) {
        double H_val = H_IPP_JD[i];
        if (H_val >= gdmin && H_val <= gdmax && std::isfinite(H_val) && 
            std::isfinite(B_IPP_JD[i]) && std::isfinite(jd_grid[i])) {
            intersections.emplace_back(B_IPP_JD[i], jd_grid[i], H_val);
        }
    }
}

void RayTracer::calculate_latitude_intersections(
    const Eigen::Vector3d& blh_station, double A, double z, double H_U,
    const Eigen::Vector2d& lat_bounds, double wdjg_rad, double gdmin, double gdmax,
    std::vector<Eigen::Vector3d>& intersections) {
    
    double B_U = blh_station[0];
    double L_U = blh_station[1];
    
    // 生成纬度网格 
    int num_lats = std::max(1, static_cast<int>((lat_bounds[1] - lat_bounds[0]) / wdjg_rad) + 1);
    Eigen::VectorXd wd_grid = Eigen::VectorXd::LinSpaced(num_lats, lat_bounds[0], lat_bounds[1]);
    
    // 数值稳定性检查
    if (std::isnan(B_U) || std::isnan(L_U) || std::isnan(A) || std::isnan(z)) {
        return;
    }
    
    for (int i = 0; i < wd_grid.size(); ++i) {
        double B_IPP = wd_grid[i];
        double b = PI/2 - B_U;
        double a = PI/2 - B_IPP;
        
        // 检查数值有效性
        if (std::abs(std::sin(a)) < 1e-15) {
            continue;
        }
        
        double sin_B = std::sin(b) / std::sin(a) * std::sin(A);
        if (std::abs(sin_B) > 1) continue;
        
        std::vector<double> B_candidates = {std::asin(sin_B), PI - std::asin(sin_B)};
        
        for (double B_val : B_candidates) {
            double A_plus_B = A + B_val;
            double A_minus_B = A - B_val;
            double a_plus_b = a + b;
            double a_minus_b = a - b;
            
            // 避免数值问题
            if (std::abs(std::cos(A_minus_B / 2)) < 1e-15 || 
                std::abs(std::tan(a_plus_b / 2)) > 1e10) {
                continue;
            }
            
            double tan_c_over_2 = std::cos(A_plus_B / 2) / std::cos(A_minus_B / 2) * std::tan(a_plus_b / 2);
            double c = 2 * std::atan(tan_c_over_2);
            
            if (std::abs(std::cos(a_plus_b / 2)) < 1e-15 || 
                std::abs(std::tan(A_plus_B / 2)) < 1e-15) {
                continue;
            }
            
            double tan_C_over_2 = std::cos(a_minus_b / 2) / std::cos(a_plus_b / 2) / std::tan(A_plus_B / 2);
            double C = 2 * std::atan(tan_C_over_2);
            
            double z_current = z - c;
            if (z_current <= 0) continue;
            
            double sin_z_current = std::sin(z_current);
            if (std::abs(sin_z_current) < 1e-15) {
                continue;
            }
            
            double H_val = ((H_U + R) * std::sin(z)) / sin_z_current - R;
            
            if (!std::isnan(H_val) && H_val >= gdmin && H_val <= gdmax) {
                double lon_val = L_U + C;
                if (!std::isnan(lon_val)) {
                    intersections.push_back(Eigen::Vector3d(B_IPP, lon_val, H_val));
                }
            }
        }
    }
}

void RayTracer::calculate_height_intersections(
    const Eigen::Vector3d& blh_station, double A, double z, double H_U,
    double gdmin, double gdmax, double gdjg,
    std::vector<Eigen::Vector3d>& intersections) {
    
    intersections.reserve(intersections.size() + INITIAL_INTERSECTION_CAPACITY);
    
    double B_U = blh_station[0];
    double L_U = blh_station[1];
    
    // 生成高度网格
    int num_heights = std::max(1, static_cast<int>((gdmax - gdmin) / gdjg) + 1);
    Eigen::VectorXd h_grid = Eigen::VectorXd::LinSpaced(num_heights, gdmin, gdmax);
    
    // 预计算常用值
    const double sin_z_val = std::sin(z);
    const double sin_B_U = std::sin(B_U);
    const double cos_B_U = std::cos(B_U);
    const double sin_A = std::sin(A);
    const double cos_A = std::cos(A);
    const double R_plus_H = R + H_U;
    
    // 向量化计算
    Eigen::VectorXd denominator = R + h_grid.array();
    Eigen::VectorXd sin_z_gdm = (R_plus_H * sin_z_val) / denominator.array();
    
    // 处理无效值
    sin_z_gdm = sin_z_gdm.unaryExpr([](double x) { 
        return (std::abs(x) > 1.0 || !std::isfinite(x)) ? 0.0 : x; 
    });
    
    Eigen::VectorXd z_gdm = sin_z_gdm.array().asin();
    Eigen::VectorXd alpha = z - z_gdm.array();
    
    Eigen::VectorXd sin_alpha = alpha.array().sin();
    Eigen::VectorXd cos_alpha = alpha.array().cos();
    
    // 计算B_IPP_GD
    Eigen::VectorXd B_IPP_GD = (cos_alpha * sin_B_U + sin_alpha * cos_B_U * cos_A).array().asin();
    Eigen::VectorXd sin_dL = sin_alpha * sin_A;
    Eigen::VectorXd cos_B_IPP = B_IPP_GD.array().cos();
    
    // 不使用OpenMP并行化（避免竞争条件和确保结果一致性）
    for (int i = 0; i < h_grid.size(); ++i) {
        if (!std::isfinite(B_IPP_GD[i]) || !std::isfinite(cos_B_IPP[i]) || 
            std::abs(cos_B_IPP[i]) < 1e-15) {
            continue;
        }
        
        if (std::abs(sin_dL[i]) > std::abs(cos_B_IPP[i])) continue;
        
        double dL = std::asin(sin_dL[i] / cos_B_IPP[i]);
        double L_IPP_GD = L_U + dL;
        
        if (std::isfinite(L_IPP_GD) && std::isfinite(B_IPP_GD[i])) {
            intersections.emplace_back(B_IPP_GD[i], L_IPP_GD, h_grid[i]);
        }
    }
}

RayResult RayTracer::trace_ray_optimized(
    const Eigen::Vector3d& station, const Eigen::Vector3d& satellite,
    double jdmin, double jdmax, double jdjg, double wdmin, double wdmax, double wdjg,
    double gdmin, double gdmax, double gdjg, double jdmin_rad, double wdmin_rad,
    double jdmax_rad, double wdmax_rad, double jdjg_rad, double wdjg_rad,
    int total_length, int nj, int nw, int ng) {
    
    RayResult result;
    
    Eigen::Vector2d lon_bounds, lat_bounds;
    bool is_valid;
    Eigen::Vector3d blh_station;
    double A, z;
    
    get_first_and_last_point_optimized(jdmin_rad, jdmax_rad, wdmin_rad, wdmax_rad, gdmin, gdmax,
                                     station, satellite, jdjg_rad, wdjg_rad,
                                     lon_bounds, lat_bounds, is_valid, blh_station, A, z);
    
    if (!is_valid) {
        return result;
    }
    
    double H_U = blh_station[2];
    std::vector<Eigen::Vector3d> intersections;
    
    // 预分配内存
    intersections.reserve(INITIAL_INTERSECTION_CAPACITY);
    
    // 计算三种类型的交点（不使用OpenMP sections，因为可能产生竞争条件）
    calculate_longitude_intersections(blh_station, A, z, H_U, lon_bounds, jdjg_rad, gdmin, gdmax, intersections);
    calculate_latitude_intersections(blh_station, A, z, H_U, lat_bounds, wdjg_rad, gdmin, gdmax, intersections);
    calculate_height_intersections(blh_station, A, z, H_U, gdmin, gdmax, gdjg, intersections);
    
    // 筛选有效点
    std::vector<Eigen::Vector3d> valid_points;
    valid_points.reserve(intersections.size());
    
    for (const auto& point : intersections) {
        if (point[0] >= wdmin_rad && point[0] <= wdmax_rad &&
            point[1] >= jdmin_rad && point[1] <= jdmax_rad &&
            point[2] >= gdmin && point[2] <= gdmax &&
            !point.hasNaN()) {
            valid_points.push_back(point);
        }
    }
    
    if (valid_points.size() < 2) {
        return result;
    }
    
    // 预分配结果内存
    result.grid_indices.reserve(valid_points.size() - 1);
    result.values.reserve(valid_points.size() - 1);
    
    // 批量计算距离和坐标转换
    std::vector<double> distances;
    distances.reserve(valid_points.size());
    
    for (const auto& point : valid_points) {
        double X, Y, Z;
        BLHtoXYZ_sphere(point[0], point[1], point[2], R, X, Y, Z);
        double dist = (Eigen::Vector3d(X, Y, Z) - station).norm();
        distances.push_back(dist);
    }
    
    // 使用索引排序，避免复制大量数据
    std::vector<size_t> indices(valid_points.size());
    std::iota(indices.begin(), indices.end(), 0);  // 现在应该可以工作了
    
    // 使用稳定排序（不使用并行排序，确保结果一致性）
    std::stable_sort(indices.begin(), indices.end(), 
        [&](size_t i, size_t j) { 
            // 如果距离非常接近，使用其他标准确保顺序一致
            if (std::abs(distances[i] - distances[j]) < 1e-10) {
                // 使用坐标作为次要排序标准
                if (std::abs(valid_points[i][0] - valid_points[j][0]) > 1e-10) {
                    return valid_points[i][0] < valid_points[j][0];
                }
                if (std::abs(valid_points[i][1] - valid_points[j][1]) > 1e-10) {
                    return valid_points[i][1] < valid_points[j][1];
                }
                return valid_points[i][2] < valid_points[j][2];
            }
            return distances[i] < distances[j]; 
        });
    
    // 批量计算中点和段长度
    for (size_t i = 0; i < indices.size() - 1; i++) {
        size_t idx1 = indices[i];
        size_t idx2 = indices[i + 1];
        
        Eigen::Vector3d mid_blh = (valid_points[idx1] + valid_points[idx2]) / 2;
        std::vector<double> mid_point_deg = {
            rad2deg(mid_blh[1]),
            rad2deg(mid_blh[0]), 
            mid_blh[2]
        };
        
        double segment_length = distances[idx2] - distances[idx1];
        
        if (segment_length > 0) {
            int grid_idx = find_grid_indices(mid_point_deg, jdmin, wdmin, gdmin, jdjg, wdjg, gdjg,
                                           jdmax, wdmax, gdmax, EPSILON, nj, nw, ng);
            
            if (grid_idx >= 1 && grid_idx <= total_length) {
                result.grid_indices.push_back(grid_idx);
                result.values.push_back(segment_length);
            }
        }
    }
    
    return result;
}