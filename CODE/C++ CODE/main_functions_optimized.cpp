#include "main_functions.h"
#include "ray_tracing_optimized.h"
#include <Eigen/Dense>
#include <omp.h>
#include <vector>

// 目标：在不改变A的计算逻辑/结果的前提下，通过：
// 1) 允许使用多线程（结果仍按原射线顺序写入）
// 2) 预统计NNZ并reserve，减少SparseMatrix反复扩容
// 3) 用move避免不必要的vector拷贝
// 来降低总耗时。

void Get_A1_trangle_simple3(const std::vector<std::vector<double>>& station,
                            const std::vector<std::vector<double>>& satellite,
                            double jdmin, double jdmax, double jdjg,
                            double wdmin, double wdmax, double wdjg,
                            double gdmin, double gdmax, double gdjg,
                            double jdmin_rad, double wdmin_rad,
                            double jdmax_rad, double wdmax_rad,
                            double jdjg_rad, double wdjg_rad,
                            int total_length, SparseMatrix& A) {

    const int n_rays = static_cast<int>(station.size());
    const int nj = static_cast<int>((jdmax - jdmin) / jdjg);
    const int nw = static_cast<int>((wdmax - wdmin) / wdjg);
    const int ng = static_cast<int>((gdmax - gdmin) / gdjg);

    // 每条射线的输出，按原index存储，保证最终写入顺序与单线程一致
    std::vector<RayResult> ray_results(n_rays);
    std::vector<unsigned char> valid(n_rays, 0);

    // 并行追踪射线
    // schedule(guided)通常比dynamic开销更小，同时保持负载均衡
    #pragma omp parallel for schedule(guided)
    for (int i = 0; i < n_rays; ++i) {
        Eigen::Vector3d sta_vec(station[i][0], station[i][1], station[i][2]);
        Eigen::Vector3d sat_vec(satellite[i][0], satellite[i][1], satellite[i][2]);

        RayResult r = RayTracer::trace_ray_optimized(
            sta_vec, sat_vec,
            jdmin, jdmax, jdjg,
            wdmin, wdmax, wdjg,
            gdmin, gdmax, gdjg,
            jdmin_rad, wdmin_rad,
            jdmax_rad, wdmax_rad,
            jdjg_rad, wdjg_rad,
            total_length, nj, nw, ng);

        if (!r.grid_indices.empty()) {
            ray_results[i] = std::move(r);
            valid[i] = 1;
        }
    }

    // 统计有效射线数量与总NNZ，用于reserve减少扩容（不影响结果）
    int valid_ray_count = 0;
    size_t total_nnz = 0;
    for (int i = 0; i < n_rays; ++i) {
        if (valid[i]) {
            ++valid_ray_count;
            total_nnz += ray_results[i].grid_indices.size();
        }
    }

    A.reserve_nnz(total_nnz);

    // 按原顺序构建稀疏矩阵（压缩无效射线的行号）
    int row_out = 0;
    for (int i = 0; i < n_rays; ++i) {
        if (!valid[i]) continue;

        const auto& cols_1based = ray_results[i].grid_indices;
        const auto& vals = ray_results[i].values;

        for (size_t k = 0; k < cols_1based.size(); ++k) {
            const int col_index = cols_1based[k] - 1; // MATLAB 1-based -> C++ 0-based
            if (col_index >= 0 && col_index < total_length) {
                A.add_value(row_out, col_index, vals[k]);
            }
        }
        ++row_out;
    }

    A.rows = valid_ray_count;
}
