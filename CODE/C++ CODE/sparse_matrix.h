#pragma once
#include <string>
#include <vector>
#include <cstddef>

struct SparseMatrix {
    int rows;
    int cols;
    std::vector<int> row_indices;
    std::vector<int> col_indices;
    std::vector<double> values;

    SparseMatrix(int r, int c) : rows(r), cols(c) {}

    inline void reserve_nnz(size_t nnz) {
        row_indices.reserve(nnz);
        col_indices.reserve(nnz);
        values.reserve(nnz);
    }

    inline void add_value(int row, int col, double value) {
        row_indices.push_back(row);
        col_indices.push_back(col);
        values.push_back(value);
    }

    // 原始慢版本（可保留）
    void save_to_file(const std::string& filename) const;

    // 快速 txt 版本（推荐使用）
    void save_to_file_fast_txt(const std::string& filename,
                               bool one_based_index = false,
                               size_t buffer_bytes = 8 * 1024 * 1024) const;
};