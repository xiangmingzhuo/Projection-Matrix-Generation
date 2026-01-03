#include "sparse_matrix.h"
#include <fstream>
#include <iostream>
#include <cstdio>
#include <string>

void SparseMatrix::save_to_file(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    for (size_t i = 0; i < values.size(); i++) {
        file << row_indices[i] << " " << col_indices[i] << " " << values[i] << "\n";
    }
}

void SparseMatrix::save_to_file_fast_txt(const std::string& filename,
                                         bool one_based_index,
                                         size_t buffer_bytes) const {
    // 用 binary 模式写文本可以避免某些平台换行转换带来的额外开销
    std::ofstream file(filename, std::ios::out | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string buffer;
    buffer.reserve(buffer_bytes);

    char line[128]; // 足够容纳 "row col value\n"
    const size_t nnz = values.size();

    for (size_t i = 0; i < nnz; ++i) {
        int r = row_indices[i];
        int c = col_indices[i];
        if (one_based_index) { r += 1; c += 1; } // 如果你希望与MATLAB一致

        // %.17g：double 友好且通常可无损往返
        int n = std::snprintf(line, sizeof(line), "%d %d %.17g\n", r, c, values[i]);
        buffer.append(line, (size_t)n);

        // 缓冲满了就写出一块
        if (buffer.size() >= buffer_bytes) {
            file.write(buffer.data(), (std::streamsize)buffer.size());
            buffer.clear();
        }
    }

    if (!buffer.empty()) {
        file.write(buffer.data(), (std::streamsize)buffer.size());
    }

    if (!file.good()) {
        std::cerr << "Error occurred while writing file: " << filename << std::endl;
    }
}