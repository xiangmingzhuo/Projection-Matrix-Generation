#include "sparse_matrix.h"
#include <fstream>
#include <iostream>

void SparseMatrix::add_value(int row, int col, double value) {
    row_indices.push_back(row);
    col_indices.push_back(col);
    values.push_back(value);
}

void SparseMatrix::save_to_file(const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    
    for (size_t i = 0; i < values.size(); i++) {
        file << row_indices[i] << " " << col_indices[i] << " " << values[i] << std::endl;
    }
    
    file.close();
}