#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <vector>
#include <string>

struct SparseMatrix {
    int rows;
    int cols;
    std::vector<int> row_indices;
    std::vector<int> col_indices;
    std::vector<double> values;
    
    SparseMatrix(int r, int c) : rows(r), cols(c) {}
    
    void add_value(int row, int col, double value);
    void save_to_file(const std::string& filename);
};

#endif