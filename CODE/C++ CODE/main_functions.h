#pragma once

#include <vector>
#include "sparse_matrix.h"

void Get_A1_trangle_simple3(const std::vector<std::vector<double>>& station,
                            const std::vector<std::vector<double>>& satellite,
                            double jdmin, double jdmax, double jdjg,
                            double wdmin, double wdmax, double wdjg,
                            double gdmin, double gdmax, double gdjg,
                            double jdmin_rad, double wdmin_rad,
                            double jdmax_rad, double wdmax_rad,
                            double jdjg_rad, double wdjg_rad,
                            int total_length, SparseMatrix& A);
