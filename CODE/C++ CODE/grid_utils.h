#pragma once
#include <vector>

int find_grid_indices(const std::vector<double>& midpoint_blh,
                      double jdmin, double wdmin, double gdmin,
                      double jdjg, double wdjg, double gdjg,
                      double jdmax, double wdmax, double gdmax,
                      double epsilon, int Nx, int Ny, int Nz);
