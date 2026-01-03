#include "grid_utils.h"
#include <cmath>
#include <vector>

int find_grid_indices(const std::vector<double>& midpoint_blh, 
                     double jdmin, double wdmin, double gdmin, 
                     double jdjg, double wdjg, double gdjg, 
                     double jdmax, double wdmax, double gdmax, 
                     double epsilon, int Nx, int Ny, int Nz) {
    
    double lon = midpoint_blh[0];
    double lat = midpoint_blh[1];
    double height = midpoint_blh[2];
    
    int indices_x = std::ceil((lon + epsilon - jdmin) / jdjg);
    int indices_y = std::ceil((lat + epsilon - wdmin) / wdjg);
    int indices_z = std::ceil((height + epsilon - gdmin) / gdjg);
    
    // 边界检查
    if (indices_x > Nx) indices_x = Nx;
    if (indices_y > Ny) indices_y = Ny;
    if (indices_z > Nz) indices_z = Nz;
    if (indices_x < 1) indices_x = 1;
    if (indices_y < 1) indices_y = 1;
    if (indices_z < 1) indices_z = 1;
    
    return (indices_z - 1) * (Nx * Ny) + (indices_y - 1) * Nx + indices_x;
}