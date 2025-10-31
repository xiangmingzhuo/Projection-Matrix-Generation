#include "math_utils.h"
#include <cmath>

double approximateNumberDown(double num, double interval) {
    return std::floor(num / interval) * interval;
}

double approximateNumberUp(double num, double interval) {
    return std::ceil(num / interval) * interval;
}

bool check_in_bounds(double value, double min_val, double max_val) {
    return (value >= min_val) && (value <= max_val);
}