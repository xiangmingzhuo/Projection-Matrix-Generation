#pragma once

#include <cmath>

static constexpr double PI = 3.1415926535897932384626433832795;
static constexpr double R  = 6371000.0;  // 球半径(米)
static constexpr double EPSILON = 1e-10;
static constexpr int INITIAL_INTERSECTION_CAPACITY = 128;

inline double deg2rad(double deg) { return deg * PI / 180.0; }
inline double rad2deg(double rad) { return rad * 180.0 / PI; }
