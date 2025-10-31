#ifndef CONFIG_H
#define CONFIG_H

#include <cmath>
#include <Eigen/Dense>

// 常量定义
const double R = 6371000.0;
const double PI = 3.14159265358979323846;
const double EPSILON = 1e-15;

// 辅助函数
inline double deg2rad(double deg) { return deg * PI / 180.0; }
inline double rad2deg(double rad) { return rad * 180.0 / PI; }

// Eigen类型别名
using Vector3d = Eigen::Vector3d;
using Matrix3d = Eigen::Matrix3d;
using VectorXd = Eigen::VectorXd;
using MatrixXd = Eigen::MatrixXd;

#endif