#include "coordinate_conversion.h"
#include "config.h"
#include <cmath>

void BLHtoXYZ_sphere(double B, double L, double H, double R, double& X, double& Y, double& Z) {
    double N = R + H;
    X = N * cos(B) * cos(L);
    Y = N * cos(B) * sin(L);
    Z = N * sin(B);
}

void XYZtoBLH_sphere(double X, double Y, double Z, double R, double& B, double& L, double& H) {
    L = atan2(Y, X);
    double P = sqrt(X * X + Y * Y);
    B = atan2(Z, P);
    H = sqrt(X * X + Y * Y + Z * Z) - R;
}

void Get_EA(double sx, double sy, double sz, double x, double y, double z, double R, double& E, double& A) {
    double sb, sl, sh;
    XYZtoBLH_sphere(sx, sy, sz, R, sb, sl, sh);  
    
    double sin_sb = sin(sb);
    double cos_sb = cos(sb);
    double sin_sl = sin(sl);
    double cos_sl = cos(sl);
    
    double deta_x = x - sx;
    double deta_y = y - sy;
    double deta_z = z - sz;
    
    double North = -sin_sb * cos_sl * deta_x - sin_sb * sin_sl * deta_y + cos_sb * deta_z;
    double East = -sin_sl * deta_x + cos_sl * deta_y;
    double Up = cos_sb * cos_sl * deta_x + cos_sb * sin_sl * deta_y + sin_sb * deta_z;
    
    E = atan2(Up, sqrt(North * North + East * East));
    A = atan2(East, North);
    
    if (A < 0) {
        A += 2 * PI;
    }
}