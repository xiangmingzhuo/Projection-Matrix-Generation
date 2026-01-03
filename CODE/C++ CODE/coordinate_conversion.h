#pragma once

void BLHtoXYZ_sphere(double B, double L, double H, double R, double& X, double& Y, double& Z);
void XYZtoBLH_sphere(double X, double Y, double Z, double R, double& B, double& L, double& H);

void Get_EA(double sx, double sy, double sz,
            double x, double y, double z,
            double R, double& E, double& A);
