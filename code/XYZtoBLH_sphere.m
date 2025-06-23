%此函数用于在球面坐标系下将xyz转为BLH
function [B, L, H] = XYZtoBLH_sphere(X, Y, Z,R)
% 将直角坐标（X、Y、Z）转换为大地经度（B）、纬度（L）和椭球高（H），假设地球是一个球体

% 计算大地经度
L = atan2(Y, X);

% 计算平面距离
P = sqrt(X.^2 + Y.^2);

% 计算大地纬度
B = atan2(Z, P);

% 计算椭球高
H = sqrt(X.^2 + Y.^2 + Z.^2) - R;