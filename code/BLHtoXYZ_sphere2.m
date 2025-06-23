%此函数在球面坐标系啊将BLH转为XYZ
function [X, Y, Z] = BLHtoXYZ_sphere2(B, L, H)
% 将大地经度（B）、纬度（L）和椭球高（H）转换为直角坐标（X、Y、Z），假设地球是一个球体
B = deg2rad(B);
L = deg2rad(L);
% 地球半径
R = 6371000;

% 计算直角坐标X
X = (R + H) .* cos(B) .* cos(L);

% 计算直角坐标Y
Y = (R + H) .* cos(B) .* sin(L);

% 计算直角坐标Z
Z = (R + H) .* sin(B);