function [E, A] = Get_EA(sx, sy, sz, x, y, z, R, sb, sl)
%Get_EA Calculates Elevation (E) and Azimuth (A) from reference to target.
%
% Usage:
%   [E,A] = Get_EA(sx,sy,sz, x,y,z, R)          % internally computes sb/sl
%   [E,A] = Get_EA(sx,sy,sz, x,y,z, R, sb, sl)  % reuse precomputed sb/sl (rad)

    % ---- Reuse sb/sl if provided, otherwise compute (same as original behavior) ----
    if nargin < 9
        [sb, sl] = XYZtoBLH_sphere(sx, sy, sz, R);
    end

    % ---- Pre-calculate trig ----
    sin_sb = sin(sb);
    cos_sb = cos(sb);
    sin_sl = sin(sl);
    cos_sl = cos(sl);

    % ---- Vector from sensor to target (avoid 3x1 array allocation) ----
    dx = x - sx;
    dy = y - sy;
    dz = z - sz;

    % ---- NEU components (same formula, just scalars instead of indexing) ----
    N = -sin_sb*cos_sl*dx - sin_sb*sin_sl*dy + cos_sb*dz; % North
    E_e = -sin_sl*dx + cos_sl*dy;                          % East
    U =  cos_sb*cos_sl*dx + cos_sb*sin_sl*dy + sin_sb*dz;  % Up

    % ---- Elevation ----
    E = atan2(U, sqrt(N^2 + E_e^2));

    % ---- Azimuth (0..2*pi) ----
    A = atan2(E_e, N);
    % if A < 0
    %     A = A + 2*pi;
    % end
    % 
    % % 避开 π 附近（极小偏移）
    % tol = 1e-12;
    % if abs(A - pi) < tol
    %     A = pi + tol;
    % end
end
