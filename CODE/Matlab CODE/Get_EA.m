function [E, A] = Get_EA(sx, sy, sz, x, y, z, R)
%Get_EA Calculates the Elevation and Azimuth angles from a reference point to a target point.
%
%   Inputs:
%       sx, sy, sz - ECEF coordinates of the reference (sensor) point.
%       x, y, z    - ECEF coordinates of the target point.
%       R          - Radius of the spherical Earth model.
%
%   Outputs:
%       E - Elevation angle in radians (positive above the horizon).
%       A - Azimuth angle in radians, measured clockwise from True North.
%
%   This function first converts the reference point's ECEF coordinates to
%   geodetic (latitude/longitude) to define the local North-East-Up (NEU)
%   coordinate system. It then transforms the vector from the reference to
%   the target into this NEU frame and calculates the spherical angles.

    % Convert the reference point's ECEF coordinates to geodetic (latitude, longitude).
    % We only need latitude (sb) and longitude (sl) to define the local NEU basis.
    [sb, sl] = XYZtoBLH_sphere(sx, sy, sz, R);
    
    % Pre-calculate sine and cosine values for efficiency.
    sin_sb = sin(sb);
    cos_sb = cos(sb);
    sin_sl = sin(sl);
    cos_sl = cos(sl);
    
    % Calculate the vector from the reference point to the target point in ECEF coordinates.
    deta_xyz = [x - sx; y - sy; z - sz];
    
    % Transform the ECEF vector into the local North-East-Up (NEU) coordinate frame.
    % This is equivalent to R_NEU * deta_xyz, but written out for speed.
    NEU = zeros(3, 1);
    NEU(1) = -sin_sb*cos_sl*deta_xyz(1) - sin_sb*sin_sl*deta_xyz(2) + cos_sb*deta_xyz(3); % North component
    NEU(2) = -sin_sl*deta_xyz(1) + cos_sl*deta_xyz(2);                                     % East component
    NEU(3) = cos_sb*cos_sl*deta_xyz(1) + cos_sb*sin_sl*deta_xyz(2) + sin_sb*deta_xyz(3);   % Up component
    
    % Calculate the Elevation angle.
    % Using atan2(Up, Horizontal_Distance) is numerically stable and avoids division by zero.
    E = atan2(NEU(3), sqrt(NEU(1)^2 + NEU(2)^2));
    
    % Calculate the Azimuth angle.
    % atan2(East, North) gives the angle measured counter-clockwise from the North axis.
    A = atan2(NEU(2), NEU(1));
    
    % Convert the azimuth to the standard [0, 2*pi) range (clockwise from North).
    if A < 0
        A = A + 2*pi;
    end
end
