function [X, Y, Z] = BLHtoXYZ_sphere(B, L, H, R)
    N = R + H;

    cB = cos(B);
    sB = sin(B);
    cL = cos(L);
    sL = sin(L);

    NcB = N .* cB;
    X = NcB .* cL;
    Y = NcB .* sL;
    Z = N .* sB;
end
