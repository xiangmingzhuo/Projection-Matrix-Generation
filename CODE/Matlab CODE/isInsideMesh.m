% Determines whether a point 'p' is inside a rectangular mesh defined by its minimum and maximum corners.
% @param p: The point to check [x, y, z].
% @param Bmin: The minimum corner of the bounding box [xmin, ymin, zmin].
% @param Bmax: The maximum corner of the bounding box [xmax, ymax, zmax].
% @return inside: 1 if the point is inside, 0 otherwise.
function inside = isInsideMesh(p, Bmin, Bmax)
    if p(1,1)>=Bmin(1,1) &&  p(1,1)<=Bmax(1,1) && p(1,2)>=Bmin(1,2) && p(1,2)<=Bmax(1,2) && p(1,3)>=Bmin(1,3) && p(1,3)<=Bmax(1,3)
        inside = 1;
    else
        inside = 0;
    end
end
