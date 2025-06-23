% 判断点是否在网格内的函数
function inside = isInsideMesh(p, Bmin,Bmax)
    if p(1,1)>=Bmin(1,1) &&  p(1,1)<=Bmax(1,1) && p(1,2)>=Bmin(1,2) && p(1,2)<=Bmax(1,2) && p(1,3)>=Bmin(1,3) && p(1,3)<=Bmax(1,3)
        inside = 1;
    else
        inside = 0;
    end
end