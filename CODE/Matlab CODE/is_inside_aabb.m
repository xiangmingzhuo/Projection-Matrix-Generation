% AABB内部检查函数
function inside = is_inside_aabb(point, Bmin, Bmax)
inside = all(point >= Bmin & point <= Bmax);
end