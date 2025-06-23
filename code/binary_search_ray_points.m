% 合并的二分搜索函数
function [entry_point, exit_point] = binary_search_ray_points(...
    Bmin, Bmax, P1, d, t_entry, t_exit, epsilon, max_iters)
% 搜索入口点
entry_point = binary_search(@(t) is_inside_aabb(P1 + t*d, Bmin, Bmax), ...
                           t_entry, t_exit, epsilon, max_iters, true);
% 搜索出口点
exit_point = binary_search(@(t) is_inside_aabb(P1 + t*d, Bmin, Bmax), ...
                          t_entry, t_exit, epsilon, max_iters, false);
entry_point = P1 + entry_point * d;  % 返回第一个在网格内的点
exit_point = P1 + exit_point * d;  % 返回最后一个在网格内的点
end