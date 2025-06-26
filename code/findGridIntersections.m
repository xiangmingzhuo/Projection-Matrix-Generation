%% 辅助函数：计算射线与网格的交点
function [entry, exit] = findGridIntersections(origin, direction, grid_min, grid_max)
    % 计算射线与网格六个面的交点
    t_enter = -inf(1,3);
    t_exit = inf(1,3);
    
    for i = 1:3
        if abs(direction(i)) > eps
            t1 = (grid_min(i) - origin(i)) / direction(i);
            t2 = (grid_max(i) - origin(i)) / direction(i);
            t_enter(i) = min(t1, t2);
            t_exit(i) = max(t1, t2);
        end
    end
    
    % 找到最近的进入点和最远的离开点
    t_entry = max(t_enter);
    t_exit = min(t_exit);
    
    % 计算交点坐标
    entry = origin + t_entry * direction;
    exit = origin + t_exit * direction;
end