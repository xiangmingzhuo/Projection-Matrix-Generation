%  二分法确定上界
function exit_point = binarySearchExitPoint(Bmin,Bmax, r0, d, t_first, t_far, epsilon, max_iters)
    % 输入参数:
    % G: 网格数据结构
    % r0: 射线起点 [3x1]
    % d: 射线方向（需归一化）[3x1]
    % t_first: 首个交点参数（进入点对应的t值）
    % t_far: 初始区间上限
    % epsilon: 容差
    % max_iters: 最大迭代次数
    %
    % 输出参数:
    % exit_point: 射线与网格的离开点 [3x1]
    if t_first < 0
        t_left = 0;
    else
        t_left = t_first;
    end
    t_right = t_far;
    
    for iter = 1:max_iters
        % 检查区间是否满足容差
        if (t_right - t_left) < epsilon
            break;
        end
        
        t_mid = (t_left + t_right) / 2;
        p = r0 + t_mid * d;  % 计算中点位置
        
        if isInsideMesh(p, Bmin,Bmax)
            t_left = t_mid;   % 交点在右半区间（继续向后搜索）
        else
            t_right = t_mid;  % 交点在左半区间（向前收缩）
        end
    end
    
    exit_point = r0 + t_left * d;  % 返回最后一个在网格内的点
end