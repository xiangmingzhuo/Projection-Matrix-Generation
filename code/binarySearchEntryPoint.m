%  二分法确定下界
function intersection_point = binarySearchEntryPoint(Bmin,Bmax, r0, d, t_near, t_far, epsilon, max_iters)
    % 输入参数:
    % G: 网格数据结构
    % r0: 射线起点 [3x1]
    % d: 射线方向（需归一化）[3x1]
    % t_near, t_far: 初始区间参数
    % epsilon: 容差
    % max_iters: 最大迭代次数
    %
    % 输出参数:
    % intersection_point: 射线与网格的进入点 [3x1]
    if t_near < 0
        t_left = 0 ;
    else
        t_left = t_near;
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
            t_right = t_mid;  % 交点在左半区间
        else
            t_left = t_mid;   % 交点在右半区间
        end
    end
    
    intersection_point = r0 + t_right * d;
end

