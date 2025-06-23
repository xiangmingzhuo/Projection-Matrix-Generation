% 通用二分搜索函数
function t = binary_search(f, t_low, t_high, epsilon, max_iters, find_entry)
if t_low < 0
    t_low = 0;
end
for i = 1:max_iters
    t_mid = (t_low + t_high) / 2;
    inside = f(t_mid);

    if find_entry
        if inside
            t_high = t_mid;
        else
            t_low = t_mid;
        end
    else
        if inside
            t_low = t_mid;
        else
            t_high = t_mid;
        end
    end

    if (t_high - t_low) < epsilon
        break;
    end
end

if find_entry
    t = t_high;
else
    t = t_low;
end
end


% function t = binary_search(f, t_low, t_high, epsilon, max_iters, find_entry)
% if t_low < 0
%     t_low = 0;
% end
% % 初始化
% best_t = NaN;
% is_inside = false;
% 
% for i = 1:max_iters
%     t_mid = (t_low + t_high) / 2;
%     current_inside = f(t_mid);
% 
%     % 更新最佳估计
%     if current_inside
%         best_t = t_mid;
%         is_inside = true;
%     end
% 
%     % 根据搜索类型调整区间
%     if find_entry
%         % 寻找进入点：当在外部时向右搜索，在内部时向左搜索
%         if current_inside
%             t_high = t_mid;
%         else
%             t_low = t_mid;
%         end
%     else
%         % 寻找退出点：当在内部时向右搜索，在外部时向左搜索
%         if current_inside
%             t_low = t_mid;
%         else
%             t_high = t_mid;
%         end
%     end
% 
%     % 检查收敛条件
%     if (t_high - t_low) < epsilon
%         break;
%     end
% end
% 
% % 返回最佳估计或边界值
% if find_entry
%     if is_inside
%         t = best_t;
%     else
%         % 如果没有找到进入点，返回t_high
%         t = t_high;
%     end
% else
%     if is_inside
%         t = best_t;
%     else
%         % 如果没有找到退出点，返回t_low
%         t = t_low;
%     end
% end
% 
% % 最终验证
% if isnan(t) || ~f(t)
%     % 如果没有有效交点，返回NaN
%     t = NaN;
% end
% end