% 辅助函数：检查值是否在范围内
function in_bounds = check_in_bounds(value, min_val, max_val)
    in_bounds = (value >= min_val) && (value <= max_val);
end