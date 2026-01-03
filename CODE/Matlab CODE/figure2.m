load time_result.mat;
% 创建图形窗口
f2 = figure('Position', [100, 100, 800, 600]); % 设置图形大小

% 绘制散点图
h1 = scatter(result(:,2), result(:,3), 80, 'filled', 'MarkerFaceColor', [0.2, 0.6, 0.8], 'MarkerEdgeColor', 'none');
hold on;
h2 = scatter(result(:,2), result(:,4), 80, 'filled', 'MarkerFaceColor', [0.8, 0.4, 0.2], 'MarkerEdgeColor', 'none');

% 设置坐标轴属性
xlabel('Total number of voxels', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Time Cost (seconds)', 'FontSize', 18, 'FontWeight', 'bold');

% 设置坐标轴范围和刻度
x_min = min(result(:,2));
x_max = max(result(:,2));
y_min = min(min(result(:,3)), min(result(:,4)));
y_max = max(max(result(:,3)), max(result(:,4)));

xlim([x_min, x_max]);
ylim([y_min * 0.95, y_max * 1.05]); % 留出一些边距

% 网格和背景
grid on;
grid minor;
set(gca, 'GridAlpha', 0.3);
set(gca, 'MinorGridAlpha', 0.1);
set(gca, 'FontSize', 16);

% 如果数据点较多，添加对数拟合趋势线
if size(result,1) > 2
    % 对数拟合：y = a * log(x) + b
    % 方法1：直接对数拟合
    log_fit1 = fit(log(result(:,2)), result(:,3), 'poly1');
    log_fit2 = fit(log(result(:,2)), result(:,4), 'poly1');
    
    % 方法2：使用多项式拟合对数关系（更稳定的方法）
    % 创建对数坐标
    log_x = log(result(:,2));
    
    % 确保没有负数或零（对数需要正数）
    if all(result(:,2) > 0)
        % 使用线性拟合对数关系：y = p(1)*log(x) + p(2)
        p1_log = polyfit(log_x, result(:,3), 1);
        p2_log = polyfit(log_x, result(:,4), 1);
        
        % 创建拟合曲线数据（对数空间）
        x_fit = linspace(x_min, x_max, 100);
        y_fit1 = p1_log(1) * log(x_fit) + p1_log(2);
        y_fit2 = p2_log(1) * log(x_fit) + p2_log(2);
        
        % 绘制趋势线并创建图例句柄
        h3 = plot(x_fit, y_fit1, '--', 'Color', [0.2, 0.6, 0.8], 'LineWidth', 2);
        h4 = plot(x_fit, y_fit2, '--', 'Color', [0.8, 0.4, 0.2], 'LineWidth', 2);
        
        % 使用归一化坐标精确控制图例位置
        legend([h1, h2, h3, h4], ...
               'SIVT', 'Spherical Trigonometry', ...
               'SIVT Log Fit', 'Spherical Trigonometry Log Fit', ...
               'Location', 'northwest', 'FontSize', 16);
    else
        % 如果数据包含非正数，使用警告并跳过拟合
        warning('数据包含非正数，无法进行对数拟合');
        legend([h1, h2], 'SIVT', 'Spherical Trigonometry', ...
               'Location', 'northwest', 'FontSize', 16);
    end
else
    % 如果没有足够的数据点，使用原始图例
    legend([h1, h2], 'SIVT', 'Spherical Trigonometry', ...
           'Location', 'northwest', 'FontSize', 16);
end

% 设置图形框线
box on;

% 设置图形背景
set(gcf, 'Color', 'white');
set(gca, 'Color', [0.98, 0.98, 0.98]);

hold off;
exportgraphics(f2,'figure2.png','Resolution',600);