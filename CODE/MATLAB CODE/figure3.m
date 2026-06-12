load time_result.mat;
% 创建图形窗口
f3 = figure('Position', [100, 100, 800, 600]); % 设置图形大小

% 绘制散点图
h1 = scatter(result(:,1), result(:,4), 80, 'filled', 'MarkerFaceColor', [0.2, 0.6, 0.8], 'MarkerEdgeColor', 'none');
hold on;
h2 = scatter(result(:,1), result(:,5), 80, 'filled', 'MarkerFaceColor', [0.8, 0.4, 0.2], 'MarkerEdgeColor', 'none');

% 设置坐标轴属性
xlabel('Total number of cleavage planes', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Time Cost (seconds)', 'FontSize', 18, 'FontWeight', 'bold');

% 设置坐标轴范围和刻度
x_min = min(result(:,1));
x_max = max(result(:,1));
y_min = min(min(result(:,4)), min(result(:,5)));
y_max = max(max(result(:,4)), max(result(:,5)));

xlim([x_min, x_max]);
ylim([y_min * 0.95, y_max * 1.05]); % 留出一些边距

% 网格和背景
grid on;
grid minor;
set(gca, 'GridAlpha', 0.3);
set(gca, 'MinorGridAlpha', 0.1);
set(gca, 'FontSize', 16);

% 如果数据点较多，添加趋势线
if size(result,1) > 2
    p1 = polyfit(result(:,1), result(:,4), 1);
    p2 = polyfit(result(:,1), result(:,5), 1);
    
    x_fit = linspace(x_min, x_max, 100);
    y_fit1 = polyval(p1, x_fit);
    y_fit2 = polyval(p2, x_fit);
    
    % 绘制趋势线并创建图例句柄
    h3 = plot(x_fit, y_fit1, '--', 'Color', [0.2, 0.6, 0.8], 'LineWidth', 2);
    h4 = plot(x_fit, y_fit2, '--', 'Color', [0.8, 0.4, 0.2], 'LineWidth', 2);
    
    % 使用归一化坐标精确控制图例位置
    legend([h1, h2, h3, h4], ...
           'SIVT', 'Spherical Trigonometry', ...
           'SIVT Trend', 'Spherical Trigonometry Trend', ...
           'Location', 'northwest', 'FontSize', 16);
    % 格式：[left, bottom, width, height]（归一化坐标，0-1之间）
else
    % 如果没有趋势线，使用原始图例
    legend([h1, h2], 'SIVT', 'Spherical Trigonometry', ...
           'Location', 'best', 'FontSize', 16);
end


% 设置图形框线
box on;

% 设置图形背景
set(gcf, 'Color', 'white');
set(gca, 'Color', [0.98, 0.98, 0.98]);

hold off;
exportgraphics(f3,'figure3.png','Resolution',600);