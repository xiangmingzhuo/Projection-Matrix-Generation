load time_result.mat;
% 创建图形窗口
f4 = figure('Position', [100, 100, 800, 600]);

% 您的数据
x_data = [18000, 36000, 72000, 144000, 288000, 576000, 1152000, 2304000, 4608000, 9216000];
y1_data = [0.32
0.48
0.55
0.91
1.00 
1.20 
1.66
1.88
2.15
3.01
]';
y2_data = [0.1
0.11
0.13
0.17
0.19
0.21
0.27
0.31
0.34
0.53
]';

% 绘制散点图
h1 = scatter(x_data, y1_data, 100, 'filled', ...
            'MarkerFaceColor', [0.2, 0.6, 0.8], ...
            'MarkerEdgeColor', 'white', ...
            'LineWidth', 1.5);
hold on;
h2 = scatter(x_data, y2_data, 100, 'filled', ...
            'MarkerFaceColor', [0.8, 0.4, 0.2], ...
            'MarkerEdgeColor', 'white', ...
            'LineWidth', 1.5);

% 设置对数坐标
set(gca, 'XScale', 'log');

% 设置横坐标刻度（10的幂次）
x_ticks = 10.^(4:7);  % 从10^4到10^7
set(gca, 'XTick', x_ticks);
set(gca, 'XTickLabel', {'10^4', '10^5', '10^6', '10^7'});

% 设置坐标轴标签
xlabel('Total number of voxels', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Time Cost (seconds)', 'FontSize', 18, 'FontWeight', 'bold');

% 设置坐标轴范围
xlim([10000, 1e7]);
ylim([0, 4]);

% 添加趋势线（在对数坐标下的线性拟合）
log_x = log10(x_data);

% SIVT趋势线（二次多项式拟合，更好地拟合曲线）
p1 = polyfit(log_x, y1_data, 2);
x_fit_log = linspace(min(log_x), max(log_x), 100);
x_fit = 10.^x_fit_log;
y_fit1 = polyval(p1, x_fit_log);
h3 = plot(x_fit, y_fit1, '--', 'Color', [0.2, 0.6, 0.8], 'LineWidth', 2.5);

% Spherical Trigonometry趋势线（线性拟合）
p2 = polyfit(log_x, y2_data, 1);
y_fit2 = polyval(p2, x_fit_log);
h4 = plot(x_fit, y_fit2, '--', 'Color', [0.8, 0.4, 0.2], 'LineWidth', 2.5);

% 网格和背景
grid on;
grid minor;
set(gca, 'GridAlpha', 0.4);
set(gca, 'MinorGridAlpha', 0.2);
set(gca, 'FontSize', 16);

% 图例
legend([h1, h2, h3, h4], ...
       'SIVT', 'Spherical Trigonometry', ...
       'SIVT Trend', 'Spherical Trigonometry Trend', ...
       'Location', 'northwest', 'FontSize', 16);


% 设置图形样式
box on;
set(gcf, 'Color', 'white');
set(gca, 'Color', [0.98, 0.98, 0.98]);


hold off;
exportgraphics(f4,'figure4.png','Resolution',600);