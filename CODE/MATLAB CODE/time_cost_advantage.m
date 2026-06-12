%% Journal-style visualization (2 figures) - ALL FONTS BOLD
clear; clc; close all;
load time_cost_advantage_data.mat;
kk = 1;
time_cost1(:,5:6) = time_cost1(:,5:6)*100;
time_cost = time_cost1';

% -------------------- Data --------------------
gridLbl = { ...
    '1','2','3', ...
    '4','5','6', ...
    '7','8','9'};

param      = time_cost(1,:);
paramSplit = time_cost(2,:);

sph        = time_cost(3,:);
sphSplit   = time_cost(4,:);

advSph      = time_cost(5,:);
advSphSplit = time_cost(6,:);

x = 1:numel(gridLbl);

% -------------------- Style settings (journal + bold) --------------------
set(0,'DefaultFigureColor','w');

set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');

set(0,'DefaultAxesFontSize',14);
set(0,'DefaultAxesLineWidth',1);
set(0,'DefaultLineLineWidth',1.8);

set(0,'DefaultAxesFontWeight','bold');
set(0,'DefaultTextFontWeight','bold');
set(0,'DefaultLegendFontSize',14);

mk1 = 'o'; mk2 = 's';

addBarLabels = @(xb,yb) arrayfun(@(i) text(xb(i), yb(i)+0.8, ...
    sprintf('%.2f', yb(i)), ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontSize',12, ...
    'FontWeight','bold'), 1:numel(yb));

% ==================== Figure 1: Spherical vs Parametric ====================
% 关键改动：加宽 figure；tiledlayout 从 2×1 改为 1×2
f1 = figure('Units','centimeters','Position',[3 3 26 10]);  % 宽一些，适合横排
tl1 = tiledlayout(f1,1,2,'TileSpacing','compact','Padding','compact');

% (a) Metric comparison - 左边
ax1 = nexttile(tl1,1);
plot(x, param, ['-' mk1], 'MarkerSize',10); hold on;
plot(x, sph,   ['-' mk2], 'MarkerSize',10);
grid on; box on;
ylabel('Time Cost(s)','FontWeight','bold');
title('Comparison: A vs C','FontWeight','bold');
lg1 = legend({'A','C'}, ...
    'Location','northwest','Box','off');
set(lg1,'FontWeight','bold');
set(ax1,'XTick',x,'XTickLabel',gridLbl,'FontWeight','bold');
xtickangle(25);

% (b) Advantage bar - 右边
ax2 = nexttile(tl1,2);
bar(x, advSph, 0.65, 'FaceAlpha',0.9);
grid on; box on;
ylabel('Advantage (%)','FontWeight','bold');
title('Advantage of C overA','FontWeight','bold');
set(ax2,'XTick',x,'XTickLabel',gridLbl,'FontWeight','bold');
xtickangle(25);
ylim([0, max(advSph)+8]);
addBarLabels(x, advSph);

set(findall(f1,'-property','FontWeight'),'FontWeight','bold');
exportgraphics(f1,['advantage',num2str(kk*2-1),'.png'],'Resolution',600);

% ================= Figure 2: With Effective Split =================
% 同样改成横排 1×2
f2 = figure('Units','centimeters','Position',[3 3 26 10]);
tl2 = tiledlayout(f2,1,2,'TileSpacing','compact','Padding','compact');

% (a) Metric comparison with split - 左边
ax3 = nexttile(tl2,1);
plot(x, paramSplit, ['-' mk1], 'MarkerSize',6); hold on;
plot(x, sphSplit,   ['-' mk2], 'MarkerSize',6);
grid on; box on;
ylabel('Time Cost(s)','FontWeight','bold');
title('Comparison: B vs E','FontWeight','bold');
lg2 = legend({'B','E'}, ...
    'Location','northwest','Box','off');
set(lg2,'FontWeight','bold');
set(ax3,'XTick',x,'XTickLabel',gridLbl,'FontWeight','bold');
xtickangle(25);

% (b) Advantage bar with split - 右边
ax4 = nexttile(tl2,2);
bar(x, advSphSplit, 0.65, 'FaceAlpha',0.9);
grid on; box on;
ylabel('Advantage (%)','FontWeight','bold');
title('Advantage of E over B','FontWeight','bold');
set(ax4,'XTick',x,'XTickLabel',gridLbl,'FontWeight','bold');
xtickangle(25);
ylim([0, max(advSphSplit)+8]);
addBarLabels(x, advSphSplit);

set(findall(f2,'-property','FontWeight'),'FontWeight','bold');
exportgraphics(f2,['advantage',num2str(kk*2),'.png'],'Resolution',600);
