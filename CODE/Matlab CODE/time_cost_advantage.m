%% Journal-style visualization (2 figures) - ALL FONTS BOLD
clear; clc; close all;
load time_cost_advantage_data.mat;
kk = 4;
time_cost4(:,5:6) = time_cost4(:,5:6)*100;
time_cost = time_cost4';
% -------------------- Data --------------------
gridLbl = { ...
    '5×5×9','5×5×90','5×5×900', ...
    '25×25×9','25×25×90','25×25×900', ...
    '250×250×9','250×250×90','250×250×900'};

param      = time_cost(1,:);
paramSplit = time_cost(2,:);

sph        = time_cost(3,:);
sphSplit   = time_cost(4,:);

advSph      = time_cost(5,:);
advSphSplit = time_cost(6,:);

x = 1:numel(gridLbl);

% -------------------- Style settings (journal + bold) --------------------
set(0,'DefaultFigureColor','w');

% Font family
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');

% Font size / axes
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultAxesLineWidth',1);
set(0,'DefaultLineLineWidth',1.8);

% Key: make ALL fonts bold by default
set(0,'DefaultAxesFontWeight','bold');
set(0,'DefaultTextFontWeight','bold');
set(0,'DefaultLegendFontSize',10);

% Marker styles
mk1 = 'o'; mk2 = 's';

% Helper: add value labels on bars (bold)
addBarLabels = @(xb,yb) arrayfun(@(i) text(xb(i), yb(i)+0.8, ...
    sprintf('%.2f', yb(i)), ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontSize',9, ...
    'FontWeight','bold'), 1:numel(yb));

% ==================== Figure 1: Spherical vs Parametric ====================
f1 = figure('Units','centimeters','Position',[3 3 18 16]); % journal-friendly size
tl1 = tiledlayout(f1,2,1,'TileSpacing','compact','Padding','compact');

% (a) Metric comparison
ax1 = nexttile(tl1,1);
plot(x, param, ['-' mk1], 'MarkerSize',6); hold on;
plot(x, sph,   ['-' mk2], 'MarkerSize',6);
grid on; box on;
ylabel('Metric value','FontWeight','bold');
title('Comparison: Parametric vs Spherical Triangle','FontWeight','bold');
lg1 = legend({'Parametric','Spherical Triangle'}, ...
    'Location','northwest','Box','off');
set(lg1,'FontWeight','bold');

set(ax1,'XTick',x,'XTickLabel',gridLbl, 'FontWeight','bold');
xtickangle(25);

% Optional: log scale if desired
% set(ax1,'YScale','log'); ylabel('Metric value (log scale)','FontWeight','bold');

% (b) Advantage bar
ax2 = nexttile(tl1,2);
bar(x, advSph, 0.65, 'FaceAlpha',0.9);
grid on; box on;
ylabel('Advantage (%)','FontWeight','bold');
xlabel('Grid size','FontWeight','bold');
title('Advantage of Spherical Triangle over Parametric','FontWeight','bold');
set(ax2,'XTick',x,'XTickLabel',gridLbl,'FontWeight','bold');
xtickangle(25);
ylim([0, max(advSph)+8]);
addBarLabels(x, advSph);

% Safety pass: force-bold any remaining text objects in the figure
set(findall(f1,'-property','FontWeight'),'FontWeight','bold');

% Export (high-res)
exportgraphics(f1,['advantage',num2str(kk*2-1),'.png'],'Resolution',600);
% exportgraphics(f1,'Fig1_Spherical_vs_Parametric.pdf','ContentType','vector');

% ================= Figure 2: With Effective Split =================
f2 = figure('Units','centimeters','Position',[3 3 18 16]);
tl2 = tiledlayout(f2,2,1,'TileSpacing','compact','Padding','compact');

% (a) Metric comparison with split
ax3 = nexttile(tl2,1);
plot(x, paramSplit, ['-' mk1], 'MarkerSize',6); hold on;
plot(x, sphSplit,   ['-' mk2], 'MarkerSize',6);
grid on; box on;
ylabel('Metric value','FontWeight','bold');
title('Comparison: Parametric+Split vs Spherical Triangle+Split','FontWeight','bold');
lg2 = legend({'Parametric + Effective split','Spherical Triangle + Effective split'}, ...
    'Location','northwest','Box','off');
set(lg2,'FontWeight','bold');

set(ax3,'XTick',x,'XTickLabel',gridLbl,'FontWeight','bold');
xtickangle(25);

% Optional: log scale if desired
% set(ax3,'YScale','log'); ylabel('Metric value (log scale)','FontWeight','bold');

% (b) Advantage bar with split
ax4 = nexttile(tl2,2);
bar(x, advSphSplit, 0.65, 'FaceAlpha',0.9);
grid on; box on;
ylabel('Advantage (%)','FontWeight','bold');
xlabel('Grid size','FontWeight','bold');
title('Advantage of Spherical Triangle+Split over Parametric+Split','FontWeight','bold');
set(ax4,'XTick',x,'XTickLabel',gridLbl,'FontWeight','bold');
xtickangle(25);
ylim([0, max(advSphSplit)+8]);
addBarLabels(x, advSphSplit);

% Safety pass: force-bold any remaining text objects in the figure
set(findall(f2,'-property','FontWeight'),'FontWeight','bold');

% % Export (high-res)
exportgraphics(f2,['advantage',num2str(kk*2),'.png'],'Resolution',600);
% exportgraphics(gcf, 'figure.png', 'ContentType', 'vector');
% exportgraphics(f2,'Fig2_SphericalSplit_vs_ParametricSplit.pdf','ContentType','vector');
% 
% disp('Done: Figures exported as high-resolution PNG + vector PDF (all fonts bold).');
