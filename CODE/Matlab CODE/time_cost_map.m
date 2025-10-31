% 创建分组柱状图比较不同算法对不同网格大小的计算成本（对数坐标版）

% 创建图形窗口，设置位置和大小
figure('Position',[500,200,850,650])

YData1=[0.470000000000000	0.370000000000000	0.450000000000000	0.480000000000000	0.340000000000000
1.08000000000000	0.940000000000000	0.990000000000000	1.02000000000000	0.910000000000000
0.520000000000000	0.440000000000000	0.490000000000000	0.510000000000000	0.430000000000000
5.01000000000000	4.87000000000000	4.89000000000000	4.92000000000000	4.82000000000000
43.7800000000000	41.8500000000000	42.5800000000000	41.2900000000000	41.2200000000000
582.390000000000	423.220000000000	545.560000000000	551.600000000000	400.830000000000];
YData2 = [0.700000000000000	0.660000000000000	0.660000000000000	0.760000000000000	0.620000000000000
1.90000000000000	1.69000000000000	1.86000000000000	1.93000000000000	1.66000000000000
0.910000000000000	0.740000000000000	0.870000000000000	0.910000000000000	0.740000000000000
14.8000000000000	13.8500000000000	13.8200000000000	13.4900000000000	13.7000000000000
79.3200000000000	76.5200000000000	78.2100000000000	76.3000000000000	75.5800000000000
834.160000000000	784.250000000000	802.430000000000	774.140000000000	743.300000000000];

% 配色方案(RGB值转换为0-1范围)
CData=[111,173,72;    % 绿色 
       92,154,215;    % 蓝色 
      255,192,1;      % 黄色 
       69,103,42;     % 深绿 
       36,94,144]./255; % 深蓝 

% 创建上子图和柱状图
ax1=subplot(2,1,1); hold on; 
hBar1=bar(YData1);           

% 创建下子图和柱状图
ax2=subplot(2,1,2); hold on; 
hBar2=bar(YData2);           

% 设置柱状图样式
for i=1:length(hBar1)
    hBar1(i).EdgeColor='none';      
    hBar1(i).FaceColor=CData(i,:);  
    hBar2(i).EdgeColor='none';      
    hBar2(i).FaceColor=CData(i,:);  
end

% 设置Y轴为对数坐标 - 这是解决大范围数据的关键
set(ax1, 'YScale', 'log');
set(ax2, 'YScale', 'log');

% 设置X轴标签
ax1.XTick=1:size(YData1,1);        
ax2.XTick=1:size(YData2,1);        
ax1.XTickLabel={'5*5*90','5*5*900','25*25*90','25*25*900','250*250*90','250*250*900'};
ax2.XTickLabel={'5*5*90','5*5*900','25*25*90','25*25*900','250*250*90','250*250*900'};

% 设置坐标轴字体样式
ax1.FontName='Cambria';    ax2.FontName='Cambria';     
ax1.FontWeight='bold';     ax2.FontWeight='bold';      
ax1.FontSize=12;           ax2.FontSize=12;            

% 网格设置
ax1.XGrid='on';            ax2.XGrid='on';             
% ax1.YGrid='on';            ax2.YGrid='on';  % 添加Y轴网格
ax1.GridAlpha=.2;          ax2.GridAlpha=.2;           

% 坐标轴边框设置
ax1.Box='on';              ax2.Box='on';               
ax1.LineWidth=1.5;         ax2.LineWidth=1.5;          

% 刻度线设置
ax1.TickLength=[0,0];      ax2.TickLength=[0,0];       

% 绘制垂直分隔线
XV=(1:size(YData1,1)-1)+.5; 
for i=1:length(XV)
    xline(ax1,XV(i),'LineWidth',1.4,'LineStyle','--','Color',[0,0,0]);
    xline(ax2,XV(i),'LineWidth',1.4,'LineStyle','--','Color',[0,0,0]);
end

% 创建图例
lgd1=legend(hBar1,'program A','program B','program C','program D','program E','FontSize',13);
lgd2=legend(hBar2,'program A','program B','program C','program D','program E','FontSize',13);

% 图例位置设置
lgd1.Location='southoutside';  
lgd2.Location='southoutside';  

% 图例排列设置
lgd1.NumColumns=length(hBar1); 
lgd2.NumColumns=length(hBar2); 

% 图例标记大小设置
lgd1.ItemTokenSize=[8,8];      
lgd2.ItemTokenSize=[8,8];      

% 图例边框设置
lgd1.Box='off';                                                             
lgd2.Box='off';                

% 调整子图位置和间距
set(ax1,'LooseInset',[.1,0,0.028,0.03],'OuterPosition',[0,1/2-1/30,1,1/2+1/30]);
set(ax2,'LooseInset',[.1,0,0.028,0.03],'OuterPosition',[0,0-1/30,1,1/2+1/30]);

% 添加子图标签(左上角) - 
text(ax1,0.6, 1000,'5:00 UT','FontSize',15,'FontWeight','bold','FontName','Cambria');
text(ax2,0.6, 1000,'11:00 UT','FontSize',15,'FontWeight','bold','FontName','Cambria');

% 添加Y轴标签
ylabel(ax1, 'time cost(s)');
ylabel(ax2, 'time cost(s)');

% 设置合适的Y轴范围（确保所有数据点可见）
ylim(ax1, [0.1, 10000]); 
ylim(ax2, [0.1, 10000]);
