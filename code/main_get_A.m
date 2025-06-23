%此函数采用球面三角的方法实现批量获取A和STEC
%日期
doy='15350';
load(['20',doy,'\Sites_Info15.mat']);
load(['SP3\','20',doy,'sp3.mat']);
%区域
jdmin=0;jdmax=25;jdjg=5;
wdmin=35;wdmax=60;wdjg=5;
gdmin=100000;gdmax=1000000;gdjg=100000;
% 获取 CPU 核心总数
numCores = feature('numcores');  
desiredCores = max(1, floor(1 * numCores));  
% 如果并行池未启动，则按计算的核心数启动
if isempty(gcp('nocreate'))
    parpool('local', desiredCores);
end
% 不同方法获取指定时间下的投影矩阵
for t = 1:23
%     %参数方程（隐式求解)
%     [A1,L1,sortedjd1,distance1]=Get_A1_initial2(doy,Sites_Info,sate,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,t);
%     %参数方程(有效分裂面隐式求解)
%     [A2,L2,sortedjd2,distance2]=Get_A1_initial1(doy,Sites_Info,sate,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,t);
%     %参数方程
%     [A3,L3,sortedjd3,distance3]=Get_A1_canshu2(doy,Sites_Info,sate,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,t);
%     %参数方程（有效分裂面计算）
%     [A4,L4,sortedjd4,distance4]=Get_A1_canshu1(doy,Sites_Info,sate,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,t);
%     %直接球面三角
%     [A5,L5,sortedjd5,distance5]=Get_A1(doy,Sites_Info,sate,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,t);
%     %二分法+球面三角（有效分裂面计算）
%     [A6,L6,sortedjd6,distance6]=Get_A1_Optimized(doy,Sites_Info,sate,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,t);
    %球面三角+球面三角（有效分裂面计算）
    [A7,L7,sortedjd7,distance7]=Get_A1_new(doy,Sites_Info,sate,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,t);
    
    % 指定保存路径(这里以第七种方案为例)
    fname = sprintf([num2str(jdjg),'_',num2str(wdjg),'_',num2str(gdjg),'_',num2str(t)]);
    % 定义子文件夹名称
    subfolder = fullfile('Result',doy);
    % 检查子文件夹是否存在，如果不存在则创建
    if ~exist(subfolder, 'dir')
        mkdir(subfolder);
    end
    % 构建完整文件路径
    filename = fullfile(subfolder, fname);
    % 保存数据
    save(filename, 'A7','L7','sortedjd7','distance7');
end