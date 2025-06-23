%日期
doy='15350';
load(['20',doy,'\Sites_Info15.mat']);
load(['SP3\','20',doy,'sp3.mat']);
%区域
jdmin=0;jdmax=25;jdjg=5;
wdmin=35;wdmax=60;wdjg=5;
gdmin=100000;gdmax=1000000;gdjg=100000;
%时间
t = 11;
% 获取 CPU 核心总数
numCores = feature('numcores');  
desiredCores = max(1, floor(1 * numCores));  
% 如果并行池未启动，则按计算的核心数启动
if isempty(gcp('nocreate'))
    parpool('local', desiredCores);
end
numRuns = 5;
runTimes = zeros(numRuns, 1);
%参数方程（隐式求解)
for i = 1:1
    tic;
    [A1,L1,sortedjd1,distance1]=Get_A1_initial2(doy,Sites_Info,sate,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,t);
    runTimes(i) = toc;
    clear A1 L1 sortedjd1 distance1;
end
averageTime1 = runTimes(i);
%参数方程(有效分裂面隐式求解)
for i = 1:1
    tic;
    [A2,L2,sortedjd2,distance2]=Get_A1_initial1(doy,Sites_Info,sate,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,t);
    runTimes(i) = toc;
    clear A2 L2 sortedjd2 distance2;
end
averageTime2 = runTimes(i);
%参数方程
numRuns = 5;
runTimes = zeros(numRuns, 1);
for i = 1:numRuns
    tic;
    [A3,L3,sortedjd3,distance3]=Get_A1_canshu2(doy,Sites_Info,sate,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,t);
    runTimes(i) = toc;
    clear A3 L3 sortedjd3 distance3;
end
averageTime3 = mean(runTimes);
%参数方程（有效分裂面计算）
numRuns = 5;
runTimes = zeros(numRuns, 1);
for i = 1:numRuns
    tic;
    [A4,L4,sortedjd4,distance4]=Get_A1_canshu1(doy,Sites_Info,sate,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,t);
    runTimes(i) = toc;
    clear A4 L4 sortedjd4 distance4;
end
averageTime4 = mean(runTimes);
%直接球面三角
numRuns = 5;
runTimes = zeros(numRuns, 1);
for i = 1:numRuns
    tic;
    [A5,L5,sortedjd5,distance5]=Get_A1(doy,Sites_Info,sate,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,t);
    runTimes(i) = toc;
    clear A5 L5 sortedjd5 distance5;
end
averageTime5 = mean(runTimes);
%二分法+球面三角（有效分裂面计算）
numRuns = 5;
runTimes = zeros(numRuns, 1);
for i = 1:numRuns
    tic;
    [A6,L6,sortedjd6,distance6]=Get_A1_Optimized(doy,Sites_Info,sate,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,t);
    runTimes(i) = toc;
    clear A6 L6 sortedjd6 distance6;
end
averageTime6 = mean(runTimes);
%球面三角+球面三角（有效分裂面计算）
numRuns = 5;
runTimes = zeros(numRuns, 1);
for i = 1:numRuns
    tic;
    [A7,L7,sortedjd7,distance7]=Get_A1_new(doy,Sites_Info,sate,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,t);
    runTimes(i) = toc;
    clear A7 L7 sortedjd7 distance7;
end
averageTime7 = mean(runTimes);