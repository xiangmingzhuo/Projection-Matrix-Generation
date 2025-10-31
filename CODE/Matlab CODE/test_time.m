for i = 1:5 
    tic;
    [A, sortedjd, distance] = Get_A1_canshu_simple1(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
    runTimes(i) = toc;
    clear A sortedjd distance;
end
averageTime(1,1) = min(runTimes);
clear runTimes;
for i = 1:5
    tic;
    [A, sortedjd, distance] = Get_A1_canshu_simple2(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
    runTimes(i) = toc;
    clear A sortedjd distance;
end
averageTime(1,2) = min(runTimes);
clear runTimes;
for i = 1:5 
    tic;
    [A, sortedjd, distance] = Get_A1_trangle_simple1(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length);
    runTimes(i) = toc;
    clear A sortedjd distance;
end
averageTime(1,3) = min(runTimes);
clear runTimes;
for i = 1:5 
    tic;
    [A, sortedjd, distance] = Get_A1_trangle_simple2(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length);
    runTimes(i) = toc;
    clear A sortedjd distance;
end
averageTime(1,4) = min(runTimes);
clear runTimes;
for i = 1:5
    tic;
    [A1, sortedjd1, distance1] = Get_A1_trangle_simple3(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length);
    runTimes(i) = toc;
    clear A sortedjd distance;
end
averageTime(1,5) = min(runTimes);
clear runTimes;
save averageTime;
% for i = 1:3 
%     tic;
%     [A, sortedjd, distance] = Get_A1_canshu_sivt_simple1(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
%     runTimes(i) = toc;
%     clear A sortedjd distance;
% end
% averageTime(2,1) = min(runTimes);
% clear runTimes;
% for i = 1:3
%     tic;
%     [A1, sortedjd, distance] = Get_A1_canshu_sivt_simple2(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,total_length);
%     runTimes(i) = toc;
%     clear A sortedjd distance;
% end
% averageTime(2,2) = min(runTimes);
% clear runTimes;
% for i = 1:3 
%     tic;
%     [A, sortedjd, distance] = Get_A1_trangle_sivt_simple1(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length);
%     runTimes(i) = toc;
%     clear A sortedjd distance;
% end
% averageTime(2,3) = min(runTimes);
% clear runTimes;
% for i = 1:3
%     tic;
%     [A, sortedjd, distance] = Get_A1_trangle_sivt_simple2(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length);
%     runTimes(i) = toc;
%     clear A sortedjd distance;
% end
% averageTime(2,4) = min(runTimes);
% clear runTimes;
% for i = 1:3
%     tic;
%     [A, sortedjd, distance] = Get_A1_trangle_sivt_simple3(station, satellite, jdmin, jdmax, jdjg, wdmin, wdmax, wdjg, gdmin, gdmax, gdjg,jdmin_rad,wdmin_rad,jdmax_rad,wdmax_rad,jdjg_rad,wdjg_rad,total_length);
%     runTimes(i) = toc;
%     clear A sortedjd distance;
% end
% averageTime(2,5) = min(runTimes);
% clear runTimes;