function [A,L,sortedjd,distance]=Get_A1_new(doy,Sites_Info,sate,jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg,t)
nj=(jdmax-jdmin)/jdjg;
nw=(wdmax-wdmin)/wdjg;
ng=(gdmax-gdmin)/gdjg;
nx=nj*nw*ng;
R = 6371000;
RDCB_ION = load(['20',doy,'\RDCB_ION.mat']);
RDCB_ION = RDCB_ION.RDCB_ION;
SDCB_REF = load(['20',doy,'\SDCB_REF.mat']);
SDCB_REF = SDCB_REF.SDCB_REF;
d=str2num(doy(3:end));
Sj=SDCB_REF.Gvalue(d,:);
Coor=Sites_Info.coor;
stations=Sites_Info.name;
r_name=unique(Sites_Info.name);
x=sate.x; y=sate.y; z=sate.z;

path_G2=['M_P4/' doy '/' 'G2'];
list_G2obs=dir([path_G2 '/*.mat']);
n_rg2=length(list_G2obs);

% 预分配并行变量
A_cell = cell(n_rg2,1);
L_cell = cell(n_rg2,1);
sortedjd_cell = cell(n_rg2,1);
distance_cell = cell(n_rg2,1);

% 频率参数
f1 = 1575.42*ones(1,32);
f2 = 1227.6*ones(1,32);
Af1=f1.^2.*f2.^2./(402800.*(f1.^2-f2.^2));


tic
parfor i=1:n_rg2
    % 初始化当前测站的变量
    site = list_G2obs(i).name(1:4);
    data = load(fullfile(path_G2, list_G2obs(i).name), 'P4G');
    P4 = data.P4G;
    
    % 查找索引
    r_g2 = find(strcmpi(site,r_name),1);
    index = find(strcmpi(site,stations),1);
    Rj = RDCB_ION(d,r_g2);
    
    if Rj==0
        A_cell{i} = [];
        L_cell{i} = [];
        continue;
    end
    
    sx=Coor(index,1); sy=Coor(index,2); sz=Coor(index,3);		
    [lat,lon,~] = XYZtoBLH_sphere(sx,sy,sz,R);
    
    % 区域检查
    if lat<deg2rad(wdmin) || lat>deg2rad(wdmax) || lon<deg2rad(jdmin) || lon>deg2rad(jdmax)
        A_cell{i} = [];
        L_cell{i} = [];
        continue;
    end
    
    k0=(t-1)*120+91;
    k2=t*120;
    
    % 预分配当前测站的结果
    max_possible = (k2-k0+1)*32;
    A_local = sparse(max_possible, nx);
    L_local = sparse(max_possible, 1);
    sortedjd_local = cell(max_possible,1);
    distance_local = cell(max_possible,1);
    valid_cnt = 0;
    
    for k=k0:k2
        for j=1:32
            if P4(k,j)==0
			    A_cell{i} = [];
				L_cell{i} = [];
                continue;
            end
            
            % STEC计算
            L1_val = (P4(k,j)-(Rj+Sj(j))*299792458/1e9)*-Af1(j);
            if L1_val<0
			    A_cell{i} = [];
				L_cell{i} = [];
                continue;               
            end
            
            % 球面三角计算
            [A1,sorted_ddzb,distances] = Get_sta_A1_new_optimized(sx,sy,sz,...
                x(k,j)*1000,y(k,j)*1000,z(k,j)*1000,...
                jdmin,jdmax,jdjg,wdmin,wdmax,wdjg,gdmin,gdmax,gdjg);
            
            if isempty(A1)
			    A_cell{i} = [];
				L_cell{i} = [];
                continue;
            end
            
            % if max(sorted_ddzb(3,:))<gdmax-gdjg
            %     continue;
            % end
            
            % 存储有效数据
            valid_cnt = valid_cnt + 1;
            A_local(valid_cnt,:) = A1;
            L_local(valid_cnt) = L1_val;
            sortedjd_local{valid_cnt} = sorted_ddzb';
            distance_local{valid_cnt} = distances';
        end
    end
    
    % 裁剪到实际有效数据
    A_cell{i} = A_local(1:valid_cnt,:);
    L_cell{i} = L_local(1:valid_cnt);
    sortedjd_cell{i} = sortedjd_local(1:valid_cnt);
    distance_cell{i} = distance_local(1:valid_cnt);
end
toc
% 合并所有测站的结果
A = vertcat(A_cell{:});
L = vertcat(L_cell{:});

% 合并结构体数据
% sortedjd = struct();
% distance = struct();
% global_idx = 1;
% for i=1:n_rg2
%     for j=1:length(sortedjd_cell{i})
%         sortedjd.(['shexian',num2str(global_idx)]) = sortedjd_cell{i}{j};
%         distance.(['shexian',num2str(global_idx)]) = distance_cell{i}{j};
%         global_idx = global_idx + 1;
%     end
% end
% 先计算总元素数量
total_elements = sum(cellfun(@numel, sortedjd_cell));
% 预分配单元格数组
all_sortedjd = cell(1, total_elements);
all_distance = cell(1, total_elements);
% 填充数据
idx = 1;
for i = 1:numel(sortedjd_cell)
    current_length = numel(sortedjd_cell{i});
    all_sortedjd(idx:idx+current_length-1) = sortedjd_cell{i}(:);
    all_distance(idx:idx+current_length-1) = distance_cell{i}(:);
    idx = idx + current_length;
end
% 生成字段名
field_names = arrayfun(@(x) ['shexian', num2str(x)], 1:total_elements, 'UniformOutput', false);
% 创建结构体
sortedjd = cell2struct(all_sortedjd, field_names, 2);
distance = cell2struct(all_distance, field_names, 2);
end