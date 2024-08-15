%%
clc
clear
a1 = xlsread('EDPs.xlsx',2);
PID_x = [a1(:,1) a1(:,4:15)]';
a2 = xlsread('EDPs.xlsx',3);
PID_y = [a2(:,1) a2(:,4:15)]';

a3 = xlsread('EDPs.xlsx',4);
PFA_x = [a3(:,1) a3(:,4:15)]';
a4 = xlsread('EDPs.xlsx',5);
PFA_y = [a4(:,1) a4(:,4:15)]';

a5 = xlsread('EDPs.xlsx',6);
LBR_x = [a5(:,1) a5(:,4:15)]';
a6 = xlsread('EDPs.xlsx',7);
LBR_y = [a6(:,1) a6(:,4:15)]';

a7 = xlsread('EDPs.xlsx',8);
EDR_x = [a7(:,1) a7(:,4:15)]';
a8 = xlsread('EDPs.xlsx',9);
EDR_y = [a8(:,1) a8(:,4:15)]';


Index = {};
n = 1;

% PID
for i = 1:length(PID_x)
    Index{n} = strcat('PID-',num2str(PID_x(1,i)),'-1' );
    Index{n+1} = strcat('PID-',num2str(PID_y(1,i)),'-2' );
        
    Demand(:,n)= PID_x(2:end,i);
    Demand(:,n+1) = PID_y(2:end,i);
    
    n = n+2;
end

% PFA
for i = 1:length(PFA_x)
    Index{n} = strcat('PFA-',num2str(PFA_x(1,i)),'-1' );
    Index{n+1} = strcat('PFA-',num2str(PFA_y(1,i)),'-2' );
    
    Demand(:,n)= PFA_x(2:end,i);
    Demand(:,n+1) = PFA_y(2:end,i);
    
    n = n+2;
end

% LBR
for i = 1:length(LBR_x)
    Index{n} = strcat('LBR-',num2str(LBR_x(1,i)),'-1' );
    Index{n+1} = strcat('LBR-',num2str(LBR_x(1,i)),'-2' );
        
    Demand(:,n)= LBR_x(2:end,i);
    Demand(:,n+1) = LBR_x(2:end,i);
    
    n = n+2;
end

% EDR
for i = 1:length(EDR_x)
    Index{n} = strcat('EDR-',num2str(EDR_x(1,i)),'-1' );
    Index{n+1} = strcat('EDR-',num2str(EDR_y(1,i)),'-2' );
        
    Demand(:,n)= EDR_x(2:end,i);
    Demand(:,n+1) = EDR_y(2:end,i);
    
    n = n+2;
end

%%
Num = (1:size(PID_x,1)-1)';
xlswrite('Demand-Sample.xlsx',Num,1,'A2')
xlswrite('Demand-Sample.xlsx',Index,1,'B1')
xlswrite('Demand-Sample.xlsx',Demand,1,'B2')


%% Example Matlab Code:
% Develop underlying statistics of the response history analysis
% clear;
num_realization = 2000;
% loading information: EDP, epistemic variability,
% EDPs =load('EDP.txt');
% EDPs = xlsread(strcat(aaa(i).name,'/R-MC.xls'));
EDPs = Demand;
% B=load(‘beta.txt’);%matrix for dispersions m ? and gm ? ,
B = [0.38;0.0];
B = repmat(B(:,1),1,size(EDPs,2));

% taking natural logarithm of the EDPs. Calling it lnEDPs
lnEDPs=log(EDPs); % Table G-2, or Table G-9
[num_rec num_var]=size(lnEDPs);
% finding the mean matrix of lnEDPs. Calling it lnEDPs_mean
lnEDPs_mean=mean(lnEDPs); % last row in Table G-2, or Table G-9
lnEDPs_mean=lnEDPs_mean';
% finding the covariance matrix of lnEDPs. Calling it lnEDPs_cov
lnEDPs_cov=cov(lnEDPs); % Table G-3, or Table G-10
% finding the rank of covariance matrix of lnEDPs. Calling it
% lnEDPs_cov_rank
lnEDPs_cov_rank=rank(lnEDPs_cov);
% inflating the variances with epistemic variability
sigma = sqrt(diag(lnEDPs_cov)); % sqrt first to avoid under/overflow
sigmap2 = sigma.*sigma;
R = lnEDPs_cov ./ (sigma*sigma');
B=B';
sigmap2=sigmap2+(B(:,1).* B(:,1)); % Inflating variance for m ?
sigmap2=sigmap2+(B(:,2).* B(:,2)); % Inflating variance for gm ?
sigma=sqrt(sigmap2);
sigma2=sigma*sigma';
lnEDPs_cov_inflated=R.*sigma2;
% finding the eigenvalues eigenvectors of the covariance matrix.Calling them D2_total and L_total
[L_total D2_total]=eig(lnEDPs_cov_inflated); % Table G-5, Table G-13
D2_total=eig(lnEDPs_cov_inflated);
% Partition L_total to L_use. L_use is the part of eigenvector matrix
% L_total that corresponds to positive eigenvalues
if lnEDPs_cov_rank >= num_var
    L_use =L_total; % Table G-5
else
    L_use =(L_total(:,num_var- lnEDPs_cov_rank+1:num_var));
    %Table G-13
end
% Partition the D2_total to D2_use. D2_use is the part of eigenvalue
%vector D2_total that corresponds to positive eigenvalues
if lnEDPs_cov_rank >= num_var
    D2_use =D2_total;
else
    D2_use =D2_total(num_var- lnEDPs_cov_rank+1:num_var);
end
% Find the square root of D2_use and call is D_use.
D_use =diag((D2_use).^0.5); %Table G-4, or Table G-12
% Generate Standard random numbers
if lnEDPs_cov_rank >= num_var
    U = randn(num_realization,num_var) ;
else
    U = randn(num_realization, lnEDPs_cov_rank) ;
end
U = U' ;
% Create Lambda = D_use . L_use
Lambda = L_use * D_use ;
% Create realizations matrix
Z = Lambda * U + lnEDPs_mean * ones(1,num_realization) ;
lnEDPs_sim_mean=mean(Z');
lnEDPs_sim_cov=cov(Z');
A=lnEDPs_sim_mean./lnEDPs_mean'; %Table G-7, or Table G-16
B=lnEDPs_sim_cov./lnEDPs_cov; %Table G-8, or Table G-17
W=exp(Z);
Demand_M = W';

%% 扩充后数据
Num_M = (1:num_realization)';
xlswrite('Demand-Sample-M.xlsx',Num_M,1,'A2')
xlswrite('Demand-Sample-M.xlsx',Index,1,'B1')
xlswrite('Demand-Sample-M.xlsx',Demand_M,1,'B2')

%% 扩充前后数据对比
DS = readmatrix('Demand-Sample.xlsx','Range', 'B2');
DS_M = readmatrix('Demand-Sample-M.xlsx','Range', 'B2');

DS_median = median(DS);
DS_M_median = median(DS_M);

PID_median = DS_median(1:2:122);
PID_M_median = DS_M_median(1:2:122);

floor = 1:61;
plot(PID_median*100,floor)
hold on
plot(PID_M_median*100,floor)

PID_M_median'

%% 数据整理
clc
clear
PID = []; PFA = []; LBR = []; EDR = [];
%%
a1 = xlsread('EDPs.xlsx',2);
PID_x = [a1(:,1) a1(:,16)];
a2 = xlsread('EDPs.xlsx',3);
PID_y = [a2(:,1) a2(:,16)];

a3 = xlsread('EDPs.xlsx',4);
PFA_x = [a3(:,1) a3(:,16)];
a4 = xlsread('EDPs.xlsx',5);
PFA_y = [a4(:,1) a4(:,16)];

a5 = xlsread('EDPs.xlsx',6);
LBR_x = [a5(:,1) a5(:,16)];
a6 = xlsread('EDPs.xlsx',7);
LBR_y = [a6(:,1) a6(:,16)];

a7 = xlsread('EDPs.xlsx',8);
EDR_x = [a7(:,1) a7(:,16)];
a8 = xlsread('EDPs.xlsx',9);
EDR_y = [a8(:,1) a8(:,16)];

%%
PID = [PID PID_x PID_y];
PFA = [PFA PFA_x PFA_y];
LBR = [LBR LBR_x LBR_y];
EDR = [EDR EDR_x EDR_y];

%%
save('EDP_all.mat','PID','PFA','LBR','EDR')

%%
clc
clear
load('EDP_all.mat')

PID_max = max(LBR(:,2:2:24));
% 计算分组数目
numGroups = floor(length(PID_max) / 2);
% 初始化结果向量
result = zeros(1, numGroups);

% 对每组进行操作
for i = 1:numGroups
    % 计算当前组的索引
    index1 = 2*i - 1;
    index2 = 2*i;
    
    % 选择当前组中的较大值
    result(i) = max(PID_max(index1), PID_max(index2));
end

result = result*100;

%% 倒塌对比 - 整理文件夹
clc
clear
File = dir('B*-*')
save('File.mat','File')

%% 倒塌对比 - 整理数据
clc
clear
load('File.mat');
Result = [];
Cost_m = [];

for k = 1:length(File)
    cd(strcat(File(k).folder,'/',File(k).name))
    DL = readmatrix('DL_summary.csv');
    
    % 提取非零数值
    nonZeroValues = DL(DL(:, 2) ~= 0, 2);    
    % 计算非零数的中位值
    Result(k,1) = length(nonZeroValues);
    Result(k,2) = sum(DL(:,5));
    Result(k,3) = sum(DL(:,6));
    Cost_m(k,1) = median(nonZeroValues);
end

%% 恢复时间-总
clc
clear

load('File.mat');
for k = 1:length(File)
    cd(strcat(File(k).folder,'/',File(k).name,'/Output_treads_R'))
    DT_s = readmatrix('DT_stepfunc_FR.csv');
    DT = DT_s(3:end,end);
    % 将数据按升序排序
    sortedData = sort(DT(:, end));
    % 找到中间两个数的索引
    numElements = numel(sortedData);
    middleIndex1 = numElements / 2;
    middleIndex2 = middleIndex1 + 1;
    % 取索引较大的那个数作为中位值中的较大值
    largerValue = max(sortedData(middleIndex1), sortedData(middleIndex2));
    DT_m(k,1) = largerValue;
end

save('DT_m.mat','DT_m')

%% 恢复时间-分解 大类
clc
clear
i = 1

load('E:\OneDrive - tongji.edu.cn\博士后归档\程序代码\SAUSAGE 61-Story - R\DT_m.mat');
DT_s = readmatrix('DT_stepfunc_FR.csv');
nonZeroValues = DT_s(DT_s(:, 4) ~= 0, 4);
N = length(nonZeroValues);
DT = DT_s(3:N+2,end);

Impeding = readmatrix('IF_delays.csv');
Impeding_I = Impeding(1:N,2);
Impeding_EP = Impeding(1:N,3)+Impeding(1:N,4);
Impeding_F = Impeding(1:N,5);
Impeding_C = max(Impeding(1:N,6:12),[],2);

Utility_s = xlsread('DT_path_FR.xlsx',5);
Utility = Utility_s(3:end,end);

index = DT==Impeding_I;
Impeding_EP(index,:) = 0;
Impeding_F(index,:) = 0;
Impeding_C(index,:) = 0;

Utility_m = median(Utility);

Tr = 854;

Imp1 = Impeding_I+max([Impeding_F Impeding_EP Impeding_C],[],2);
Imp2 = DT_s(N+3:end,end) - Tr;
Imp = [Imp1; Imp2];
Impeding_m = median(Imp);

% Repair = DT-max([Utility Imp],[],2);
Repair1 = DT-max([Imp1],[],2);
Repair2  = repmat(Tr, length(Imp2), 1) ;
Repair = [Repair1; Repair2];
Repair_m = median(Repair);

Result = [Utility_m Impeding_m Repair_m]
median_value = DT_m(i);
% 找到中位值所在的位置
median_position = find(DT == median_value);
Result1 = [Utility(median_position) Imp(median_position) Repair(median_position)]
%% 恢复时间-分解 RS
clc
clear

RT_s = readmatrix('RT_RSeq_FR.csv'); %没有考虑多层同时修复
RS1 = sum(RT_s(:,2:7:end),2);
RS2 = sum(RT_s(:,3:7:end),2);
RS3 = sum(RT_s(:,4:7:end),2);
RS4 = sum(RT_s(:,5:7:end),2);
RS5 = sum(RT_s(:,6:7:end),2);
RS6 = sum(RT_s(:,7:7:end),2);
RS7 = sum(RT_s(:,8:7:end),2);
    
%%
clc
clear
a1 = xlsread('RT_stepfunc_FR.xlsx',1);
RS1 = a1(:,end);
a2 = xlsread('RT_stepfunc_FR.xlsx',2);
RS2 = a2(:,end);
a3 = xlsread('RT_stepfunc_FR.xlsx',3);
RS3 = a3(:,end);
a4 = xlsread('RT_stepfunc_FR.xlsx',4);
RS4 = a4(:,end);
a5 = xlsread('RT_stepfunc_FR.xlsx',5);
RS5 = a5(:,end);
a6 = xlsread('RT_stepfunc_FR.xlsx',6);
RS6 = a6(:,end);
a7 = xlsread('RT_stepfunc_FR.xlsx',7);
RS7 = a7(:,end);
    
Result = median([RS1 RS2 RS3 RS4 RS5 RS6 RS7]);
%% 
nonZeroValues = DT_s(DT_s(:, 4) ~= 0, 4);
N = length(nonZeroValues);
DT = DT_s(3:N+2,end);

Impeding = readmatrix('IF_delays.csv');
Impeding_I = Impeding(1:N,2);
Impeding_EP = Impeding(1:N,3)+Impeding(1:N,4);
Impeding_F = Impeding(1:N,5);
Impeding_C = max(Impeding(1:N,6:12),[],2);

Utility_s = xlsread('DT_path_FR.xlsx',5);
Utility = Utility_s(3:end,end);

index = DT==Impeding_I;
Impeding_EP(index,:) = 0;
Impeding_F(index,:) = 0;
Impeding_C(index,:) = 0;

Utility_m = median(Utility);

Tr = 854;

Imp1 = Impeding_I+max([Impeding_F Impeding_EP Impeding_C],[],2);
Imp2 = DT_s(N+3:end,end) - Tr;
Imp = [Imp1; Imp2];
Impeding_m = median(Imp);

% Repair = DT-max([Utility Imp],[],2);
Repair1 = DT-max([Imp1],[],2);
Repair2  = repmat(Tr, length(Imp2), 1) ;
Repair = [Repair1; Repair2];
Repair_m = median(Repair);

Result = [Utility_m Impeding_m Repair_m]

%% 恢复时间-分解 Utility起控制作用
clc
clear

DT_s = readmatrix('DT_stepfunc_FR.csv');
nonZeroValues = DT_s(DT_s(:, 4) ~= 0, 4);
N = length(nonZeroValues);
DT = DT_s(3:N+2,end);

Impeding = readmatrix('IF_delays.csv');
Impeding_I = Impeding(1:N,2);
Impeding_EP = Impeding(1:N,3)+Impeding(1:N,4);
Impeding_F = Impeding(1:N,5);
Impeding_C = max(Impeding(1:N,6:12),[],2);

Utility_s = xlsread('DT_path_FR.xlsx',5);
Utility = Utility_s(3:end,end);

U = DT-Utility(1:N,1);
nonZeroValuesU = U(U < 0, 1);
NU = length(nonZeroValuesU)/2000;

index = DT==Impeding_I;
Impeding_EP(index,:) = 0;
Impeding_F(index,:) = 0;
Impeding_C(index,:) = 0;

Utility_m = median(Utility);

Tr = 854;

Imp1 = Impeding_I+max([Impeding_F Impeding_EP Impeding_C],[],2);
Imp2 = DT_s(N+3:end,end) - Tr;
Imp = [Imp1; Imp2];
Impeding_m = median(Imp);

% Repair = DT-max([Utility Imp],[],2);
Repair1 = DT-max([Imp1],[],2);
Repair2  = repmat(Tr, length(Imp2), 1) ;
Repair = [Repair1; Repair2];
Repair_m = median(Repair);

Result = [Utility_m Impeding_m Repair_m]