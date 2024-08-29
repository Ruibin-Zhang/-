%% ά�޷�����ֵ
clc
clear

currentDir = pwd; % ��ȡ��ǰĿ¼·��
subfolder = dir('*-*');

cost_median = [];
for f = 1:length(subfolder)
    
    cd(subfolder(f).name)
    DL = readmatrix('Output_pelicun\DL_summary.csv');
    
    cost_median = [cost_median; median(DL(:,2))];
    cd(currentDir)
end

%% ά�޷���84%��֤�����ֵ
clc
clear

currentDir = pwd; % ��ȡ��ǰĿ¼·��
subfolder = dir('*-MCE');

cost_84 = [];
for f = 1:length(subfolder)
    
    cd(subfolder(f).name)
    DL = readmatrix('Output_pelicun\DL_summary.csv');
    data = DL(:,2);
    % ������̬�ֲ����
    pd = fitdist(data, 'Lognormal');
    
    % ���� 84% ��֤�ʵ����ֵ
    p = 0.84; % 84% �ٷ�λ
    value_at_84_percent = icdf(pd, p);
    cost_84 = [cost_84; value_at_84_percent];
    cd(currentDir)
end

%% 84%��λ��
clc
clear

currentDir = pwd; % ��ȡ��ǰĿ¼·��
subfolder = dir('*-MCE');

cost_84 = [];
for f = 1:length(subfolder)
    
    cd(subfolder(f).name)
    DL = readmatrix('Output_pelicun\DL_summary.csv');
    data = DL(:,2);
    % ����� 84 �ٷ�λ��
    p = 84; % �ٷ�λ
    value_at_84_percentile = prctile(data, p);
    cost_84 = [cost_84; value_at_84_percentile];
    cd(currentDir)
end

    
%% ��ʧ�ֲ�-dmg
clc
clear

currentDir = pwd; % ��ȡ��ǰĿ¼·��
subfolder = dir('*-*');

for f = 1:length(subfolder)
    
    cd(subfolder(f).name)
    % ���� CSV �ļ�
    filename = 'Output_preprocessing\loss_sample_dmg.csv';  % CSV �ļ��� !!!��Ҫ�޸�
    delimiter = ',';  % �ָ��������ţ�
    dataTable = readtable(filename, 'Delimiter', delimiter);
    
    DL = readmatrix('Output_pelicun\DL_summary.csv');
    
    % ��ȡ�ı�����
    % ��ȡ���� "COST" ����
    % columnsWithCost = dataTable(:, startsWith(dataTable.Properties.VariableNames, 'COST'));
    % columnsWithCost = columnsWithCost(:,1:end-1);
    % dataTable.Properties.VariableDescriptions
    % ��ȡ x y ֵ
    x = dataTable.Properties.VariableNames(2:end);
    % ��ȡ y ֵ
    y = dataTable{:, 2:end};
    
    % % ��ȡ��Ԫ��Ĵ�С
    [numRows, numCols] = size(y);
    % % ����һ���뵥Ԫ���С��ͬ������
    % y_numericData = zeros(numRows, numCols);
    % % ������Ԫ�񣬲������Ԫ��ת��Ϊ����
    % for i = 1:numRows
    %     for j = 1:numCols
    %         y_numericData(i, j) = str2double(y{i, j});
    %     end
    % end
    total_cost = 373000000;
    y_repairable = [];
    for i = 1:numRows
        if DL(i,5)+DL(i,6) == 0
            y_repairable = [y_repairable; y(i,:)];
        else
            y_repairable = [y_repairable; repmat(total_cost/numCols, 1, 37)];
        end
    end
    
    plot_x = x';
    plot_y = y_repairable';
    plot_y_m = median(y_repairable)';
    
    output_folder = 'Output_plot';  % ָ�������ļ���Ŀ¼
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    output_file = fullfile(output_folder, 'dmg.mat');  % ��������·��
    save(output_file, 'DL', 'plot_x', 'plot_y', 'plot_y_m');
    
    cd(currentDir)
end

%% ��ͬ�����Թ���ά�޷��öԱ�
clc
clear

currentDir = pwd; % ��ȡ��ǰĿ¼·��
% subfolder = dir('Baseline*');
dd = 'dmg.mat';

aFOE = load(strcat(currentDir,'\Baseline building-FOE\Output_plot\',dd));
aDBE = load(strcat(currentDir,'\Baseline building-DBE\Output_plot\',dd));
aMCE = load(strcat(currentDir,'\Baseline building-MCE\Output_plot\',dd));

plot_median_x = aMCE.plot_x;
plot_median_y = [aFOE.plot_y_m aDBE.plot_y_m aMCE.plot_y_m];

%% �ṹ��ǽṹ����ռ��

plot_str = sum(plot_median_y(1:8,:))';
plot_nostr = sum(plot_median_y(9:end,:))';
    

%%
%Originģ��ľ���·�����мǱ�д������Ῠ��������
dd = 'E:\OneDrive - tongji.edu.cn\��ʿ��鵵\��������\02-ͼ���ļ�';
Path_Origin = strcat(dd,'\������ʧ-dmg - FOE.opju');
%ͨ��COM�ӿڵ���Origin�������ӻ�����
originObj = actxserver('Origin.ApplicationSI');                                    
invoke(originObj, 'Execute', 'doc -mc 1;');                                           
invoke(originObj, 'IsModified', 'false');                                                
%��Origion ģ��
invoke(originObj, 'Load', Path_Origin);
%������д��Origion��[Book1]Data��
invoke(originObj, 'PutWorksheet', '[All]Sheet1', plot_x);
invoke(originObj, 'PutWorksheet', '[All]Sheet2', plot_y);
invoke(originObj, 'PutWorksheet', '[All]Sheet3', plot_y_m);

%�ͷ�Origion ,����Origion�޷��ر�
release(originObj);

%% ��ʧ�ֲ�-loc
clc
clear

currentDir = pwd; % ��ȡ��ǰĿ¼·��
subfolder = dir('*-*');

for f = 1:length(subfolder)
    
    cd(subfolder(f).name)
    % ���� CSV �ļ�
    filename = 'Output_preprocessing\loss_sample_loc.csv';  % CSV �ļ��� !!!��Ҫ�޸�
    delimiter = ',';  % �ָ��������ţ�
    dataTable = readtable(filename, 'Delimiter', delimiter);
    
    DL = readmatrix('Output_pelicun\DL_summary.csv');
    
    % ��ȡ�ı�����
    % ��ȡ���� "COST" ����
    % columnsWithCost = dataTable(:, startsWith(dataTable.Properties.VariableNames, 'COST'));
    % columnsWithCost = columnsWithCost(:,1:end-1);
    % dataTable.Properties.VariableDescriptions
    % ��ȡ x y ֵ
    x = dataTable{1, 2:end};
    % ��ȡ y ֵ
    y = dataTable{2:end, 2:end};
    
    % % ��ȡ��Ԫ��Ĵ�С
    [numRows, numCols] = size(y);
    % % ����һ���뵥Ԫ���С��ͬ������
    % y_numericData = zeros(numRows, numCols);
    % % ������Ԫ�񣬲������Ԫ��ת��Ϊ����
    % for i = 1:numRows
    %     for j = 1:numCols
    %         y_numericData(i, j) = str2double(y{i, j});
    %     end
    % end
    total_cost = 373000000;
    y_repairable = [];
    for i = 1:numRows
        if DL(i,5)+DL(i,6) == 0
            y_repairable = [y_repairable; y(i,:)];
        else
            y_repairable = [y_repairable; repmat(total_cost/numCols, 1, 61)];
        end
    end
    
    plot_x = x';
    plot_y = y_repairable';
    plot_y_m = median(y_repairable)';
    
    output_folder = 'Output_plot';  % ָ�������ļ���Ŀ¼
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    output_file = fullfile(output_folder, 'loc.mat');  % ��������·��
    save(output_file, 'DL', 'plot_x', 'plot_y', 'plot_y_m');
    
    cd(currentDir)
end

%% ��ͬ¥��ά�޷��öԱ�
clc
clear

currentDir = pwd; % ��ȡ��ǰĿ¼·��
% subfolder = dir('Baseline*');
dd = 'loc.mat';

aFOE = load(strcat(currentDir,'\Baseline building-FOE\Output_plot\',dd));
aDBE = load(strcat(currentDir,'\Baseline building-DBE\Output_plot\',dd));
aMCE = load(strcat(currentDir,'\Baseline building-MCE\Output_plot\',dd));

plot_median_x = aMCE.plot_x;
plot_median_y = [aFOE.plot_y_m aDBE.plot_y_m aMCE.plot_y_m];

%%
%Originģ��ľ���·�����мǱ�д������Ῠ��������
dd = 'E:\OneDrive - tongji.edu.cn\��ʿ��鵵\��������\02-ͼ���ļ�';
Path_Origin = strcat(dd,'\������ʧ-loc - FOE.opju');
%ͨ��COM�ӿڵ���Origin�������ӻ�����
originObj = actxserver('Origin.ApplicationSI');                                    
invoke(originObj, 'Execute', 'doc -mc 1;');                                           
invoke(originObj, 'IsModified', 'false');                                                
%��Origion ģ��
invoke(originObj, 'Load', Path_Origin);
%������д��Origion��[Book1]Data��
invoke(originObj, 'PutWorksheet', '[All]Sheet1', [plot_x plot_y_m]);
invoke(originObj, 'PutWorksheet', '[All]Sheet2', plot_y);
invoke(originObj, 'PutWorksheet', '[All]Sheet3', plot_y_m);

%�ͷ�Origion ,����Origion�޷��ر�
release(originObj);

%% �ָ�����-����
clc
clear

load('File.mat');
k = 5
cd(strcat(File(k).folder,'/',File(k).name,'/Output_treads_R'))
    
    
DT_stepfunc_FR = readmatrix('DT_stepfunc_FR.csv');

plot_x = DT_stepfunc_FR(3:end,2:end);
plot_y = DT_stepfunc_FR(1:2,2:end);

%% �ָ�����-��ͼ
close all
figure

hold on
set(groot, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'Bold', 'DefaultAxesLineWidth', 0.5);

for i = 1:size(plot_x,1)
    if plot_x(i,3) ~= 0
        stairs(plot_x(i,:),plot_y(1,:),'color', [0.5 0.5 0.5])
    else
        stairs(plot_x(i,:),plot_y(2,:),'color', [0.5 0.5 0.5])
    end
end
plot_1 = stairs(plot_x(1,:),plot_y(1,:),'color', [0.5 0.5 0.5])

% �����ݰ���������
sortedData = sort(plot_x(:, end));
% �ҵ��м�������������
numElements = numel(sortedData);
middleIndex1 = numElements / 2;
middleIndex2 = middleIndex1 + 1;
% ȡ�����ϴ���Ǹ�����Ϊ��λֵ�еĽϴ�ֵ
largerValue = max(sortedData(middleIndex1), sortedData(middleIndex2));
% �ҵ���λ�����ڵ��к�
rowIndices = find(plot_x(:, end) == largerValue);

plot_medain = stairs(plot_x(rowIndices,:),plot_y(1,:), 'LineWidth', 1, 'color', [0 0 0]);

% ����ͼ�γߴ�
width = 17;  % ͼ���ȣ�cm��
height = 7;  % ͼ��߶ȣ�cm��
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, width, height]);
% axes('Parent',gcf,'Position',[0.14 0.14 0.54 0.78]);
% ��������������С
set(gca, 'FontSize', 7);


% ���� x �᷶Χ�ͼ��
xlim([0, 1800]);  % ���� x �᷶ΧΪ 3 �� 7
xticks(0:200:1800);  % ���� x ����Ϊ 0.5
ylim([0, 1.0]);  % ���� x �᷶ΧΪ 3 �� 7
yticks(0:0.2:1.0);  % ���� x ����Ϊ 0.5

% ���ͼ��
legend([plot_1 plot_medain],{'Individal','Medain'},'Location','northeast', 'FontName', 'Times New Roman', 'FontSize', 7)


% ���ú������ơ�����ʹ�С
xlabel('Day', 'FontName', 'Times New Roman', 'FontSize', 10);
% �����������ơ�����ʹ�С
ylabel('Building usability', 'FontName', 'Times New Roman', 'FontSize', 10);
% Get handle to the current axes
ax = gca;
ax.Box = 'on';
ax.Layer = 'top';

% You might want to adjust additional properties to ensure visibility
set(ax, 'LineWidth', 0.5); % Makes the axis box lines thicker
% Specify the desired resolution in DPI
resolution = 1200; % You can change this value as needed

% Specify the filename and path
filename = 'treads.png';

% Save the figure with the specified resolution
print(gcf, filename, '-dpng', ['-r' num2str(resolution)]);
%% ��ȡ����
clc
clear

currentDir = pwd; % ��ȡ��ǰĿ¼·��
subfolder = dir('*-*');

f = 1
cd(subfolder(f).name)

DT_stepfunc_FR = readmatrix('Output_treads\DT_stepfunc_FR.csv');

plot_x = DT_stepfunc_FR(3:end,2:end);
plot_y = DT_stepfunc_FR(1:2,2:end);
%% �ָ�����-��ͼ ����
close all
figure

hold on
set(groot, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontWeight', 'Bold', 'DefaultAxesLineWidth', 0.5);

for i = 1:size(plot_x,1)
    if plot_x(i,3) ~= 0
        stairs(plot_x(i,:),plot_y(1,:),'color', [0.5 0.5 0.5])
    else
        stairs(plot_x(i,:),plot_y(2,:),'color', [0.5 0.5 0.5])
    end
end
plot_1 = stairs(plot_x(1,:),plot_y(1,:),'color', [0.5 0.5 0.5])

% �����ݰ���������
sortedData = sort(plot_x(:, end));
% �ҵ��м�������������
numElements = numel(sortedData);
middleIndex1 = numElements / 2;
middleIndex2 = middleIndex1 + 1;
% ȡ�����ϴ���Ǹ�����Ϊ��λֵ�еĽϴ�ֵ
largerValue = max(sortedData(middleIndex1), sortedData(middleIndex2));
% �ҵ���λ�����ڵ��к�
rowIndices = find(plot_x(:, end) == largerValue);

plot_medain = stairs(plot_x(rowIndices,:),plot_y(1,:), 'LineWidth', 1, 'color', [0 0 0]);

% ����ͼ�γߴ�
width = 17;  % ͼ���ȣ�cm��
height = 7;  % ͼ��߶ȣ�cm��

% Get handle to the current axes
ax = gca;
ax.Box = 'on';
ax.Layer = 'top';

set(gcf, 'Units', 'centimeters', 'Position', [5, 5, width, height]);
% axes('Parent',gcf,'Position',[0.14 0.14 0.54 0.78]);
% ��������������С
set(gca, 'FontSize', 7);


% ���� x �᷶Χ�ͼ��
xlim([0, 1800]);  % ���� x �᷶ΧΪ 3 �� 7
xticks(0:200:1800);  % ���� x ����Ϊ 0.5
ylim([0, 1.0]);  % ���� x �᷶ΧΪ 3 �� 7
yticks(0:0.2:1.0);  % ���� x ����Ϊ 0.5

% ���ͼ��
legend([plot_1 plot_medain],{'�����켣','�ж�ʱ����λֵ��Ӧ�Ĺ켣'},'Location','northeast', 'FontName', '����', 'FontSize', 7)


% ���ú������ơ�����ʹ�С
xlabel('��', 'FontName', '����', 'FontSize', 9);
% �����������ơ�����ʹ�С
ylabel('������ʹ����', 'FontName', '����', 'FontSize', 9);

output_folder = 'Output_plot';
output_file = fullfile(output_folder, 'DT_FR��');  % ��������·��

% saveas(gcf, 'DT_FR1-��.png');
print(gcf,output_file,'-dpng','-r1200')
cd(currentDir)

%%
% % ��������
% Tr = 854*2;
% rowIndices = 404
% area_x = [plot_x(:,2:end) repmat(Tr, length(plot_x), 1), repmat(Tr, length(plot_x), 1)];
% area_y = [plot_y(:,2:end) [1;1] [0;0]];
% R_median_p = polyarea(area_x(rowIndices,:),area_y(1,:))/Tr;

%%
saveas(gcf, 'DT_FR1-��.png');
print(gcf,'DT_FR1-��1','-dpng','-r600')
%% �ж�ʱ����ֵ
clc
clear

currentDir = pwd; % ��ȡ��ǰĿ¼·��
subfolder = dir('*-*');

dt_median = [];
for f = 1:length(subfolder)
    
    cd(subfolder(f).name)
    DT = readmatrix('Output_treads\DT_summary.csv');
    
    dt_median = [dt_median; median(DT(3,2))];
    cd(currentDir)
end
%% ��������
clc
clear

currentDir = pwd; % ��ȡ��ǰĿ¼·��
subfolder = dir('*-*');
for f = 1:length(subfolder)
    
    cd(subfolder(f).name)
    
    DT_stepfunc_FR = readmatrix('Output_treads\DT_stepfunc_FR.csv');
    
    plot_x = DT_stepfunc_FR(3:end,2:end);
    plot_y = DT_stepfunc_FR(1:2,2:end);
    
    Tr = 854*2;
    index1 = plot_x(:,3) == 0;
    index2 = plot_x(:,end) > Tr;
    
    area_x = [plot_x(:,2:end) repmat(Tr, length(plot_x), 1), repmat(Tr, length(plot_x), 1)];
    area_y = [plot_y(:,2:end) [1;1] [0;0]];
    
    % i = 3
    % plot(area_x(i,:),area_y(1,:))
    % Area(i,1) = polyarea(area_x(i,:),area_y(1,:));
    for i = 1:size(plot_x,1)
        if plot_x(i,3) ~= 0
            Area(i,:) = polyarea(area_x(i,:),area_y(1,:));
        else
            Area(i,:) = polyarea(area_x(i,:),area_y(2,:));
        end
    end
    
    R = Area./Tr;
    R(index1,1) = 0;
    R(index2,1) = 0;
    R_median(f) = median(R);
    
    cd(currentDir)
end
%%
% R_median_p = polyarea(area_x(rowIndices,:),area_y(1,:))/Tr;
Tr = 854*2;
area_x = [plot_x(:,2:end) repmat(Tr, length(plot_x), 1), repmat(Tr, length(plot_x), 1)];
area_y = [plot_y(:,2:end) [1;1] [0;0]];
R_median_p = polyarea(area_x(rowIndices,:),area_y(1,:))/Tr;


%% ����������������
save('R.mat', 'rowIndices', 'R', 'R_median', 'R_median_p');

%%
plot_all = [x' y'];R_median_p
plot_medain = [plot_y' plot_x(rowIndices,:)'];
%Originģ��ľ���·�����мǱ�д������Ῠ��������
dd = 'E:\OneDrive - tongji.edu.cn\��ʿ��鵵\��������\02-ͼ���ļ�';
Path_Origin = strcat(dd,'\������ʧ-dmg.opju');
%ͨ��COM�ӿڵ���Origin�������ӻ�����
originObj = actxserver('Origin.ApplicationSI');                                    
invoke(originObj, 'Execute', 'doc -mc 1;');                                           
invoke(originObj, 'IsModified', 'false');                                                
%��Origion ģ��
invoke(originObj, 'Load', Path_Origin);
%������д��Origion��[Book1]Data��
invoke(originObj, 'PutWorksheet', '[All]Sheet3', plot_all);
invoke(originObj, 'PutWorksheet', '[Medain]Sheet3', plot_medain);

%�ͷ�Origion ,����Origion�޷��ر�
release(originObj);

%% �޸�ʱ�� - �ֽ�


clc
clear
RT_RSeq_FR = readmatrix('RT_RSeq_FR.csv');

RT_RSeq = RT_RSeq_FR(:,2:end);

RT_RS1 = sum(RT_RSeq(:,1:7:size(RT_RSeq, 2)),2);
RT_RS2 = sum(RT_RSeq(:,2:7:size(RT_RSeq, 2)),2);
RT_RS3 = sum(RT_RSeq(:,3:7:size(RT_RSeq, 2)),2);
RT_RS4 = sum(RT_RSeq(:,4:7:size(RT_RSeq, 2)),2);
RT_RS5 = sum(RT_RSeq(:,5:7:size(RT_RSeq, 2)),2);
RT_RS6 = sum(RT_RSeq(:,6:7:size(RT_RSeq, 2)),2);
RT_RS7 = sum(RT_RSeq(:,7:7:size(RT_RSeq, 2)),2);

RT_RS = [RT_RS1 RT_RS2 RT_RS3 RT_RS4 RT_RS5 RT_RS6 RT_RS7];
RT_RS_m = median(RT_RS)';

% Ӱ�����طֽ�
N = length(RT_RS1);
Impeding = readmatrix('IF_delays.csv');
Impeding_I = Impeding(1:N,2);
Impeding_E = Impeding(1:N,3);
Impeding_P = Impeding(1:N,4);
Impeding_F = Impeding(1:N,5);
Impeding_C1 = Impeding(1:N,6);
Impeding_C2 = Impeding(1:N,7);
Impeding_C3 = Impeding(1:N,8);
Impeding_C4 = Impeding(1:N,9);
Impeding_C5 = Impeding(1:N,10);
Impeding_C6 = Impeding(1:N,11);
Impeding_C7 = Impeding(1:N,12);

IP = [Impeding_I Impeding_E Impeding_P Impeding_F Impeding_C1 Impeding_C2 Impeding_C3 Impeding_C4 Impeding_C5 Impeding_C6 Impeding_C7];
IP_m = median(IP)';

Result = [IP_m; RT_RS_m];

%% �ָ�ʱ��-�ֽ� ����
clc
clear

currentDir = pwd; % ��ȡ��ǰĿ¼·��
subfolder = dir('Base*');

Result = [];
for f = 1:length(subfolder)
    
    cd(subfolder(f).name)
    % load('E:\OneDrive - tongji.edu.cn\��ʿ��鵵\�������\SAUSAGE 61-Story - R\DT_m.mat');
    DT_s = readmatrix('Output_treads\DT_stepfunc_FR.csv');
    nonZeroValues = DT_s(DT_s(:, 4) ~= 0, 4);
    N = length(nonZeroValues);
    DT = DT_s(3:N+2,end);
    
    Impeding = readmatrix('Output_treads\IF_delays.csv');
    Impeding_I = Impeding(1:N,2);
    Impeding_EP = Impeding(1:N,3)+Impeding(1:N,4);
    Impeding_F = Impeding(1:N,5);
    Impeding_C = max(Impeding(1:N,6:12),[],2);
    
    Utility_s = readmatrix('Output_treads\DT_path_FR.xlsx','Sheet','utility');
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
    
    Result = [Result; Utility_m Impeding_m Repair_m]
    % median_value = DT_m(i);
    % % �ҵ���λֵ���ڵ�λ��
    % median_position = find(DT == median_value);
    % Result1 = [Utility(median_position) Imp(median_position) Repair(median_position)]
    cd(currentDir)
end