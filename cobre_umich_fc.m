%%
%%  COBRE and Umich Psychosis Data
%%
%%

clear all; clc; close all;

load mycmap.mat
load COBRE_FC_FIN.mat;

%% communities
Power14indsr=[13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41 254 ...
    42 43 44 45 46 ...
    47 48 49 50 51 52 53 54 55 56 57 58 59 60 ...
    61 62 63 64 65 66 67 68 69 70 71 72 73 ...
    74  75  76  77  78  79  80  81  82  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 136 138 ...
    132 133 134 135 220 ...
    142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 ...
    173 174 175 176 177 178 179 180 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 ...
    202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 ...
    221 222 223 224 225 226 227 228 229 230 231 232 233 ...
    137 234 235 236 237 238 239 240 241 ...
    250 251 255 256 257 258 259 260 261 262 263 ...
    242 243 244 245 ...
    1   2   3   4   5   6   7   8   9  10  11  12  83  84 131 139 140 141 181 182 183 184 246 247 248 249 252 253];
Power14lbls=[ones(30,1);ones(5,1)*2;ones(14,1)*3;ones(13,1)*4;ones(57,1)*5;ones(5,1)*6;ones(31,1)*7;ones(25,1)*8;ones(18,1)*9;ones(13,1)*10;ones(9,1)*11;ones(11,1)*12;ones(4,1)*13;ones(28,1)*14];
Power14indsr=Power14indsr';

[xxy,yyy] = grid_communities(Power14lbls);
N = 263;
lcol = [0.2 0.2 0.2];

test=zeros(263, 263);
figure, 
imagesc(test(Power14indsr,Power14indsr),[-0.6 0.8])
hold on
plot(xxy,yyy,'color',lcol,'linewidth',2);
set(gca,'Colormap',flipud(mycmap),'XTick',[],'YTick',[])
%xlabel('Nodes','FontSize',14)
axis square
cb = colorbar; set(cb,'Ticks',[-0.6 0 0.8],'FontSize',12)
set(gca, 'fontweight', 'bold')
set(gcf, 'color', 'w')


figure, 
subplot(131)
imagesc(HC(Power14indsr,Power14indsr),[-0.6 0.8])
hold on
plot(xxy,yyy,'color',lcol,'linewidth',2);
set(gca,'Colormap',flipud(mycmap),'XTick',[],'YTick',[])
xlabel('Nodes','FontSize',14)
axis square
cb = colorbar; set(cb,'Ticks',[-0.6 0 0.8],'FontSize',12)
title('HC','FontSize',16)
set(gca, 'fontweight', 'bold')
set(gcf, 'color', 'w')
subplot(132)
imagesc(SZ(Power14indsr,Power14indsr),[-0.6 0.8])
hold on
plot(xxy,yyy,'color',lcol,'linewidth',2);
set(gca,'Colormap',flipud(mycmap),'XTick',[],'YTick',[])
xlabel('Nodes','FontSize',14)
axis square
cb = colorbar; set(cb,'Ticks',[-0.6 0 0.8],'FontSize',12)
title('SZ','FontSize',16)
set(gca, 'fontweight', 'bold')
set(gcf, 'color', 'w')
dd=HC-SZ;
subplot(133)
imagesc(dd(Power14indsr,Power14indsr),[-0.2 0.2])
hold on
plot(xxy,yyy,'color',lcol,'linewidth',2);
set(gca,'Colormap',flipud(mycmap),'XTick',[],'YTick',[])
xlabel('Nodes','FontSize',14)
axis square
cb = colorbar; set(cb,'Ticks',[-0.2 0 0.2],'FontSize',12)
title('HC-SZ','FontSize',16)
set(gca, 'fontweight', 'bold')
set(gcf, 'color', 'w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot as tree
% generate minimum spanning tree
mst = graphminspantree(sparse(max(HC(:)) - HC),'method','kruskal');
% node properties
rngSz = [5,50];                     % range of node sizes
nodeSz = fcn_sz(sum(HC,2),rngSz);    % node sizes
% visualization parameters
mapSz = 1000;    % size of map
sigma = 10;      % radius of smoothed colors
% generate thresholded network and union with mst
thr = 0.02;
HC = double(threshold_proportional(HC,thr) | mst | mst');
% make figure
load layoutData.mat
[f,ax,nodeHandles,edgeHandles] = fcn_community_boundaries(HC,Power14lbls,cmap,nodeSz,mapSz,sigma);
set(gcf, 'color', 'w')


% generate minimum spanning tree
mst = graphminspantree(sparse(max(SZ(:)) - SZ),'method','kruskal');
nodeSz = fcn_sz(sum(SZ,2),rngSz);    % node sizes
SZ = double(threshold_proportional(SZ,thr) | mst | mst');
[f,ax,nodeHandles,edgeHandles] = fcn_community_boundaries(SZ,Power14lbls,cmap,nodeSz,mapSz,sigma);
set(gcf, 'color', 'w')


% generate minimum spanning tree
mst = graphminspantree(sparse(max(dd(:)) - dd),'method','kruskal');
nodeSz = fcn_sz(sum(dd,2),rngSz);    % node sizes
dd = double(threshold_proportional(dd,thr) | mst | mst');
[f,ax,nodeHandles,edgeHandles] = fcn_community_boundaries(dd,Power14lbls,cmap,nodeSz,mapSz,sigma);
set(gcf, 'color', 'w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defining groups of nodes
A_group = Power14lbls;
A_group_cols = [
    0.9290, 0.6940, 0.1250; % Yellow
    0, 0.4470, 0.7410;      % Blue
    0.8500, 0.3250, 0.0980; % Orange
    0.4940, 0.1840, 0.5560; % Purple
    0.4660, 0.6740, 0.1880; % Green
    0.3010, 0.7450, 0.9330; % Light Blue
    0.6350, 0.0780, 0.1840; % Dark Red
    0, 0.75, 0.75;          % Teal
    1, 0, 0;                % Red
    0, 1, 0;                % Bright Green
    0, 0, 1;                % Bright Blue
    0.75, 0.75, 0;          % Olive
    0.75, 0, 0.75;          % Dark Magenta
    0, 0.75, 0              % Dark Green
];
figure, circos_multilayer(HC(Power14indsr,Power14indsr), 'group', A_group, 'group_color', A_group_cols);
figure, circos_multilayer(SZ(Power14indsr,Power14indsr), 'group', A_group, 'group_color', A_group_cols);
figure, circos_multilayer(dd(Power14indsr,Power14indsr), 'group', A_group, 'group_color', A_group_cols);

%%You need cocoanCORE-master
load('Power_L.mat')
view_parcellation(parcels, 'L', 'inflated');set(gcf, 'color', 'w')
load('Power_R.mat')
view_parcellation(parcels, 'R', 'inflated');set(gcf, 'color', 'w')

load('Power_L.mat')
view_parcellation(parcels, 'L', 'sphere');set(gcf, 'color', 'w')
load('Power_R.mat')
view_parcellation(parcels, 'R', 'sphere');set(gcf, 'color', 'w')

figure;
nuPoints = 100;
sizes = linspace(10, 100, numPoints);
theta = linspace(0, 2*pi, numPoints);
cmap = jet(numPoints);
scatter(theta, ones(1, numPoints), sizes, 1:numPoints, 'filled', 'MarkerEdgeColor', 'k');
colormap(cmap);
c = colorbar;
caxis([1, numPoints]);
c.Ticks = linspace(1, numPoints, 6);
c.TickLabels = arrayfun(@(x) sprintf('%.2f', x/numPoints), c.Ticks, 'UniformOutput', false);
ax = gca;
ax.Visible = 'off';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load COBRE_Seg.mat
data1=[SHS, SHW, SHB];
data2=[SSS, SSW, SSB];

figure,
gretna_plot_violin({SHS,SSS},{'HC','SZ'},{' '},'dotfill')
ylabel('Brain System Segregation Score')
ylim([0.2 0.6])
set(gca, 'fontweight', 'bold')
set(gcf, 'color', 'w')

gretna_plot_violin({data1,data2},{'HC','SZ'},{'system segreation', 'mean correlation with community', 'mean correlation between community'},'boxfill')
ylabel('brain system segregation score')
ylim([0 0.6])
set(gca, 'fontweight', 'bold')
set(gcf, 'color', 'w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load COBRE_Seg.mat
data1=[SHS, SHW, SHB];
data2=[SSS, SSW, SSB];

% Your supply IDs
supplyIDs = [ 40000 40002 40003 40004 40006 40008 40010 40012 40013 40014 40015 40017 40018 40019 40020 40021 40022 40023 40025 40026 40027 40028 40029 40030 40031 ...
    40032 40033 40034 40035 40037 40038 40041 40042 40043 40044 40045 40046 40047 40048 40049 40050 40051 40052 40053 40054 40055 ...
    40056 40057 40058 40059 40061 40062 40063 40065 40066 40067 40068 40069 40072 40074 40076 40077 40078 40079 40080 40081 40082 40085 40086 ...
    40087 40090 40091 40092 40093 40094 40095 40096 40098 40099 40100 40101 40102 40103 40104 40106 40107 40108 40109 40110 40112 40113 40114 ...
    40115 40116 40117 40118 40119 40120 40121 40122 40123 40124 40125 40126 40127 40128 40129 40130 40131 40132 40133 40134 40135 40137 40138 ...
    40139 40140 40141 40142 40143 40144 40145 40146 40147];

% Your labels
labels = [1  1  1  1  1  1  1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1  1 -1 -1  1 -1  1 -1  1 -1  1  1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1 ...
         -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1  1  1  1  1  1 -1 -1 -1 -1  1 -1  1 -1  1  1  1  1  1 -1  1 -1  1 -1  1  1  1  1 -1 -1 ...
         -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1 -1 -1];

dataTable = table(supplyIDs', labels', 'VariableNames', {'SupplyID', 'Label'});

filteredTable = dataTable(dataTable.Label == -1, :);
SegCore = data1(:,1);
filteredTable.SegCore = SegCore;
MCI = data1(:,2);
filteredTable.MCI = MCI;
MCB = data1(:,3);
filteredTable.MCB = MCB;

Age = [34 31 30 ...
     47 ...
     44 ...
     22 ...
     48 ...
     44 ...
     48 ...
     43 ...
     43 ...
     31 ...
     30 ...
     53 ...
     38 ...
     30 ...
     31 ...
     36 ...
     23 ...
     22 ...
     24 ...
     52 ...
     30 ...
     27 ...
     36 ...
     27 ...
     18 ...
     50 ...
     37 ...
     22 ...
     62 ...
     33 ...
     24 ...
     58 ...
     34 ...
     52 ...
     65 ...
     27 ...
     18 ...
     24 ...
     25 ...
     40 ...
     40 ...
     26 ...
     33 ...
     20 ...
     23 ...
     27 ...
     47 ...
     34 ...
     44 ...
     26 ...
     21 ...
     48 ...
     35 ...
     48 ...
     26 ...
     22 ...
     23 ...
     47 ...
     47 ...
     42 ...
     24 ...
     28 ...
     28 ...
     53 ...
     35 ...
     54 ...
     39 ...
     34];
filteredTable.Age = Age';
%disp(filteredTable);
HC_Seg=filteredTable;

filteredTables = dataTable(dataTable.Label == 1, :);

SegCores = data2(:,1);
filteredTables.SegCores = SegCores;
MCIs = data2(:,2);
filteredTables.MCIs = MCIs;
MCBs = data2(:,3);
filteredTables.MCBs = MCBs;
Ages = [20 ...
     19 ... 
     28 ...
     55 ...
     53 ...
     28 ...
     52 ...
     32 ...
     47 ...
     21 ...
     23 ...
     33 ...
     64 ...
     20 ...
     31 ...
     29 ...
     24 ...
     62 ...
     40 ...
     48 ...
     18 ...
     37 ...
     44 ...
     25 ...
     25 ...
     28 ...
     57 ...
     26 ...
     26 ...
     43 ...
     50 ...
     22 ...
     49 ...
     57 ...
     22 ...
     35 ...
     38 ...
     35 ...
     50 ...
     40 ...
     46 ...
     29 ...
     33 ...
     43 ...
     42 ...
     19 ...
     22 ...
     41 ...
     50 ...
     20 ...
     21 ...
     23 ...
     52 ...
     19];
filteredTables.Ages = Ages';
%disp(filteredTable);
SZ_Seg=filteredTables;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,
gretna_plot_regression_spec(table2array(SZ_Seg(:, 6)), table2array(SZ_Seg(:, 3)), 1, 'off');
hold on,
gretna_plot_regression_spec1(table2array(SZ_Seg(:, 6)), table2array(SZ_Seg(:, 4)), 1, 'off');
hold on,
gretna_plot_regression_spec2(table2array(SZ_Seg(:, 6)), table2array(SZ_Seg(:, 5)), 1, 'off');
xlim([18, 64]);
xticks(18:3:64);
ylim([0, 0.6])
yticks(0:0.1:0.6);
xlabel('Age in Years');
ylabel('brain system segregation score')
set(gca, 'fontweight', 'bold');
set(gcf, 'color', 'w');

figure,
gretna_plot_regression_spec(table2array(HC_Seg(:, 6)), table2array(HC_Seg(:, 3)), 1, 'off');
hold on,
gretna_plot_regression_spec1(table2array(HC_Seg(:, 6)), table2array(HC_Seg(:, 4)), 1, 'off');
hold on,
gretna_plot_regression_spec2(table2array(HC_Seg(:, 6)), table2array(HC_Seg(:, 5)), 1, 'off');
xlim([18, 65]);
xticks(18:3:65);
ylim([0, 0.6])
yticks(0:0.1:0.6);
xlabel('Age in Years');
ylabel('brain system segregation score')
set(gca, 'fontweight', 'bold');
set(gcf, 'color', 'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PNAS, Wang et al,. https://www.pnas.org/doi/10.1073/pnas.2022288118
[Clus_numHC,Clus_sizeHC,FCHC] = Functional_HP(HC,263);
[Clus_numSZ,Clus_sizeSZ,FCSZ] = Functional_HP(SZ,263);

[HinHC,HseHC,pHC] =Balance(HC,263,Clus_sizeHC,Clus_numHC);
[HinSZ,HseSZ,pSZ] =Balance(SZ,263,Clus_sizeSZ,Clus_numSZ);

load COBRE_FC_FIN_Subjs.mat
for i=1:size(HC,3)
    i
    [Clus_numHC,Clus_sizeHC,FCHC] = Functional_HP(HC(:,:,i),263);
    [HinHCs(:,i),HseHCs(:,i),~] =Balance(HC(:,:,i),263,Clus_sizeHC,Clus_numHC);
end

for i=1:size(SZ,3)
    i
    [Clus_numSZ,Clus_sizeSZ,FCSZ] = Functional_HP(SZ(:,:,i),263);
    [HinSZs(:,i),HseSZs(:,i),~] =Balance(SZ(:,:,i),263,Clus_sizeSZ,Clus_numSZ);
end

figure,
gretna_plot_regression_spec(table2array(HC_Seg(:, 6)), HinSZ', 1, 'off');  % HERE PUT HinSZ IS RIGHT, JUST MY SAVE NAME PROBLEM!!!
hold on;
gretna_plot_regression_spec1(table2array(SZ_Seg(:, 6)), HinHC', 1, 'off');  % HERE PUT HinSZ IS RIGHT, JUST MY SAVE NAME PROBLEM!!!
xlim([18, 65]);
xticks(18:3:65);
ylim([0, 0.03])
yticks(0:0.01:0.03);
xlabel('Age in Years');
ylabel('brain system integration score')
set(gca, 'fontweight', 'bold');
set(gcf, 'color', 'w');

figure,
gretna_plot_violin({HinHCs',HinSZs'},{'SZ','HC'},{'system integration'},'boxfill')
ylabel('brain system integration score')
ylim([0 0.03])
set(gca, 'fontweight', 'bold')
set(gcf, 'color', 'w')

figure,
gretna_plot_violin({HseHCs',HseSZs'},{'SZ','HC'},{'system segregation'},'boxfill')
ylabel('brain system segregation score')
ylim([0 0.03])
set(gca, 'fontweight', 'bold')
set(gcf, 'color', 'w')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Calibration
load COBRE_FC_FIN.mat;
[Hin_c,Hse_c] = Stable_correct(HC,HinHCs,HseHCs,263);
[sHin_c,sHse_c] = Stable_correct(SZ,HinSZs,HseSZs,263);

figure,
gretna_plot_violin({sHin_c',Hin_c'},{'HC', 'SZ'},{'system integration'},'boxfill')
ylabel('brain system integration score')
ylim([0 0.01])
yticks([0 0.01])
set(gca, 'fontweight', 'bold')
set(gcf, 'color', 'w')

figure,
gretna_plot_violin({sHse_c',Hse_c'},{'HC', 'SZ'},{'system segregation'},'boxfill')
ylabel('brain system segregation score')
ylim([0.1 0.4])
yticks([0.1 0.4])
set(gca, 'fontweight', 'bold')
set(gcf, 'color', 'w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Calibration
h_b=abs(table2array(HC_Seg(:, 3)) - HinSZ');
s_b=abs(table2array(SZ_Seg(:, 3)) - HinHC');

gretna_plot_violin({h_b, s_b},{'HC','SZ'},{'system balance'},'boxfill')
ylabel('brain system integration score')
ylim([0 0.6])
set(gca, 'fontweight', 'bold')
set(gcf, 'color', 'w')

figure,
gretna_plot_regression_spec(table2array(HC_Seg(:, 6)), h_b, 1, 'off');  % HERE PUT HinSZ IS RIGHT, JUST MY SAVE NAME PROBLEM!!!
hold on;
gretna_plot_regression_spec1(table2array(SZ_Seg(:, 6)), s_b, 1, 'off');  % HERE PUT HinSZ IS RIGHT, JUST MY SAVE NAME PROBLEM!!!
xlim([18, 65]);
xticks(18:3:65);
ylim([0, 0.6])
yticks(0:0.3:0.6);
xlabel('Age in Years');
ylabel('brain system balance score')
set(gca, 'fontweight', 'bold');
set(gcf, 'color', 'w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
