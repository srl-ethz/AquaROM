% ------------------------------------------------------------------------ 
% B_shape_optimisation.m
%
% Description: 3D shape optimisation of a fish, to obtain best swimming
% performance.
% 
% Last modified: 04/02/2026, Mathieu Dubied, ETH Zurich
% ------------------------------------------------------------------------
clear; 
close all; 
clc
if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

%% PREPARE MODELS _________________________________________________________                                                   

% load material parameters
load('parameters.mat') 

% specify and create FE mesh
n_elements = 8086;
filename = strcat('InputFiles/3d_rectangle_', num2str(n_elements), 'el');

[MeshNominal, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
[Lx, Ly, Lz] = mesh_dimensions_3D(nodes);
n_elements = size(elements,1);

if n_elements == 8086
    expName1 = "SO1";
    expName2 = "SO2";
    expName3 = "SO3"; 
elseif n_lemenets == 16009
    expName1 = "SO4";
    expName2 = "SO5";
    expName3 = "SO6"; 
else
    fprintf("n_elements must be either 8086 or 16009 to replicate the results.")
end

%% SHAPE VARIATIONS _______________________________________________________

[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish, xz_concaveTail, y_fin...
    y_bumpBack, y_bumpFront] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);

% plot an example, with SO1 parameters
U = [z_tail,z_head,y_thinFish];
xiPlot = [-0.5;0.5,;-0.5];
U = [y_fin];
xiPlot = [0.5];
f1 = figure('units','centimeters','position',[3 3 10 7],'name','Shape-varied mesh');
elementPlot = elements(:,1:4); hold on 
v1 = reshape(U*xiPlot, 3, []).';
S = 1.0;
hf = PlotFieldonDeformedMesh_ext(nodes, elementPlot, v1, 'factor', 1.0,'lineWidth',0.2);
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];
plotcube(L,O,.05,[0 0 0]);
axis equal; grid on; box on; 
set(f1,'PaperUnits','centimeters');
set(f1,'Units','centimeters');

%% OPTIMIZATION PARAMETERS
h = 0.02;
tmax = 2.0;
kActu = 3.1*1e5;    % multiplicative factor for the actuation forces (8086 elements)
kActu = 7.4*1e4;
% kActu = 2.8*1e5;       % value for 16009 elements

%% OPTIMISATION SO1 _______________________________________________________

% Preparation (shape variation and opt. constraints)
U_1 = [z_tail,z_head,y_thinFish]; 
A = [1 0 0 ;
    -1 0 0;
    0 1 0;
    0 -1 0;
    0 0 1;
    0 0 -1];
b = [0.5;0.5;0.5;0.5;0.15;0.15];
barrierParam = 10*ones(1,length(b));

% Optimisation
tStart = tic;
out1 = optimise_shape_3D(myElementConstructor,nsetBC, ...
    nodes,elements,muscleBoundaries,kActu,U_1,h,tmax,A,b, ...
    'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, ...
    'maxIteration',40, ...
    'convCrit',0.004, ... % 0.004 for 8086, 0.01 for 16009
    'convCritCost',1.0, ... % 1.0 for 8086, 2.0 for 16009
    'barrierParam',barrierParam, ...
    'gStepSize',0.0004, ...  % 0.0004 for 8086, 0.0002 for 16609
    'nRebuild',6, ...
    'rebuildThreshold',0.15,...
    'wSize', 5, ...
    'USEJULIA',1, ...
    'paramInit', [0.4;0.2;0.1]);
out1.tOpti = toc(tStart)/60;   % unit is minute
out1.tPerIt = out1.tOpti/out1.nIt;

% Print stats
fprintf('Computation time: %.2fmin\n',out1.tOpti)
fprintf('Number of built models and solved EoMs: %3d\n',out1.nIt)
fprintf('Computation time per models/EoMs: %.2f\n',out1.tPerIt)

% Save results
filename = sprintf('Results/Data/B_shape_optimization/%s_%d_el_kActu_%.3f.mat', ...
                    expName1, n_elements, kActu);
save(filename,'out1')

%% OPTIMISATION SO2 _______________________________________________________

% Preparation (shape variation and opt. constraints)
U_2 = [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish];
% U_2 = [z_tail,z_head,y_thinFish,y_head,y_ellipseFish];
nParam = 5;
A = zeros(2 * nParam, nParam);
for i = 1:nParam
    A(2*i-1:2*i,i) =[1;-1];
end
% yTotConstr = [0 0 1 0 1;0 0  -1 0 -1;
%               0 0 0 1 1;0 0  0 -1 -1];
          
% A = [A;yTotConstr];
% b = [0.5;0.5;
%     0.5;0.5;
%     0.4;0.4;
%     0.5;0.5;
%     0.5;0.5;   
%     0.8;0.8;
%     0.8;0.8];
b = [0.5;0.5;
    0.5;0.5;
    0.4;0.4;
    0.5;0.5;
    0.5;0.5];
% b = [0.5;0.5;
%     0.5;0.5;
%     0.15;0.15;
%     0.3;0.3;
%     0.3;0.3];
barrierParam = 2*ones(1,length(b));

% Optimisation
tStart = tic;
out2 = optimise_shape_3D(myElementConstructor,nsetBC, ...
    nodes,elements,muscleBoundaries,kActu,U_2,h,tmax,A,b,...
    'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, ...
    'maxIteration',50, ... 
    'convCrit',0.01, ...   % set to 0.01 for 8086el, 0.01 for 16609
    'convCritCost',0.5, ... % set to 0.5 for 4270 el, 0.5 for 16609
    'barrierParam',barrierParam, ...
    'gStepSize',0.0004,...   % set to 0.0005 for 8086, 0.0002 for 16609
    'nRebuild',5, ...
    'rebuildThreshold',0.15,...
    'wSize', 5, ...
    'USEJULIA',1);
out2.tOpti = toc(tStart)/60;   % unit is minute
out2.tPerIt = out2.tOpti/out2.nIt;

% Print stats
fprintf('Computation time: %.2fmin\n',out2.tOpti)
fprintf('Number of built models and solved EoMs: %3d\n',out2.nIt)
fprintf('Computation time per models/EoMs: %.2f\n',out2.tPerIt)

% Save results
filename = sprintf('Results/Data/B_shape_optimization/%s_%d_el_kActu_%.3f.mat', ...
                    expName2, n_elements, kActu);
save(filename,'out2')

%% OPTIMISATION SO3 _______________________________________________________

% Preparation (shape variation and opt. constraints)
U_3 = [z_tail,z_head,...
    y_linLongTail,y_head,y_ellipseFish,...
    z_smallFish, z_notch, xz_concaveTail];
% U_3 = [z_tail,z_head,y_thinFish,y_head,y_ellipseFish,y_linLongTail,z_notch,xz_concaveTail];


nParam = 8;
A = zeros(2 * nParam, nParam);
for i = 1:nParam
    A(2*i-1:2*i,i) =[1;-1];
end
b = [0.5;0.5;
    0.5;0.5;
    0.3;0.3;
    0.4;0.4;
    0.3;0.3;
    0.3;0.3 ;
    0.3;0.3;
    0.2;0.01];
% b = [0.5;0.5;
%     0.5;0.5;
%     0.35;0.35;  
%     0.3;0.3;
%     0.4;0.4;    
%     0.35;0.35;
%     0.3;0.3;
%     0.2; 0.01];  % concave tail only in one direction 

barrierParam = 2*ones(1,length(b)); % 3 for 8086, 1  for 16009

% Optimisation
tStart = tic;
out3 = optimise_shape_3D(myElementConstructor,nsetBC, ...
    nodes,elements,muscleBoundaries,kActu,U_3,h,tmax,A,b, ...
    'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, ...
    'maxIteration',40, ...
    'convCrit',0.015, ...    % set to 0.015 for 8086el, 0.015 for 16609
    'convCritCost',5.0, ...% set to 5.0 for 8086, 8.0 for 16609
    'barrierParam',barrierParam, ...
    'gStepSize',0.0003, ... % set to 0.0003 for 8086, 0.0003 for 16609
    'nRebuild',5, ...
    'rebuildThreshold',0.15, ...
    'wSize', 5, ...
    'USEJULIA',1);
out3.tOpti = toc(tStart)/60;   % unit is minute
out3.tPerIt = out3.tOpti/out3.nIt;

% Print stats
fprintf('Computation time: %.2fmin\n',out3.tOpti)
fprintf('Number of built models and solved EoMs: %3d\n',out3.nIt)
fprintf('Computation time per models/EoMs: %.2f\n',out3.tPerIt)

% Save results
filename = sprintf('Results/Data/B_shape_optimization/%s_%d_el_kActu_%.3f.mat', ...
                    expName3, n_elements, kActu);
save(filename,'out3')

%% RETRIEVE RESULTS _______________________________________________________
filename = sprintf('Results/Data/B_shape_optimization/%s_%d_el_kActu_%.3f.mat', ...
                    expName1, n_elements, kActu);
load(filename)
filename = sprintf('Results/Data/B_shape_optimization/%s_%d_el_kActu_%.3f.mat', ...
                    expName2, n_elements, kActu);
load(filename)
filename = sprintf('Results/Data/B_shape_optimization/%s_%d_el_kActu_%.3f.mat', ...
                    expName3, n_elements, kActu);
load(filename)
%% PLOT OPTIMAL SHAPE _____________________________________________________

SOIdx =   1;
switch SOIdx
    case 1
        xiStar = out1.xiStar;
        U = U_1; 
    case 2
        xiStar = out2.xiStar;
        U = U_2;
    case 3
        xiStar = out3.xiStar;
        U = U_3;
end

f_opt_shape = figure('units','centimeters','position',[3 3 9 4]);
elementPlot = elements(:,1:4); hold on
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];

t = tiledlayout(1, 2, 'TileSpacing', 'None', 'Padding', 'tight', 'InnerPosition',[0.01, 0.01 ,0.98,0.98]);
v = reshape(U*xiStar, 3, []).';

% View 1
ax1 = nexttile(t);
PlotFieldonDeformedMesh_ext(nodes, elementPlot, v, 'factor', 1, 'lineWidth',0.2);
plotcube(L,O,.05,[0 0 0]);

% View 2
ax2 = nexttile(t);
PlotFieldonDeformedMesh_ext(nodes, elementPlot, v, 'factor', 1, 'lineWidth',0.2);
plotcube(L,O,.05,[0 0 0]);
view(ax2,[23.28863587057058,35.063216347984103])

axis([ax1 ax2],[-0.23 0 -0.02 0.02 -0.1 0.1])

fig_title = sprintf('Results/Figures/SO%d_opt_shape_%d_el_kActu_%.3f', SOIdx, n_elements, kActu);
% exportgraphics(f_opt_shape,strcat(fig_title,'.pdf'),'Resolution',1200)

%% PRINT RESULTS AND PERFORMANCE STATS ____________________________________
SOIdx = 3;
switch SOIdx
    case 1
        xiStar = out1.xiStar;
        tOpti = out1.tOpti;
        nIt = out1.nIt;
        tPerIt = out1.tPerIt;    
    case 2
        xiStar = out2.xiStar;
        tOpti = out2.tOpti;
        nIt = out2.nIt;
        tPerIt = out2.tPerIt;   
    case 3
        xiStar = out3.xiStar;
        tOpti = out3.tOpti;
        nIt = out3.nIt;
        tPerIt = out3.tPerIt; 
end
fprintf('Computation time: %.2f min\n', tOpti);
fprintf('Number of built models and solved EoMs: %d\n', nIt);
fprintf('Computation time per built models/EoMs: %.2f min\n', tPerIt)
disp('Optimal parameters:')
disp(xiStar)
%% PLOT COST FUNCTION FOR THE 3 EXPERIMENTS _______________________________
f_cost = figure('units','centimeters','position',[3 3 9 6]);
hold on
plot(out1.LwoBEvo,'LineWidth',1)
plot(out2.LwoBEvo,'--','LineWidth',1)
plot(out3.LwoBEvo,'-.','LineWidth',1)
grid on
ylabel('$$L$$')
xlabel('Iterations')
legend('SO1','SO2','SO3')
hold off

fig_title = sprintf('Results/Figures/B_shape_optimization/SO_cost_evolution_%d_el.pdf', n_elements);
exportgraphics(f_cost,fig_title,'Resolution',1200)

%% PLOT SWIMMING DISTANCE FOR THE 3 EXPERIMENTS ___________________________
f_dist = figure('units','centimeters','position',[3 3 9 6]);
hold on
p1 = plot(out1.xEvo*100,'LineWidth',1);
p2 = plot(out2.xEvo*100,'--','LineWidth',1);
p3 = plot(out3.xEvo*100,'-.','LineWidth',1);
grid on
ylabel('Swimming distance [cm]')
xlabel('Iterations')
legend([p1,p2,p3],'SO1','SO2','SO3', 'Location', 'southeast')
hold off

fig_title = sprintf('Results/Figures/SO_dist_evolution_%d_el.pdf', n_elements);
% exportgraphics(f_dist,fig_title,'Resolution',1200)


%% PLOT NORMALIZED SWIMMING DISTANCE FOR THE 3 EXPERIMENTS ________________
% normalized by the quantity of energy available, i.e., the size of the
% muscles
f_norm_dist = figure('units','centimeters','position',[3 3 9 6]);
hold on
p1 = plot(out1.xEnergyEvo*100,'LineWidth',1);
p2 = plot(out1.xEnergyEvo*100,'--','LineWidth',1);
p3 = plot(out1.xEnergyEvo*100,'-.','LineWidth',1);
grid on
ylabel('Normalized swimming distance')
xlabel('Iterations')
legend([p1,p2,p3], 'SO1','SO2','SO3', 'Location', 'northwest')
hold off

fig_title = sprintf('Results/Figures/SO_normalized_dist_evolution_%d_el.pdf', n_elements);
% exportgraphics(f_norm_dist,fig_title,'Resolution',1200)
