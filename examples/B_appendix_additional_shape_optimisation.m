% ------------------------------------------------------------------------ 
% B_appendix_additional_shape_optimisation.m
%
% Description: 3D shape optimisation of a fish, to obtain best swimming
% performance. Analysis of the shape variations presented in the appendix.
% 
% Last modified: 08/02/2026, Mathieu Dubied, ETH Zurich
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

%% SHAPE VARIATIONS _______________________________________________________

[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish, xz_concaveTail, y_fin, ...
    y_bumpBack, y_bumpFront] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);

[z_tail_1, z_tail_2] = ...
    shape_variations_3D_extra(nodes,Lx,Ly,Lz);


% plot an example, with SO1 parameters
U = [z_tail,z_head,y_thinFish];
xiPlot = [-0.5;0.5,;-0.5];
U = [y_bumpBack];
xiPlot = [-0.2];
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
kActu = 6.0*1e4;    % multiplicative factor for the actuation forces (8086 elements)


%% OPTIMISATION SO7 _______________________________________________________

% Preparation (shape variation and opt. constraints)
U_7 = [z_tail,z_head,y_thinFish, y_fin, y_bumpBack]; 
nParam = 5;
A = zeros(2 * nParam, nParam);
for i = 1:nParam
    A(2*i-1:2*i,i) =[1;-1];
end
b = [0.5;0.5;
    0.5;0.5;
    0.15;0.15;
    0.3;0.02;
    0.3;0.3];
barrierParam = 5*ones(1,length(b));


% Optimisation
tStart = tic;
out7 = optimise_shape_3D(myElementConstructor,nsetBC, ...
    nodes,elements,muscleBoundaries,kActu,U_7,h,tmax,A,b, ...
    'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, ...
    'maxIteration',50, ...
    'convCrit',0.004, ... 
    'convCritCost',1.0, ...
    'barrierParam',barrierParam, ...
    'gStepSize',0.0003, ... 
    'nRebuild',8, ...
    'rebuildThreshold',0.15,...
    'wSize', 5, ...
    'USEJULIA',1);
out7.tOpti = toc(tStart)/60;   % unit is minute
out7.tPerIt = out7.tOpti/out7.nIt;

% Print stats
fprintf('Computation time: %.2fmin\n',out7.tOpti)
fprintf('Number of built models and solved EoMs: %3d\n',out7.nIt)
fprintf('Computation time per models/EoMs: %.2f\n',out7.tPerIt)

% Save results
filename = sprintf('Results/Data/Appendix/Additional_shape_optimization/SO7_%d_el_kActu_%.3f.mat', ...
                   n_elements, kActu);
save(filename,'out7')


%% OPTIMISATION SO8 _______________________________________________________

% Preparation (shape variation and opt. constraints)
U_8 = [z_tail,z_head,y_thinFish, y_linLongTail, y_head]; 
nParam = 5;
A = zeros(2 * nParam, nParam);
for i = 1:nParam
    A(2*i-1:2*i,i) =[1;-1];
end
b = [0.5;0.5;
    0.5;0.5;
    0.15;0.15;
    0.3;0.3;
    0.3;0.3];
barrierParam = 5*ones(1,length(b));


% Optimisation
tStart = tic;
out8 = optimise_shape_3D(myElementConstructor,nsetBC, ...
    nodes,elements,muscleBoundaries,kActu,U_8,h,tmax,A,b, ...
    'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, ...
    'maxIteration',50, ...
    'convCrit',0.004, ... 
    'convCritCost',1.0, ...
    'barrierParam',barrierParam, ...
    'gStepSize',0.0003, ... 
    'nRebuild',8, ...
    'rebuildThreshold',0.15,...
    'wSize', 5, ...
    'USEJULIA',1);
out8.tOpti = toc(tStart)/60;   % unit is minute
out8.tPerIt = out8.tOpti/out8.nIt;

% Print stats
fprintf('Computation time: %.2fmin\n',out8.tOpti)
fprintf('Number of built models and solved EoMs: %3d\n',out8.nIt)
fprintf('Computation time per models/EoMs: %.2f\n',out8.tPerIt)

% Save results
filename = sprintf('Results/Data/Appendix/Additional_shape_optimization/SO8_%d_el_kActu_%.3f.mat', ...
                   n_elements, kActu);
save(filename,'out8')




%% OPTIMISATION SO9 _______________________________________________________

% Preparation (shape variation and opt. constraints)
U_9 = [z_tail,z_tail_1,z_tail_2, y_linLongTail, y_ellipseFish, z_head];
nParam = 6;
A = zeros(2 * nParam, nParam);
for i = 1:nParam
    A(2*i-1:2*i,i) =[1;-1];
end
A = [A;
    1 1 1 0 0 0;
    -1 -1 -1 0 0 0];
b = [0.5;0.5;
    0.5;0.5;
    0.5;0.5;
    0.4;0.4;
    0.4;0.4;
    0.5;0.5;
    0.8;0.8];
barrierParam = 2*ones(1,length(b));
barrierParam(1, end-1:end) = 0.05;

% Optimisation
tStart = tic;
out9 = optimise_shape_3D(myElementConstructor,nsetBC, ...
    nodes,elements,muscleBoundaries,kActu,U_9,h,tmax,A,b,...
    'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, ...
    'maxIteration',50, ... 
    'convCrit',0.01, ...   
    'convCritCost',0.5, ... 
    'barrierParam',barrierParam, ...
    'gStepSize',0.0003,...   
    'nRebuild',5, ...
    'rebuildThreshold',0.15,...
    'wSize', 5, ...
    'USEJULIA',1);
out9.tOpti = toc(tStart)/60;   % unit is minute
out9.tPerIt = out9.tOpti/out9.nIt;

% Print stats
fprintf('Computation time: %.2fmin\n',out9.tOpti)
fprintf('Number of built models and solved EoMs: %3d\n',out9.nIt)
fprintf('Computation time per models/EoMs: %.2f\n',out9.tPerIt)

% Save results
filename = sprintf('Results/Data/Appendix/Additional_shape_optimization/SO9_%d_el_kActu_%.3f.mat', ...
                   n_elements, kActu);
save(filename,'out9')


%% RETRIEVE RESULTS _______________________________________________________
filename = sprintf('Results/Data/Appendix/Additional_shape_optimization/SO7_%d_el_kActu_%.3f.mat', ...
                   n_elements, kActu);
load(filename)
filename = sprintf('Results/Data/Appendix/Additional_shape_optimization/SO8_%d_el_kActu_%.3f.mat', ...
                   n_elements, kActu);
load(filename)

filename = sprintf('Results/Data/Appendix/Additional_shape_optimization/SO9_%d_el_kActu_%.3f.mat', ...
                   n_elements, kActu);
load(filename)

%% PLOT OPTIMAL SHAPE _____________________________________________________

SOIdx =   9;
switch SOIdx
    case 7
        xiStar = out7.xiStar;
        U = U_7;
    case 8
        xiStar = out8.xiStar;
        U = U_8;
    case 9
        xiStar = out9.xiStar;
        U = U_9;
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

fig_title = sprintf('Results/Figures/Appendix/Additional_shape_optimization/SO%d_opt_shape_%d_el_kActu_%.3f', SOIdx, n_elements, kActu);
exportgraphics(f_opt_shape,strcat(fig_title,'.pdf'),'Resolution',1200)

%% PRINT RESULTS AND PERFORMANCE STATS ____________________________________
SOIdx = 7;
switch SOIdx
    case 7
        xiStar = out7.xiStar;
        tOpti = out7.tOpti;
        nIt = out7.nIt;
        tPerIt = out7.tPerIt;    
    case 9
        xiStar = out9.xiStar;
        tOpti = out9.tOpti;
        nIt = out9.nIt;
        tPerIt = out9.tPerIt;   

end
fprintf('Computation time: %.2f min\n', tOpti);
fprintf('Number of built models and solved EoMs: %d\n', nIt);
fprintf('Computation time per built models/EoMs: %.2f min\n', tPerIt)
disp('Optimal parameters:')
disp(xiStar)
%% PLOT COST FUNCTION FOR THE 3 EXPERIMENTS _______________________________
f_cost = figure('units','centimeters','position',[3 3 9 6]);
hold on
plot(out7.LwoBEvo,'LineWidth',1)
plot(out8.LwoBEvo,'--','LineWidth',1)
grid on
ylabel('$$L$$')
xlabel('Iterations')
legend('SO7','SO8')
hold off

fig_title = sprintf('Results/Figures/Appendix/Additinal_shape_optimization/SO_add_cost_evolution_%d_el.pdf', n_elements);
exportgraphics(f_cost,fig_title,'Resolution',1200)

%% PLOT SWIMMING DISTANCE FOR THE 3 EXPERIMENTS ___________________________
f_dist = figure('units','centimeters','position',[3 3 9 6]);
hold on
p1 = plot(out7.xEvo*100,'LineWidth',1);
p2 = plot(out8.xEvo*100,'--','LineWidth',1);
grid on
ylabel('Swimming distance [cm]')
xlabel('Iterations')
legend([p1,p2],'SO7','SO8', 'Location', 'southeast')
hold off

fig_title = sprintf('Results/Figures/Appendix/Additional_shape_optimization/SO_add_dist_evolution_%d_el.pdf', n_elements);
exportgraphics(f_dist,fig_title,'Resolution',1200)

