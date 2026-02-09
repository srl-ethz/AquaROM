% ------------------------------------------------------------------------ 
% C_co_optimisation.m
%
% Description:  optimization of both shape and actuation signal. Use of a
%               multiobjective costs to trade energy vs. swimming distance.
% 
% Last modified: 09/02/2026, Mathieu Dubied, ETH Zurich
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

% plot an example, with SO1 parameters
U = [z_tail,z_head,y_thinFish];
xiPlot = [-0.5;0.5,;-0.5];
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
kActu = 6.0*1e4;

%% CO-OPTIMISATION  _______________________________________________________

% Preparation (shape variation and opt. constraints). Need to provide U as
% we solve a PROM. We only use the gradient from the actuation however.
% Preparation (shape variation and opt. constraints)
U_2 = [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish];
nParam = 6;
A = zeros(2 * nParam, nParam);
for i = 1:nParam
    A(2*i-1:2*i,i) =[1;-1];
end
b = [0.5;0.5;
    0.4;0.4;
    0.4;0.4;
    0.5;0.5;
    0.5;0.5;
    0.40;-0.05];
barrierParam = 1*ones(1,length(b));

w1 = 0.98;   % weight on distance objective
w2 = 1-w1;  % weight on energy objective

% Optimisation
tStart = tic;
out = optimise_shape_actuation(myElementConstructor,nsetBC, ...
    nodes,elements,muscleBoundaries,kActu,U_2,h,tmax,A,b, ...
    'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, ...
    'maxIteration',50, ...
    'convCrit',0.002, ... 
    'convCritCost',1.0, ... 
    'barrierParam',barrierParam, ...
    'gStepSize',0.0006, ...  % 0.0005 ok
    'nRebuild',15, ...
    'rebuildThreshold',0.3,...
    'wSize', 5, ...
    'USEJULIA',1, ...
    'paramInit', [0.2], ... % for actuation signal (amplitude)
    'w1', w1, ...
    'w2', w2,...
    'alphaActu', 800, ...
    'rebuildThresholdActu', 0.05, ...
    'learningRateActu', 0.006);%0.001    
out.tOpti = toc(tStart)/60;   % unit is minute
out.tPerIt = out.tOpti/out.nIt;

% Print stats
fprintf('Computation time: %.2fmin\n',out.tOpti)
fprintf('Number of built models and solved EoMs: %3d\n',out.nIt)
fprintf('Computation time per models/EoMs: %.2f\n',out.tPerIt)

% Save results
filename = sprintf('Results/Data/C_co_optimization/C_w1_%.2f_w2_%.2f.mat', ...
                   w1, w2);
save(filename,'out')

%% RESULTS ANALYSIS _______________________________________________________
w1 = 0.5;
filename = sprintf('Results/Data/C_co_optimization/C_w1_%.2f_w2_%.2f.mat', ...
                   w1, 1-w1);
load(filename)
out_1 = out;

w1 = 0.7;
filename = sprintf('Results/Data/C_co_optimization/C_w1_%.2f_w2_%.2f.mat', ...
                    w1, 1-w1);
load(filename)
out_2 = out;

w1 = 0.2;
filename = sprintf('Results/Data/C_co_optimization/C_w1_%.2f_w2_%.2f.mat', ...
                    w1, 1-w1);
load(filename)
out_3 = out;

w1 = 0.9;
filename = sprintf('Results/Data/C_co_optimization/C_w1_%.2f_w2_%.2f.mat', ...
                    w1, 1-w1);
load(filename)
out_4 = out;

w1 = 0.98;
filename = sprintf('Results/Data/C_co_optimization/C_w1_%.2f_w2_%.2f.mat', ...
                    w1, 1-w1);
load(filename)
out_5 = out;


%%
% Put all trajectories in cell arrays
X = {out_1.LObj2Evo, out_2.LObj2Evo, out_3.LObj2Evo, out_5.LObj2Evo};   % add more if needed
Y = {out_1.LObj1Evo, out_2.LObj1Evo, out_3.LObj1Evo, out_5.LObj1Evo};


% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'linear')
f = figure('units','centimeters','position',[3 3 9 9]);
defaultColors = get(gca, 'ColorOrder');
hold on

for k = 1:numel(X)

    % plot trajectory and capture handle
    p = plot(X{k}, Y{k}, '-o', 'LineWidth', 1, 'Color', defaultColors(k,:));

    % use same color for markers
    c = p.Color;

    % end marker
    plot(X{k}(end), Y{k}(end), ...
        'o', 'MarkerFaceColor', c, ...
        'Color', c, ...
        'MarkerSize', 6, 'LineWidth', 2);

end

% start marker
plot(X{1}(1), Y{1}(1), ...
    'x', 'Color', 'k', ...
    'MarkerSize', 12, 'LineWidth', 2);

xlabel('$L_{energy}$','Interpreter','latex')
ylabel('$L_{distance}$','Interpreter','latex')
grid on
fig_title = sprintf('Results/Figures/C_co_optimization/cost_trajectories_co_optimization.pdf');
exportgraphics(f,fig_title,'Resolution',1200)





