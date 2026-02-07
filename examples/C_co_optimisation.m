% ------------------------------------------------------------------------ 
% C_co_optimisation.m
%
% Description:  optimization of both shape and actuation signal. Use of a
%               multiobjective costs to trade energy vs. swimming distance.
% 
% Last modified: 07/02/2026, Mathieu Dubied, ETH Zurich
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
kActu = 3.1*1e5;    % multiplicative factor for the actuation forces (8086 elements)
kActu = 5.0*1e4;   %7.4
% kActu = 2.8*1e5;       % value for 16009 elements


%% OPTIMISATION TEST 1 ____________________________________________________

% Preparation (shape variation and opt. constraints). Need to provide U as
% we solve a PROM. We only use the gradient from the actuation however.
U_1 = [z_tail,z_head,y_thinFish]; 
nParam = 4; % 3 shape parameters and 1 actuation paramter
A = zeros(2 * nParam, nParam);
for i = 1:nParam
    A(2*i-1:2*i,i) =[1;-1];
end
b = [0.5;0.5;
    0.5;0.5;
    0.15;0.15;
    0.4;0];
barrierParam = 5*ones(1,length(b));

% Optimisation
tStart = tic;
out = optimise_shape_actuation(myElementConstructor,nsetBC, ...
    nodes,elements,muscleBoundaries,kActu,U_1,h,tmax,A,b, ...
    'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, ...
    'maxIteration',40, ...
    'convCrit',0.002, ... 
    'convCritCost',1.0, ... 
    'barrierParam',barrierParam, ...
    'gStepSize',0.0001, ...  
    'nRebuild',6, ...
    'rebuildThreshold',0.05,...
    'wSize', 5, ...
    'USEJULIA',1, ...
    'paramInit', [0.2]);    % for actuation signal (amplitude)
out.tOpti = toc(tStart)/60;   % unit is minute
out.tPerIt = out.tOpti/out.nIt;

% Print stats
fprintf('Computation time: %.2fmin\n',out.tOpti)
fprintf('Number of built models and solved EoMs: %3d\n',out.nIt)
fprintf('Computation time per models/EoMs: %.2f\n',out.tPerIt)

% Save results
filename = sprintf('Results/Data/C_actuation_optimization/actuation_only_%d_el_kActu_%.3f.mat', ...
                   n_elements, kActu);
save(filename,'out')

%% RESULTS ANALYSIS _______________________________________________________
filename = sprintf('Results/Data/C_actuation_optimization/actuation_only_%d_el_kActu_%.3f.mat', ...
                   n_elements, kActu);
load(filename)

f_dist = figure('units','centimeters','position',[3 3 9 6]);
hold on
p1 = plot(out.xEvo*100,'LineWidth',1);
grid on
ylabel('Swimming distance [cm]')
xlabel('Iterations')
legend([p1],'AO1', 'Location', 'southeast')
hold off



