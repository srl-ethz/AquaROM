% ------------------------------------------------------------------------ 
% B_shape_optimisation.m
%
% Description: 3D shape optimisation of a fish, to obtain best swimming
% performance.
% 
% Last modified: 17/05/2025, Mathieu Dubied, ETH Zurich
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
filename ='InputFiles/3d_rectangle_8086el';     % or 16009el, 8086

[MeshNominal, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
[Lx, Ly, Lz] = mesh_dimensions_3D(nodes);
n_elements = size(elements,1);

%% SHAPE VARIATIONS _______________________________________________________

[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish, xz_concaveTail] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);

% plot an example, with SO1 parameters
U = [z_tail,z_head,y_thinFish];
xiPlot = [-0.5;0.5,;0.5];
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
b = [0.5;0.5;0.5;0.5;0.3;0.3];
barrierParam = 10*ones(1,length(b));

% Optimisation
tStart = tic;
[xiStar_1,xiEvo_1,LEvo_1,LwoBEvo_1,xEvo_1,xEnergyEvo_1,volEvo_1,nIt_1,rebuildIdx_1] = ...
    optimise_shape_3D(myElementConstructor,nsetBC, ...
    nodes,elements,muscleBoundaries,kActu,U_1,h,tmax,A,b, ...
    'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, ...
    'maxIteration',40, ...
    'convCrit',0.01, ... % 0.004 for 8086, 0.01 for 16009
    'convCritCost',2.0, ... % 1.0 for 8086, 2.0 for 16009
    'barrierParam',barrierParam, ...
    'gStepSize',0.0002, ...  % 0.0004 for 8086, 0.0002 for 16609
    'nRebuild',6, ...
    'rebuildThreshold',0.15,...
    'wSize', 5, ...
    'USEJULIA',1);
tOpti_1 = toc(tStart)/60;   % unit is minute
tPerIt_1 = tOpti_1/nIt_1;

% Print stats
fprintf('Computation time: %.2fmin\n',tOpti_1)
fprintf('Number of built models and solved EoMs: %3d\n',nIt_1)
fprintf('Computation time per models/EoMs: %.2f\n',tPerIt_1)

% Save results
filename = sprintf('Results/Data/SO1_4_results_%d_el_kActu_%.3f.mat', n_elements, kActu);
save(filename,'xiStar_1','xiEvo_1','LEvo_1','LwoBEvo_1','tOpti_1', ...
                'xEvo_1','xEnergyEvo_1', 'volEvo_1','nIt_1','rebuildIdx_1','tPerIt_1')

%% OPTIMISATION SO2 _______________________________________________________

% Preparation (shape variation and opt. constraints)
U_2 = [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish];
nParam = 5;
A = zeros(2 * nParam, nParam);
for i = 1:nParam
    A(2*i-1:2*i,i) =[1;-1];
end
yTotConstr = [0 0 1 0 1;0 0  -1 0 -1;
              0 0 0 1 1;0 0  0 -1 -1];
A = [A;yTotConstr];
b = [0.5;0.5;
    0.5;0.5;
    0.4;0.4;
    0.5;0.5;
    0.5;0.5;   
    0.8;0.8;
    0.8;0.8];
barrierParam = ones(1,length(b));

% Optimisation
tStart = tic;
[xiStar_2,xiEvo_2,LEvo_2,LwoBEvo_2,xEvo_2,xEnergyEvo_2,volEvo_2,nIt_2,rebuildIdx_2] = ...
    optimise_shape_3D(myElementConstructor,nsetBC, ...
    nodes,elements,muscleBoundaries,kActu,U_2,h,tmax,A,b,...
    'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, ...
    'maxIteration',50, ... 
    'convCrit',0.01, ...   % set to 0.01 for 8086el, 0.01 for 16609
    'convCritCost',0.5, ... % set to 0.5 for 4270 el, 0.5 for 16609
    'barrierParam',barrierParam, ...
    'gStepSize',0.0002,...   % set to 0.0005 for 8086, 0.0002 for 16609
    'nRebuild',5, ...
    'rebuildThreshold',0.15,...
    'wSize', 5, ...
    'USEJULIA',1);
tOpti_2 = toc(tStart)/60;   % unit: minutes
tPerIt_2 = tOpti_2/nIt_2;  

% Print stats
fprintf('Computation time: %.2fmin\n',tOpti_2)
fprintf('Number of built models and solved EoMs: %5d\n',nIt_2)
fprintf('Computation time per models/EoMs: %.2f\n',tPerIt_2)

% Save results
filename = sprintf('Results/Data/SO2_5_results_%d_el_kActu_%.3f.mat', n_elements, kActu);
save(filename,'xiStar_2','xiEvo_2','LEvo_2','LwoBEvo_2','tOpti_2', ...
                'xEvo_2','xEnergyEvo_2','volEvo_2','nIt_2','rebuildIdx_2','tPerIt_2')

%% OPTIMISATION SO3 _______________________________________________________

% Preparation (shape variation and opt. constraints)
U_3 = [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish,...
    z_smallFish, z_notch, xz_concaveTail];

nParam = 8;
A = zeros(2 * nParam, nParam);
for i = 1:nParam
    A(2*i-1:2*i,i) =[1;-1];
end
b = [0.5;0.5;
    0.5;0.5;
    0.4;0.4;  
    0.3;0.3;
    0.4;0.4;    
    0.3;0.3;
    0.3;0.3;
    0.2; 0.01];  % concave tail only in one direction 
barrierParam = 1*ones(1,length(b)); % 3 for 8086, 1  for 16009

% Optimisation
tStart = tic;
[xiStar_3,xiEvo_3,LEvo_3,LwoBEvo_3,xEvo_3,xEnergyEvo_3,volEvo_3,nIt_3,rebuildIdx_3] = ...
    optimise_shape_3D(myElementConstructor,nsetBC, ...
    nodes,elements,muscleBoundaries,kActu,U_3,h,tmax,A,b, ...
    'FORMULATION',FORMULATION, ...
    'VOLUME',VOLUME, ...
    'maxIteration',40, ...
    'convCrit',0.015, ...    % set to 0.015 for 8086el, 0.015 for 16609
    'convCritCost',8.0, ...% set to 5.0 for 8086, 8.0 for 16609
    'barrierParam',barrierParam, ...
    'gStepSize',0.0003, ... % set to 0.0003 for 8086, 0.0003 for 16609
    'nRebuild',5, ...
    'rebuildThreshold',0.15, ...
    'wSize', 5, ...
    'USEJULIA',1);
tOpti_3 = toc(tStart)/60; % unit: minutes
tPerIt_3 = tOpti_3/nIt_3;
xiStar = xiStar_3;  

% Print stats
fprintf('Computation time: %.2fmin\n',tOpti_3)
fprintf('Number of built models and solved EoMs: %5d\n',nIt_3)
fprintf('Computation time per models/EoMs: %.2f\n',tPerIt_3)

% Save results
filename = sprintf('Results/Data/SO3_6_results_%d_el_kActu_%.3f.mat', n_elements, kActu);
save(filename,'xiStar_3','xiEvo_3','LEvo_3','LwoBEvo_3','tOpti_3', ...
                'xEvo_3','xEnergyEvo_3','volEvo_3','nIt_3','rebuildIdx_3','tPerIt_3')

%% RETRIEVE RESULTS _______________________________________________________
filename = sprintf('Results/Data/SO1_4_results_%d_el_kActu_%.3f.mat', n_elements, kActu);
load(filename)
filename = sprintf('Results/Data/SO2_5_results_%d_el_kActu_%.3f.mat', n_elements, kActu);
load(filename)
filename = sprintf('Results/Data/SO3_6_results_%d_el_kActu_%.3f.mat', n_elements, kActu);
load(filename)
%% PLOT OPTIMAL SHAPE _____________________________________________________

SOIdx = 1;
switch SOIdx
    case 1
        xiStar = xiStar_1;
        U = U_1; 
    case 2
        xiStar = xiStar_2;
        U = U_2;
    case 3
        xiStar = xiStar_3;
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
exportgraphics(f_opt_shape,strcat(fig_title,'.pdf'),'Resolution',1200)

%% PRINT RESULTS AND PERFORMANCE STATS ____________________________________
SOIdx = 3;
switch SOIdx
    case 1
        xiStar = xiStar_1;
        tOpti = tOpti_1;
        nIt = nIt_1;
        tPerIt = tPerIt_1;    
    case 2
        xiStar = xiStar_2;
        tOpti = tOpti_2;
        nIt = nIt_2;
        tPerIt = tPerIt_2; 
    case 3
        xiStar = xiStar_3;
        tOpti = tOpti_3;
        nIt = nIt_3;
        tPerIt = tPerIt_3; 
end
fprintf('Computation time: %.2f min\n', tOpti);
fprintf('Number of built models and solved EoMs: %d\n', nIt);
fprintf('Computation time per built models/EoMs: %.2f min\n', tPerIt)
disp('Optimal parameters:')
disp(xiStar)
%% PLOT COST FUNCTION FOR THE 3 EXPERIMENTS _______________________________
f_cost = figure('units','centimeters','position',[3 3 9 6]);
hold on
plot(LwoBEvo_1,'LineWidth',1)
plot(LwoBEvo_2,'--','LineWidth',1)
plot(LwoBEvo_3,'-.','LineWidth',1)
grid on
ylabel('$$L$$')
xlabel('Iterations')
legend('SO1','SO2','SO3')
hold off

fig_title = sprintf('Results/Figures/SO_cost_evolution_%d_el.pdf', n_elements);
exportgraphics(f_cost,fig_title,'Resolution',1200)

%% PLOT SWIMMING DISTANCE FOR THE 3 EXPERIMENTS ___________________________
f_dist = figure('units','centimeters','position',[3 3 9 6]);
hold on
p1 = plot(xEvo_1*100,'LineWidth',1);
% plot(rebuildIdx_1, xEvo_1(rebuildIdx_1) * 100, 'x', ...
%      'MarkerSize', 8, 'LineWidth', 1.5, 'Color', p1.Color);
p2 = plot(xEvo_2*100,'--','LineWidth',1);
% plot(rebuildIdx_2, xEvo_2(rebuildIdx_2) * 100, 'x', ...
%      'MarkerSize', 8, 'LineWidth', 1.5, 'Color', p2.Color);
p3 = plot(xEvo_3*100,'-.','LineWidth',1);
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
p1 = plot(xEnergyEvo_1*100,'LineWidth',1);
% plot(rebuildIdx_1+1, xEnergyEvo_1(rebuildIdx_1+1) * 100, 'x', ...
%      'MarkerSize', 8, 'LineWidth', 1.0, 'Color', p1.Color);
p2 = plot(xEnergyEvo_2*100,'--','LineWidth',1);
% plot(rebuildIdx_2+1, xEnergyEvo_2(rebuildIdx_2+1) * 100, 'x', ...
%      'MarkerSize', 8, 'LineWidth', 1.0, 'Color', p2.Color)
p3 = plot(xEnergyEvo_3*100,'-.','LineWidth',1);
% plot(rebuildIdx_3+1, xEnergyEvo_3(rebuildIdx_3+1) * 100, 'x', ...
%      'MarkerSize', 8, 'LineWidth', 1.0, 'Color', p3.Color)
grid on
ylabel('Normalized swimming distance')
xlabel('Iterations')
legend([p1,p2,p3], 'SO1','SO2','SO3', 'Location', 'northwest')
hold off

fig_title = sprintf('Results/Figures/SO_normalized_dist_evolution_%d_el.pdf', n_elements);
exportgraphics(f_norm_dist,fig_title,'Resolution',1200)




%% ------------------------------------------------------------------------










%% SIMULATION/ANIMATION
U_sim = U_3;
xiStar_3=[-0.49;0.26;0.4;0.28;0.39;0.3;0.29;0.15];
% Find dorsal nodes on nominal mesh
[spineNodes, spineElements, spineElementWeights, nodeIdxPosInElements] = find_spine_TET4(elements,nodes);
[~,matchedDorsalNodesIdx,dorsalNodesElementsVec,matchedDorsalNodesZPos] = ....
            find_dorsal_nodes(elements, nodes, spineElements, nodeIdxPosInElements);

dorsalNodesStructFromUser.matchedDorsalNodesIdx = matchedDorsalNodesIdx;
dorsalNodesStructFromUser.dorsalNodesElementsVec = dorsalNodesElementsVec;
dorsalNodesStructFromUser.matchedDorsalNodesZPos = matchedDorsalNodesZPos;

% shape varied mesh
df = U_sim*xiStar_3;                       % displacement fields introduced by defects
ddf = [df(1:3:end) df(2:3:end) df(3:3:end)]; 
nodes_defected = nodes + ddf;    % nominal + d ---> defected 
svMesh = Mesh(nodes_defected);
svMesh.create_elements_table(elements,myElementConstructor);
for l=1:length(nsetBC)
    svMesh.set_essential_boundary_condition([nsetBC{l}],1:3,0)   
end


% Build PROM
[V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuLeft,actuRight] = ...
                 build_PROM_3D(svMesh,nodes,elements,muscleBoundaries,U_sim,USEJULIA,VOLUME,FORMULATION,'dorsalNodes',dorsalNodesStructFromUser);


             %%
% Simulate PROM
TI_NL_PROM = solve_EoMs_and_sensitivities(V, PROM_Assembly, tensors_PROM, tailProperties, spineProperties, dragProperties, actuLeft, actuRight, kActu, h, 8.0); 

%% PLOT RESULTS
tmax=8.0;
uTail = zeros(3,tmax/h);
timePlot = linspace(0,tmax-h,tmax/h);
x0Tail = min(nodes(:,1));
for a=1:tmax/h
    uTail(:,a) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*TI_NL_PROM.Solution.q(:,a);
end

figure
subplot(2,1,1);
plot(timePlot,x0Tail+uTail(1,:),'DisplayName','k=0')
hold on
xlabel('Time [s]')
ylabel('x-position tail node')
legend('Location','northwest')
subplot(2,1,2);
plot(timePlot,uTail(2,:))
hold on
xlabel('Time [s]')
ylabel('y-position tail node')
legend('Location','southwest')
hold off

%% ANIMATION

elementPlot = elements(:,1:4); 
nel = size(elements,1);
% top muscle
topMuscle = zeros(nel,1);

for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
    if elementCenterY>0.00 &&  elementCenterX < -Lx*propRigid && elementCenterX > -Lx*0.9
        topMuscle(el) = 1;
    end    
end

% bottom muscle
bottomMuscle = zeros(nel,1);
for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
    if elementCenterY<0.00 &&  elementCenterX < -Lx*propRigid && elementCenterX > -Lx*0.9
        bottomMuscle(el) = 1;
    end    

end

actuationValues = zeros(size(TI_NL_PROM.Solution.q,2),1);
for t=1:size(TI_NL_PROM.Solution.q,2)
    actuationValues(t) = 0;
end

actuationValues2 = zeros(size(TI_NL_PROM.Solution.q,2),1);
for t=1:size(TI_NL_PROM.Solution.q,2)
    actuationValues2(t) = 0;
end


sol = V*TI_NL_PROM.Solution.q(:,1:end);
AnimateFieldonDeformedMeshActuation2Muscles(nodes_defected, elementPlot,topMuscle,actuationValues,...
    bottomMuscle,actuationValues2,sol, ...
    'factor',1,'index',1:3,'filename','result_video','framerate',1/h)


%% OTHER/PREVIOUS PLOTS ___________________________________________________








%% PLOT PARAMETERS' EVOLUTION OVER ITERATIONS _____________________________
figure('Position',[100,100,600,200])
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
subplot(1,3,1)
plot(xiEvo_3(5,:));
ylabel('$$\xi_1$$','Interpreter','latex')
xlabel('Iterations')
grid on
subplot(1,3,2)
plot(xiEvo(2,:));
ylabel('$$\xi_2$$','Interpreter','latex')
xlabel('Iterations')
grid on
subplot(1,3,3)
plot(xiEvo(3,:));
ylabel('$$\xi_3$$','Interpreter','latex')
xlabel('Iterations')
grid on
