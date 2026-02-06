% ------------------------------------------------------------------------ 
% B_appendix_PROM_local_validity.m
%
% Description:  Compare the PROM at deformed configuration with a ROM
%               built at these deformed configuration
% 
% Last modified: 06/02/2026, Mathieu Dubied, ETH Zurich
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

%% LOAD CASE SO1 __________________________________________________________                                                   

% load material parameters
load('parameters.mat') 

% parameters
n_elements = 8086;
kActu = 6.0*1e4;
h = 0.02;
tmax = 2.0;
filename = strcat('InputFiles/3d_rectangle_', num2str(n_elements), 'el');
base_data_output_dir = 'Results/Data/Appendix/PROM_local_validity/';
base_fig_output_dir = 'Results/Figures/Appendix/PROM_local_validity/';

% load optimization results
results_filename = sprintf('Results/Data/B_shape_optimization/SO1_8086_el_kActu_%.3f.mat',kActu);
load(results_filename)

%% PREPARE NOMINAL MODEL FOR REPRODUCIBILITY ______________________________
% nominal model
[MeshNominal, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
[Lx, Ly, Lz] = mesh_dimensions_3D(nodes);
volVector = compute_nominal_vol_per_element(MeshNominal,size(elements,1));
headNode = find_node(0, 0, 0, nodes);
for l=1:length(nsetBC)
    MeshNominal.set_essential_boundary_condition([nsetBC{l}],1:3,0)   
end

% shape variation
[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish, xz_concaveTail, y_fin, ...
    y_bumpBack, y_bumpFront] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);
U = [z_tail,z_head,y_thinFish];

% build PROM
fprintf('____________________\n')
fprintf('Building nominal PROM ... \n')

[V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
build_PROM_3D(MeshNominal,nodes,elements,muscleBoundaries,U,USEJULIA,VOLUME,FORMULATION, 'nomVolVector', volVector); 

% store dorsal nodes for future use
dorsalNodesStructFromUser.matchedDorsalNodesIdx = spineProperties.dorsalNodeIdx;
dorsalNodesStructFromUser.dorsalNodesElementsVec = spineProperties.dorsalNodesElementVec;
dorsalNodesStructFromUser.matchedDorsalNodesZPos = spineProperties.zPos;
tailNode = tailProperties.tailNode;

%% VISUALIZE OPTIMIZATION TRAJECTORY FOR PROM REBUILD _____________________
% figure of cost function
f_cost = figure('units','centimeters','position',[3 3 9 6]);
ax = gca;
hold on
y = out1.LwoBEvo(2:end);        % first index is the initialization
rebuildIdx = out1.rebuildIdx;
rebuildIdx = [1,rebuildIdx];    % add the first step for figure
plot(y,'LineWidth',1,'HandleVisibility','off')
plot(rebuildIdx, y(rebuildIdx), 'o', ...
    'MarkerSize',6, ...
    'MarkerFaceColor',ax.ColorOrder(1,:), ...
    'MarkerEdgeColor',ax.ColorOrder(1,:), ...
    'DisplayName','PROM (re)build');
grid on
ylabel('$$L$$')
xlabel('Iterations')
legend('show')
hold off

% select iteration step to analyse
rebuildIdx = rebuildIdx(2:4);
middleIdx = [3,6,10];

% parameters at selected rebuild idx
xiAtRebuild = out1.xiEvo(:,rebuildIdx);
xiRebuildAtRebuild = out1.xiRebuildEvo(:,rebuildIdx);

% parameters at selectedmiddle idx
xiAtMiddle = out1.xiEvo(:,middleIdx);
xiRebuildAtMiddle = out1.xiRebuildEvo(:,middleIdx);

% figure with cases to analysis
f_overview = figure('units','centimeters','position',[3 3 9 6]);
hold on
y = out1.LwoBEvo(2:end);        % first index is the initialization
plot(y,'LineWidth',1,'HandleVisibility','off')

% rebuild idx
plot(rebuildIdx, y(rebuildIdx), 'o', ...
    'MarkerSize',6, ...
    'MarkerFaceColor',ax.ColorOrder(1,:), ...
    'MarkerEdgeColor',ax.ColorOrder(1,:), ...
    'DisplayName','PROM (re)build');
letters = arrayfun(@(k) char('A'+k-1), 1:numel(rebuildIdx), ...
                   'UniformOutput', false);
text(rebuildIdx+0.5, y(rebuildIdx), letters, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'Color', ax.ColorOrder(1,:), ...
    'FontWeight','bold', ...
    'DisplayName','PROM (re)build');

% middle idx
plot(middleIdx, y(middleIdx), 'x', ...
    'MarkerSize',6, ...
    'MarkerFaceColor',ax.ColorOrder(2,:), ...
    'MarkerEdgeColor',ax.ColorOrder(2,:), ...
    'DisplayName','PROM (re)build');
letters = arrayfun(@(k) char('a'+k-1), 1:numel(middleIdx), ...
                   'UniformOutput', false);
text(middleIdx+0.5, y(middleIdx), letters, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'Color', ax.ColorOrder(2,:), ...
    'FontWeight','bold', ...
    'DisplayName','Middle configuration');

grid on
ylabel('$$L$$')
xlabel('Iterations')
legend('show')
hold off

%% RUN MODELS SEGMENT a-A _________________________________________________

xi_0 = out1.xiEvo(:,1);
xi_a = xiAtMiddle(:,1);
xi_A = xiAtRebuild(:,1);

% PROM used in optimization: xi=xi_0
[eta_0, S_0, V] = get_results_PROM(nodes, elements,muscleBoundaries, ...
                                    myElementConstructor, nsetBC, ...                            
                                    U, xi_0, 'nominal', ...
                                    kActu, tmax, h, ...
                                    USEJULIA,VOLUME,FORMULATION, ...
                                    dorsalNodesStructFromUser, ...
                                    volVector);
% Approximation for a)
eta_approx_a = eta_0 + double(ttv(S_0,xi_a,2));
[uHead_approx_a, uTail_approx_a] = compute_head_tail_motion(eta_approx_a, V, headNode, tailNode);

% Aproximation for A)
eta_approx_A = eta_0 + double(ttv(S_0,xi_A,2));             
[uHead_approx_A, uTail_approx_A] = compute_head_tail_motion(eta_approx_A, V, headNode, tailNode);

% Exact solution at a)
[eta_a, S_a, V] = get_results_PROM(nodes, elements,muscleBoundaries, ...
                                    myElementConstructor, nsetBC, ...                            
                                    U, xi_a, 'a)', ...
                                    kActu, tmax, h, ...
                                    USEJULIA,VOLUME,FORMULATION, ...
                                    dorsalNodesStructFromUser, ...
                                    volVector);
[uHead_a, uTail_a] = compute_head_tail_motion(eta_a, V, headNode, tailNode);                              

% Exact solution at A)
[eta_A, S_A, V] = get_results_PROM(nodes, elements,muscleBoundaries, ...
                                    myElementConstructor, nsetBC, ...                            
                                    U, xi_A, 'A)', ...
                                    kActu, tmax, h, ...
                                    USEJULIA,VOLUME,FORMULATION, ...
                                    dorsalNodesStructFromUser, ...
                                    volVector);
[uHead_A, uTail_A] = compute_head_tail_motion(eta_A, V, headNode, tailNode); 

%% Store results
outfile = strcat(base_data_output_dir, 'aA.mat');
outaA = wrap_results(uHead_approx_a, uTail_approx_a, S_0, ...
                     uHead_a, uTail_a, S_a, ...
                     uHead_approx_A, uTail_approx_A, ...
                     uHead_A, uTail_A, S_A);

save(outfile, 'outaA');

%% RUN MODELS SEGMENT b-B _________________________________________________

xi_A = xiAtRebuild(:,1);
xi_b = xiAtMiddle(:,2);
xi_B = xiAtRebuild(:,2);
xi_Rebuild_b = xiRebuildAtMiddle(:,2);
xi_Rebuild_B = xiRebuildAtRebuild(:,2);

% PROM used in optimization: xi=xi_0
[eta_A, S_A, V] = get_results_PROM(nodes, elements,muscleBoundaries, ...
                                    myElementConstructor, nsetBC, ...                            
                                    U, xi_A, 'nominal in A', ...
                                    kActu, tmax, h, ...
                                    USEJULIA,VOLUME,FORMULATION, ...
                                    dorsalNodesStructFromUser, ...
                                    volVector);
% Approximation for b)
eta_approx_b = eta_A + double(ttv(S_A,xi_Rebuild_b,2));
[uHead_approx_b, uTail_approx_b] = compute_head_tail_motion(eta_approx_b, V, headNode, tailNode);

% Aproximation for B)
eta_approx_B = eta_A + double(ttv(S_A,xi_Rebuild_B,2));             
[uHead_approx_B, uTail_approx_B] = compute_head_tail_motion(eta_approx_B, V, headNode, tailNode);

% Exact solution at b)
[eta_b, S_b, V] = get_results_PROM(nodes, elements,muscleBoundaries, ...
                                    myElementConstructor, nsetBC, ...                            
                                    U, xi_b, 'b)', ...
                                    kActu, tmax, h, ...
                                    USEJULIA,VOLUME,FORMULATION, ...
                                    dorsalNodesStructFromUser, ...
                                    volVector);
[uHead_b, uTail_b] = compute_head_tail_motion(eta_b, V, headNode, tailNode);                              

% Exact solution at B)
[eta_B, S_B, V] = get_results_PROM(nodes, elements,muscleBoundaries, ...
                                    myElementConstructor, nsetBC, ...                            
                                    U, xi_B, 'B)', ...
                                    kActu, tmax, h, ...
                                    USEJULIA,VOLUME,FORMULATION, ...
                                    dorsalNodesStructFromUser, ...
                                    volVector);
[uHead_B, uTail_B] = compute_head_tail_motion(eta_B, V, headNode, tailNode); 

% Store results
outfile = strcat(base_data_output_dir, 'bB.mat');
outbB = wrap_results(uHead_approx_b, uTail_approx_b, S_A, ...
                     uHead_b, uTail_b, S_b, ...
                     uHead_approx_B, uTail_approx_B, ...
                     uHead_B, uTail_B, S_B);

save(outfile, 'outbB');

%% ANALYSIS
% Load result
savedfile = strcat(base_data_output_dir, 'aA.mat');
load(savedfile);
savedfile = strcat(base_data_output_dir, 'bB.mat');
load(savedfile);

%% ANALYSIS RESULTS SEGMENT a-A ___________________________________________
timePlot = linspace(0,tmax-h,tmax/h+1);
% a)
fig = fig_comparison_PROM_approximation(timePlot, ...
                                outaA.uHead_1, outaA.uHead_approx_1, ...
                                outaA.uTail_1, outaA.uTail_approx_1, ...
                                Lx, Ly, [3 3 9 7.0]);
title('a')
% A)
fig = fig_comparison_PROM_approximation(timePlot, ...
                                outaA.uHead_2, outaA.uHead_approx_2, ...
                                outaA.uTail_2, outaA.uTail_approx_2, ...
                                Lx, Ly, [13 3 9 7.0]);
title('A')

%% ANALYSIS RESULTS SEGMENT b-B ___________________________________________
timePlot = linspace(0,tmax-h,tmax/h+1);
% b)
fig = fig_comparison_PROM_approximation(timePlot, ...
                                outbB.uHead_1, outbB.uHead_approx_1, ...
                                outbB.uTail_1, outbB.uTail_approx_1, ...
                                Lx, Ly, [3 3 9 7.0]);
title('b')
% B)
fig = fig_comparison_PROM_approximation(timePlot, ...
                                outbB.uHead_2, outbB.uHead_approx_2, ...
                                outbB.uTail_2, outbB.uTail_approx_2, ...
                                Lx, Ly, [13 3 9 7.0]);
title('B')


%% UTILS __________________________________________________________________

% Wrapper with eta and S as output. MeshNominal is the nominal mesh.
function [eta, S, V]= get_results_PROM(nodes,elements,muscleBoundaries, ...
                                        myElementConstructor, nsetBC, ...                            
                                        U, xi, name, ...
                                        kActu, tmax, h, ...
                                        USEJULIA,VOLUME,FORMULATION, ...
                                        dorsalNodesStructFromUser, ...
                                        volVector)                
                                    
    % update shape-varied mesh nodes
    df = U*xi;                       % displacement fields introduced by shape variations/defects
    ddf = [df(1:3:end) df(2:3:end) df(3:3:end)]; 
    nodes_defected = nodes + ddf;    % nominal + d ---> defected 
    svMesh = Mesh(nodes_defected);
    svMesh.create_elements_table(elements,myElementConstructor);
    for l=1:length(nsetBC)
        svMesh.set_essential_boundary_condition([nsetBC{l}],1:3,0)   
    end
    fprintf('____________________\n')
    fprintf('Building PROM for %s...\n', name) 
    [V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
         build_PROM_3D(svMesh,nodes_defected,elements,muscleBoundaries,U,USEJULIA,VOLUME,FORMULATION,...
                        'dorsalNodes',dorsalNodesStructFromUser, 'nomVolVector', volVector);

    fprintf('____________________\n')
    fprintf('Solving EoMs and sensitivity for %s ...\n', name) 
    TI_NL_PROM = solve_EoMs_and_sensitivities(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,kActu,h,tmax);                        
    eta = TI_NL_PROM.Solution.q;
    S = TI_NL_PROM.Solution.s;
    S=tensor(S);
                      
end

% Compute uHead and uTail
function [uHead, uTail] = compute_head_tail_motion(eta, V, headNode, tailNode)
    nStep = size(eta,2);
    for t = 1:nStep              
        uHead(:, t) = V(headNode*3-2:headNode*3, :)*eta(:,t)*100;
        uTail(:,t) = V(tailNode*3-2:tailNode*3,:)*eta(:,t)*100;
    end
end

% Results wrapper
function out = wrap_results(uHead_approx_1, uTail_approx_1, S_0, ...
                            uHead_1, uTail_1, S_1, ...
                            uHead_approx_2, uTail_approx_2, ...
                            uHead_2, uTail_2, S_2)
    out = struct();
    out.uHead_approx_1 = uHead_approx_1;
    out.uTail_approx_1 = uTail_approx_1;
    out.S_0            = S_0;

    out.uHead_1 = uHead_1;
    out.uTail_1 = uTail_1;
    out.S_1     = S_1;

    out.uHead_approx_2 = uHead_approx_2;
    out.uTail_approx_2 = uTail_approx_2;

    out.uHead_2 = uHead_2;
    out.uTail_2 = uTail_2;
    out.S_2     = S_2;

end
                   
