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
results_filename = sprintf('Results/Data/B_shape_optimization/SO2_8086_el_kActu_%.3f.mat',kActu);
load(results_filename)

%% PREPARE NOMINAL MODEL FOR REPRODUCIBILITY ______________________________
% nominal model
[MeshNominal, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
[Lx, Ly, Lz] = mesh_dimensions_3D(nodes);

%%
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
U = [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish];

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
y = out2.LwoBEvo(2:end);        % first index is the initialization
rebuildIdx = out2.rebuildIdx;
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
rebuildIdx = rebuildIdx(1:end);
middleIdx = [3,7,10, 13];

% parameters at selected rebuild idx
xiAtRebuild = out2.xiEvo(:,rebuildIdx);
xiRebuildAtRebuild = out2.xiRebuildEvo(:,rebuildIdx);

% parameters at selectedmiddle idx
xiAtMiddle = out2.xiEvo(:,middleIdx);
xiRebuildAtMiddle = out2.xiRebuildEvo(:,middleIdx);

% figure with cases to analysis
f_overview = figure('units','centimeters','position',[3 3 9 6]);
hold on
y = out2.LwoBEvo(2:end);        % first index is the initialization
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
    'DisplayName','midsteps');
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

fig_title = sprintf('Results/Figures/Appendix/PROM_local_validity/cost_function_with_markers.pdf');
exportgraphics(f_cost,fig_title,'Resolution',1200)

%% RUN MODELS SEGMENT A-a-B _________________________________________________

xi_A = xiAtRebuild(:,1);
xi_a = xiAtMiddle(:,1);
xi_B = xiAtRebuild(:,2);

% PROM used in optimization: xi=xi_0
[eta_A, S_A, V_A] = get_results_PROM(nodes, elements,muscleBoundaries, ...
                                    myElementConstructor, nsetBC, ...                            
                                    U, xi_A, 'nominal', ...
                                    kActu, tmax, h, ...
                                    USEJULIA,VOLUME,FORMULATION, ...
                                    dorsalNodesStructFromUser, ...
                                    volVector);
% Approximation for a)
eta_approx_a = eta_A + double(ttv(S_A,xi_a,2));
[uHead_approx_a, uTail_approx_a] = compute_head_tail_motion(eta_approx_a, V_A, headNode, tailNode);

% Aproximation for B)
eta_approx_B = eta_A + double(ttv(S_A,xi_B,2));             
[uHead_approx_B, uTail_approx_B] = compute_head_tail_motion(eta_approx_B, V_A, headNode, tailNode);

% Exact solution at a)
[eta_a, S_a, V_a] = get_results_PROM(nodes, elements,muscleBoundaries, ...
                                    myElementConstructor, nsetBC, ...                            
                                    U, xi_a, 'a)', ...
                                    kActu, tmax, h, ...
                                    USEJULIA,VOLUME,FORMULATION, ...
                                    dorsalNodesStructFromUser, ...
                                    volVector);
[uHead_a, uTail_a] = compute_head_tail_motion(eta_a, V_a, headNode, tailNode);                              

% Exact solution at B)
[eta_B, S_B, V_B] = get_results_PROM(nodes, elements,muscleBoundaries, ...
                                    myElementConstructor, nsetBC, ...                            
                                    U, xi_B, 'B)', ...
                                    kActu, tmax, h, ...
                                    USEJULIA,VOLUME,FORMULATION, ...
                                    dorsalNodesStructFromUser, ...
                                    volVector);
[uHead_B, uTail_B] = compute_head_tail_motion(eta_B, V_B, headNode, tailNode); 

%% Store results
outfile = strcat(base_data_output_dir, 'AaB.mat');
outAaB = wrap_results(uHead_approx_a, uTail_approx_a, S_A, ...
                     uHead_a, uTail_a, S_a, ...
                     uHead_approx_B, uTail_approx_B, ...
                     uHead_B, uTail_B, S_B);

save(outfile, 'outAaB');

%% RUN MODELS SEGMENT B-b-C (B already ok)_________________________________

xi_b = xiAtMiddle(:,2);
xi_C = xiAtRebuild(:,3);
xi_Rebuild_b = xiRebuildAtMiddle(:,2);
xi_Rebuild_C = xiRebuildAtRebuild(:,2);

% Approximation for b)
eta_approx_b = eta_B + double(ttv(S_B,xi_Rebuild_b,2));
[uHead_approx_b, uTail_approx_b] = compute_head_tail_motion(eta_approx_b, V_B, headNode, tailNode);

% Aproximation for C)
eta_approx_C = eta_B + double(ttv(S_B,xi_Rebuild_B,2));             
[uHead_approx_C, uTail_approx_C] = compute_head_tail_motion(eta_approx_C, V_B, headNode, tailNode);

% Exact solution at b)
[eta_b, S_b, V_b] = get_results_PROM(nodes, elements,muscleBoundaries, ...
                                    myElementConstructor, nsetBC, ...                            
                                    U, xi_b, 'b)', ...
                                    kActu, tmax, h, ...
                                    USEJULIA,VOLUME,FORMULATION, ...
                                    dorsalNodesStructFromUser, ...
                                    volVector);
[uHead_b, uTail_b] = compute_head_tail_motion(eta_b, V_b, headNode, tailNode);                              

% Exact solution at C)
[eta_C, S_C, V_C] = get_results_PROM(nodes, elements,muscleBoundaries, ...
                                    myElementConstructor, nsetBC, ...                            
                                    U, xi_C, 'C)', ...
                                    kActu, tmax, h, ...
                                    USEJULIA,VOLUME,FORMULATION, ...
                                    dorsalNodesStructFromUser, ...
                                    volVector);
[uHead_C, uTail_C] = compute_head_tail_motion(eta_C, V_C, headNode, tailNode); 

% Store results
outfile = strcat(base_data_output_dir, 'BbC.mat');
outBbC = wrap_results(uHead_approx_b, uTail_approx_b, S_B, ...
                     uHead_b, uTail_b, S_b, ...
                     uHead_approx_C, uTail_approx_C, ...
                     uHead_C, uTail_C, S_C);

save(outfile, 'outBbC');

%% RUN MODELS SEGMENT C-c-D (C already ok)_________________________________

xi_c = xiAtMiddle(:,3);
xi_D = xiAtRebuild(:,4);
xi_Rebuild_c = xiRebuildAtMiddle(:,3);
xi_Rebuild_D = xiRebuildAtRebuild(:,3);

% Approximation for c)
eta_approx_c = eta_C + double(ttv(S_C,xi_Rebuild_c,2));
[uHead_approx_c, uTail_approx_c] = compute_head_tail_motion(eta_approx_c, V_C, headNode, tailNode);

% Aproximation for D)
eta_approx_D = eta_C + double(ttv(S_C,xi_Rebuild_D,2));             
[uHead_approx_D, uTail_approx_D] = compute_head_tail_motion(eta_approx_D, V_C, headNode, tailNode);

% Exact solution at c)
[eta_c, S_c, V_c] = get_results_PROM(nodes, elements,muscleBoundaries, ...
                                    myElementConstructor, nsetBC, ...                            
                                    U, xi_c, 'c)', ...
                                    kActu, tmax, h, ...
                                    USEJULIA,VOLUME,FORMULATION, ...
                                    dorsalNodesStructFromUser, ...
                                    volVector);
[uHead_c, uTail_c] = compute_head_tail_motion(eta_c, V_c, headNode, tailNode);                              

% Exact solution at D)
[eta_D, S_D, V_D] = get_results_PROM(nodes, elements,muscleBoundaries, ...
                                    myElementConstructor, nsetBC, ...                            
                                    U, xi_D, 'D)', ...
                                    kActu, tmax, h, ...
                                    USEJULIA,VOLUME,FORMULATION, ...
                                    dorsalNodesStructFromUser, ...
                                    volVector);
[uHead_D, uTail_D] = compute_head_tail_motion(eta_D, V_D, headNode, tailNode); 

% Store results
outfile = strcat(base_data_output_dir, 'CcD.mat');
outCcD = wrap_results(uHead_approx_c, uTail_approx_c, S_C, ...
                     uHead_c, uTail_c, S_c, ...
                     uHead_approx_D, uTail_approx_D, ...
                     uHead_D, uTail_D, S_D);

save(outfile, 'outCcD');

%% RUN MODELS SEGMENT D-d-E (D already ok)_________________________________

xi_d = xiAtMiddle(:,4);
xi_E = xiAtRebuild(:,5);
xi_Rebuild_d = xiRebuildAtMiddle(:,4);
xi_Rebuild_E = xiRebuildAtRebuild(:,4);

% Approximation for d)
eta_approx_d = eta_D + double(ttv(S_D,xi_Rebuild_d,2));
[uHead_approx_d, uTail_approx_d] = compute_head_tail_motion(eta_approx_d, V_D, headNode, tailNode);

% Aproximation for E)
eta_approx_E = eta_D + double(ttv(S_D,xi_Rebuild_D,2));             
[uHead_approx_E, uTail_approx_E] = compute_head_tail_motion(eta_approx_E, V_D, headNode, tailNode);

% Exact solution at d)
[eta_d, S_d, V_d] = get_results_PROM(nodes, elements,muscleBoundaries, ...
                                    myElementConstructor, nsetBC, ...                            
                                    U, xi_d, 'd)', ...
                                    kActu, tmax, h, ...
                                    USEJULIA,VOLUME,FORMULATION, ...
                                    dorsalNodesStructFromUser, ...
                                    volVector);
[uHead_d, uTail_d] = compute_head_tail_motion(eta_d, V_d, headNode, tailNode);                              

% Exact solution at E)
[eta_E, S_E, V_E] = get_results_PROM(nodes, elements,muscleBoundaries, ...
                                    myElementConstructor, nsetBC, ...                            
                                    U, xi_D, 'E)', ...
                                    kActu, tmax, h, ...
                                    USEJULIA,VOLUME,FORMULATION, ...
                                    dorsalNodesStructFromUser, ...
                                    volVector);
[uHead_E, uTail_E] = compute_head_tail_motion(eta_E, V_E, headNode, tailNode); 

% Store results
outfile = strcat(base_data_output_dir, 'DdE.mat');
outDdE = wrap_results(uHead_approx_d, uTail_approx_d, S_D, ...
                     uHead_d, uTail_d, S_d, ...
                     uHead_approx_E, uTail_approx_E, ...
                     uHead_E, uTail_E, S_E);

save(outfile, 'outDdE');

%% ANALYSIS
% Load result
savedfile = strcat(base_data_output_dir, 'AaB.mat');
load(savedfile);
savedfile = strcat(base_data_output_dir, 'BbC.mat');
load(savedfile);
savedfile = strcat(base_data_output_dir, 'CcD.mat');
load(savedfile);
savedfile = strcat(base_data_output_dir, 'DdE.mat');
load(savedfile);


%% ANALYSIS RESULTS SEGMENT A-a-B _________________________________________
timePlot = linspace(0,tmax-h,tmax/h+1);
% a)
fig = fig_comparison_PROM_approximation(timePlot, ...
                                outAaB.uHead_1, outAaB.uHead_approx_1, ...
                                outAaB.uTail_1, outAaB.uTail_approx_1, ...
                                Lx, Ly, [3 3 9 7.0]);
title('a')
% B)
fig = fig_comparison_PROM_approximation(timePlot, ...
                                outAaB.uHead_2, outAaB.uHead_approx_2, ...
                                outAaB.uTail_2, outAaB.uTail_approx_2, ...
                                Lx, Ly, [13 3 9 7.0]);
title('B')

%% ANALYSIS RESULTS SEGMENT B-b-C _________________________________________
timePlot = linspace(0,tmax-h,tmax/h+1);
% b)
fig = fig_comparison_PROM_approximation(timePlot, ...
                                outBbC.uHead_1, outBbC.uHead_approx_1, ...
                                outBbC.uTail_1, outBbC.uTail_approx_1, ...
                                Lx, Ly, [3 3 9 7.0]);
title('b')
% C)
fig = fig_comparison_PROM_approximation(timePlot, ...
                                outBbC.uHead_2, outBbC.uHead_approx_2, ...
                                outBbC.uTail_2, outBbC.uTail_approx_2, ...
                                Lx, Ly, [13 3 9 7.0]);
title('C')

%% ANALYSIS RESULTS SEGMENT C-c-D _________________________________________
timePlot = linspace(0,tmax-h,tmax/h+1);
% c)
fig = fig_comparison_PROM_approximation(timePlot, ...
                                outCcD.uHead_1, outCcD.uHead_approx_1, ...
                                outCcD.uTail_1, outCcD.uTail_approx_1, ...
                                Lx, Ly, [3 3 9 7.0]);
title('c')
% D)
fig = fig_comparison_PROM_approximation(timePlot, ...
                                outCcD.uHead_2, outCcD.uHead_approx_2, ...
                                outCcD.uTail_2, outCcD.uTail_approx_2, ...
                                Lx, Ly, [13 3 9 7.0]);
title('D')

%% ANALYSIS RESULTS SEGMENT D-d-E _________________________________________
timePlot = linspace(0,tmax-h,tmax/h+1);
% d)
fig = fig_comparison_PROM_approximation(timePlot, ...
                                outDdE.uHead_1, outDdE.uHead_approx_1, ...
                                outDdE.uTail_1, outDdE.uTail_approx_1, ...
                                Lx, Ly, [3 3 9 7.0]);
title('d')
% E)
fig = fig_comparison_PROM_approximation(timePlot, ...
                                outDdE.uHead_2, outDdE.uHead_approx_2, ...
                                outDdE.uTail_2, outDdE.uTail_approx_2, ...
                                Lx, Ly, [13 3 9 7.0]);
title('E')


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
                   
