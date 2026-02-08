% ------------------------------------------------------------------------ 
% A_appendix.m
%
% Description: Accuracy and computational speed comparison between the FOM
% and the (P)ROM formulations. Additional results presented in the
% appendix of the paper.
%
% Last modified: 08/02/2026, Mathieu Dubied, ETH Zurich
% ------------------------------------------------------------------------
clear; 
close all; 
clc
if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
load('parameters.mat')

%% EXTRA RESULT - MAX STRAIN COMPUTATION (StVK MODEL VALIDITY)_____________
num_elements = 8086;
kActu = 7.0*1e4;

filename = strcat('InputFiles/3d_rectangle_', num2str(num_elements), 'el');  % Construct filena
[Mesh_FOM, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
% muscleBoundaries = [0.9, 0.25]
[Lx, Ly, Lz] = mesh_dimensions_3D(nodes);
    
% Set boundary conditions for FOM and ROM
for l = 1:length(nsetBC)
    Mesh_FOM.set_essential_boundary_condition([nsetBC{l}], 2:3, 0);
end

% Simulate FOM for half an oscillation period
h = 0.02;
tmax = 0.3;
tStartFOM = tic;
fprintf('\nBuilding and solving FOM for %d elements, kActu: %.3f...\n', num_elements, kActu);
[Assembly, tailProperties, spineProperties, dragProperties, actuLeft, actuRight] = ...
    build_FOM_3D(Mesh_FOM, nodes, elements, muscleBoundaries);
TI_NL_FOM = solve_EoMs_FOM(Assembly, elements, tailProperties, spineProperties, dragProperties, actuLeft, actuRight, kActu, h, tmax);
timeFOM = toc(tStartFOM);
sol_FOM = Assembly.unconstrain_vector(TI_NL_FOM.Solution.q);
fprintf("\nBuilding and solving FOM in %.2f\n",timeFOM);

% Obtain moment with maximal amplitude and store the deformation
idx = tailProperties.tailNode*3-2 : tailProperties.tailNode*3;
uTail_FOM = sol_FOM(idx, 1:tmax/h) * 100;
[maxAmplitude, stepMaxAmp] = max(uTail_FOM(2,:));
fprintf("\nMaximum amplitude %.2f reached at %.2f s\n", maxAmplitude, stepMaxAmp*h);

% Get maximum strain in the deformed configuration
out = Assembly.max_strain_in_structure(sol_FOM(:,stepMaxAmp));
fprintf("\nMax strain = %.3f in element %d (GP %d)\n", out.smax, out.element, out.gaussPoint);

% Figure showing the critical element
figure
elHighlightWeights = zeros(num_elements,1);
elHighlightWeights(out.element) = 1;
h = HighlightElement(nodes,elements,elHighlightWeights,"red");

% Figure showing the deformation state at max oscillation
figure
elementPlot = elements(:,1:4);    
u = reshape(sol_FOM(:,stepMaxAmp), 3, []).';
h2 = PlotFieldonDeformedMesh_ext(nodes, elementPlot, u, 'lineWidth',0.2);

%% APPENDIX C.1 - ROB SELECTION - ITERATIVE TESTING OF MULTIPLE CASES _____

% Vector of element counts
elements_vec = [4270,8086, 16009]; % Number of elements for each input file
kActu_values = {[8.5, 7]*1e4 ...
                [8.5, 7]*1e4, ...
                [6.6, 5.5]*1e4
                };
            
% Set simulation parameters
h = 0.02;
tmax = 2.0;

% Set number of VMs and if using MDs
n_VMs = 1;      % 1 or 2
MDs_flag = 0;   % 1 (0) for (no) MDs in ROB   

% Loop over each element count and each actuation value
for elem_idx = 1:length(elements_vec)
    num_elements = elements_vec(elem_idx); % Get current number of elements
    kValues_for_n_el = kActu_values{elem_idx};
    % Pre-allocate matrix to store results. 7 columns
    % Columns:  [num_elements, kActu, ...
    %           max_uHead_ROM, max_uTail_ROM, ...
    %           timeROMBuild, timeROMSolve, timeROM];
    results_matrix = zeros(length(kValues_for_n_el), 7); 
    
    % Load the FE mesh for FOM and ROM
    filename = strcat('InputFiles/3d_rectangle_', num2str(num_elements), 'el');  % Construct filena
    [Mesh_ROM, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
    [Lx, Ly, Lz] = mesh_dimensions_3D(nodes);

    % Set boundary conditions for FOM and ROM
    for l = 1:length(nsetBC)
        Mesh_ROM.set_essential_boundary_condition([nsetBC{l}], 1:3, 0);
    end
    
    for k_idx = 1:length(kValues_for_n_el)
        kActu = kValues_for_n_el(k_idx);  % Get current actuation value
    
        % ROM Simulation __________________________________
        tStartROM = tic;
        fprintf('\nBuilding and solving ROM for %d elements, kActu: %.3f...\n', num_elements, kActu);
        [V, ROM_Assembly, tensors_ROM, tailProperties, spineProperties, dragProperties, actuLeft, actuRight] = ...
            build_ROM_3D(Mesh_ROM, nodes, elements, muscleBoundaries, USEJULIA, n_VMs, MDs_flag);
        timeROMBuild = toc(tStartROM);
        tStartROMSolve = tic;
        TI_NL_ROM = solve_EoMs(V, ROM_Assembly, tensors_ROM, tailProperties, spineProperties, dragProperties, actuLeft, actuRight, kActu, h, tmax);
        timeROMSolve = toc(tStartROMSolve);
        timeROM = toc(tStartROM);       
        
        % Data Analysis ___________________________________
        fprintf('\nPost-processing results for %d elements, kActu: %.3f...\n\n\n', num_elements, kActu)
        timePlot = linspace(0, tmax-h, tmax/h);
        x0Tail = min(nodes(:, 1));
        
        headNode = find_node(0, 0, 0, nodes);
        
        uHead_ROM = zeros(3, tmax/h);
        uTail_ROM = zeros(3, tmax/h);
        
        for t = 1:tmax/h                   
            % ROM displacement, convert to cm
            uHead_ROM(:, t) = V(headNode*3-2:headNode*3, :) * TI_NL_ROM.Solution.q(:, t)*100;
            uTail_ROM(:, t) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3, :) * TI_NL_ROM.Solution.q(:, t)*100;  
        end
        
        % Calculate maximum displacements (lateral, 2, and horizontal, 1)  
        max_uHead_ROM = max(uHead_ROM(1, :));
        max_uTail_ROM = max(uTail_ROM(2, :));
     
        results_matrix(k_idx, :) = [num_elements, kActu, ...
              max_uHead_ROM, max_uTail_ROM, ... 
              timeROMBuild, timeROMSolve, timeROM];        
    end
    
    % Save summary table for each mesh
    results_table_filename = sprintf('Results/Data/Appendix/ROB_selection/A_results_app_%del_%dVMs_%dMDs.csv', num_elements, n_VMs, MDs_flag);
    csvwrite(results_table_filename,results_matrix);
end

%% APPENDIX - ROB SELECTION - RESULTS _____________________________________
nElementsForResult = [4270, 8086, 16009]; % Number of elements for each input file

% Read files
for i=1:length(nElementsForResult)
    nElements = nElementsForResult(i);
    filename_main = sprintf('Results/Data/A_FOM_PROM/A_results_%del.csv', nElements);
    filename_VM1_MD0 = sprintf('Results/Data/Appendix/ROB_selection/A_results_app_%del_%dVMs_%dMDs.csv', nElements, 1, 0);
    filename_VM2_MD1 = sprintf('Results/Data/Appendix/ROB_selection/A_results_app_%del_%dVMs_%dMDs.csv', nElements, 2, 1);
    if i==1
        res_main = readmatrix(filename_main);
        res_VM1_MD0 = readmatrix(filename_VM1_MD0);
        res_VM2_MD1 = readmatrix(filename_VM2_MD1);
    else
        res_main = [res_main;readmatrix(filename_main)];
        res_VM1_MD0 = [res_VM1_MD0;readmatrix(filename_VM1_MD0)];
        res_VM2_MD1 = [res_VM2_MD1;readmatrix(filename_VM2_MD1)];
    end   
end

% Select and reorder rows
% main: nElements, actu, max_uHead_FOM, max_uTail_FOM, max_uHead_ROM
% VMX_MDY: nElements, actu, max_uHead_ROM,
validPairs = [  4270 8.5e4
                4270 7.0e4
                8086 8.5e4
                8086 7.0e4
                16009 6.6e4
                16009 5.5e4];

mask = ismember(res_main(:,1:2), validPairs, 'rows');
res_main = res_main(mask, [1,2,3,4,5]);
res_VM1_MD0 = res_VM1_MD0(:,[1,2,3]);
res_VM2_MD1 = res_VM2_MD1(:,[1,2,3]);

key_main = res_main(:,1:2);
[~, idx1] = ismember(key_main, res_VM1_MD0(:,1:2), 'rows');
mask = idx1 > 0;
res_VM1_MD0 = res_VM1_MD0(idx1(mask), :);
[~, idx2] = ismember(key_main, res_VM2_MD1(:,1:2), 'rows');
mask = idx2 > 0;
res_VM2_MD1 = res_VM2_MD1(idx2(mask), :);

% Compute differences
res_main = [res_main,(res_main(:, 5)-res_main(:, 3))./(res_main(:,3))*100] 
res_VM1_MD0 = [res_VM1_MD0,(res_VM1_MD0(:, 3)-res_main(:, 3))./(res_main(:,3))*100]
res_VM2_MD1 = [res_VM2_MD1,(res_VM2_MD1(:, 3)-res_main(:, 3))./(res_main(:,3))*100]

% Create result table
tab_VM1_MD0 = [res_main(:,[1,4,6]), res_VM1_MD0(:,end)];
colNames = {'nEl','y-oscillation','rel. error main [%]',' rel. error VM1 MD0 [%]'};
T = array2table(tab_VM1_MD0, 'VariableNames', colNames)

tab_VM2_MD1 = [res_main(:,[1,4,6]), res_VM2_MD1(:,end)];
colNames = {'nEl','y-oscillation','rel. error main [%]',' rel. error VM2 MD1 [%]'};
T = array2table(tab_VM2_MD1, 'VariableNames', colNames)

%%  APPENDIX - ROB SELECTION - FIGURE 2ND VIBRATION MODE __________________                                                  
num_elements_fig = 8086;
filename = strcat('InputFiles/3d_rectangle_', num2str(num_elements_fig), 'el');
[Mesh_ROM, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);

% set boundary conditions
for l=1:length(nsetBC)
    Mesh_ROM.set_essential_boundary_condition([nsetBC{l}],1:3,0)  % all DOFs constrained to get VMs. Rigid body modes are added in build_ROM
end  

% FIGURE IN APPENDIX of the paper: 2nd VM
f = fig_VM(Mesh_ROM, nodes, elements,muscleBoundaries, esetBC, 2);
fig_filename = sprintf('Setup/Figures/Appendix/FOMVM_2_%del.pdf', Mesh_ROM.nElements);
exportgraphics(f, fig_filename, 'Resolution', 1400);


%% APPENDIX - CONVERGENCE ANALYSIS OF PROM ________________________________
% Set simulation parameters
hVec = [0.02, 0.01, 0.005];
tmax = 2.0;
kActu = 7.0*1e4;
num_elements = 8086;

% Loop over each element count and each actuation value
for h_idx = 1:length(hVec)
    
    % Time step
    h = hVec(h_idx);
    
    fprintf('__________________________\n')
    fprintf('Analysis for h = %.2f s\n', h)
    
    % Load the FE mesh for FOM and ROM
    filename = strcat('InputFiles/3d_rectangle_', num2str(num_elements), 'el');  % Construct filena
    [Mesh_ROM, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
    [Lx, Ly, Lz] = mesh_dimensions_3D(nodes);
    
    % shape variations for PROM
    [y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
        y_tail,y_head,y_linLongTail,y_ellipseFish, xz_concaveTail] = ...
            shape_variations_3D(nodes,Lx,Ly,Lz);
    
    U = [z_tail,z_head,y_thinFish]; % 3 parameters
    U = [z_tail,z_head,y_thinFish]; % 3 parameters

    % Set boundary conditions
    for l = 1:length(nsetBC)
        Mesh_ROM.set_essential_boundary_condition([nsetBC{l}], 1:3, 0);
    end
    
    % PROM Simulation _________________________________
    tStartPROM = tic;
    volVector = compute_nominal_vol_per_element(Mesh_ROM,size(elements,1));
    fprintf('\nBuilding and solving PROM (3 param., w/ sens.)...\n')      
    [V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuLeft,actuRight] = ...
        build_PROM_3D(Mesh_ROM, nodes, elements, muscleBoundaries, U, USEJULIA, VOLUME, FORMULATION, 'nomVolVector', volVector); 

    timePROMBuild = toc(tStartPROM);
    tStartPROMSolve = tic;
    TI_NL_PROM = solve_EoMs_and_sensitivities(V, PROM_Assembly, tensors_PROM, tailProperties, spineProperties, dragProperties, actuLeft, actuRight, kActu, h, tmax); 
    timePROMSolve = toc(tStartPROMSolve);
    timePROM = toc(tStartPROM);
    fprintf('Time for PROM 3p: %.2f min\n', timePROM/60)
    
    % Data Analysis ___________________________________
    fprintf('\nPost-processing results for time step %.3f ...\n\n', h)
    timePlot = linspace(0, tmax-h, tmax/h);
    x0Tail = min(nodes(:, 1));
    headNode = find_node(0, 0, 0, nodes);
   
    uHead = zeros(3, int8(tmax/h));
    uTail = zeros(3, int8(tmax/h));
    
    for t = 1:tmax/h
        % PROM (3 params) displacement, convert to cm
        uHead(:, t) = V(headNode*3-2:headNode*3, :) * TI_NL_PROM.Solution.q(:, t)*100;
        uTail(:, t) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3, :) * TI_NL_PROM.Solution.q(:, t)*100;  
    end
    
    % Calculate maximum displacements (lateral, 2, and horizontal, 1)
    max_uHead_PROM_3 = max(uHead(1, :));
    max_uTail_PROM_3 = max(uTail(2, :)); 
    
    % Output struct
    out = struct();
    out.h = h;
    out.uHead = uHead;
    out.uTail = uTail;
    out.uHead_max = max_uHead_PROM_3;
    out.uTail_max = max_uTail_PROM_3;
    out.S = TI_NL_PROM.Solution.s;
    out.soltime = TI_NL_PROM.Solution.soltime;
    out.itVec = TI_NL_PROM.Solution.itVec;
    out.epsilonVec = TI_NL_PROM.Solution.epsilonVec;
    
    % Save results
    results_filename = sprintf('Results/Data/Appendix/Convergence/A_convergence_h_%.3f.mat', h);
    save(results_filename, 'out');
end

%% APPENDIX - CONVERGENCE ANALYSIS OF PROM - RESULTS ______________________

% Load results
filename = sprintf('Results/Data/Appendix/Convergence/A_convergence_h_0.020.mat');
load(filename);
out1 = out;
filename = sprintf('Results/Data/Appendix/Convergence/A_convergence_h_0.010.mat');
load(filename);
out2 = out;
filename = sprintf('Results/Data/Appendix/Convergence/A_convergence_h_0.005.mat');
load(filename);
out3 = out;

%%
% 1) Plot swimming dynamics
uHead_h1 = out1.uHead(1,:); 
uHead_h2 = out2.uHead(1,:);
uHead_h3 = out3.uHead(1,:);
uTail_h1 = out1.uTail(2,:); 
uTail_h2 = out2.uTail(2,:);
uTail_h3 = out3.uTail(2,:);
n1 = numel(u1);
t1 = (0:n1-1) * out1.h;
u2 = out2.uHead(1,:); 
n2 = numel(u2);
t2 = (0:n2-1) * out2.h;
u3 = out3.uHead(1,:); 
n3 = numel(u3);
t3 = (0:n3-1) * out3.h;

fig = fig_trajectory_time_step_refinement(t1,t2,t3, ...
                                uHead_h1, uHead_h2, uHead_h3, ...
                                uTail_h1, uTail_h2, uTail_h3, ...
                                Lx, Ly, [3 3 9 7.0]);
fig_filename = sprintf('Results/Figures/Appendix/Convergence/A_convergence_trajectories.pdf');
exportgraphics(fig, fig_filename, 'Resolution', 600);                           

% 2) Focus on head position: successive errors and refinement ratio
E_12 = norm(uHead_h1 - uHead_h2(1:2:end))
E_23 = norm(uHead_h2 - uHead_h3(1:2:end))
ref_ratio = E_12/E_23


% 3) Sensitivities
S1 = out1.S;
S2 = out2.S;
S3 = out3.S;

% Compute Frobenius norm at each time step
E_12_vec = [];
E_23_vec = [];

for i=1:length(S1(1,1,:))
    E_12_vec = [E_12_vec, norm(S1(:,:,i) - S2(:,:,2*i-1), 'fro')];
end

for i=1:length(S2(1,1,:))
    E_23_vec = [E_23_vec, norm(S2(:,:,i) - S3(:,:,2*i-1), 'fro')];
end

% Compute refinement ratio 
E_12 = norm(E_12_vec)
E_23 = norm(E_23_vec)
ref_ratio = E_12/E_23

