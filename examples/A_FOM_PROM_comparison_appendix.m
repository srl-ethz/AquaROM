% ------------------------------------------------------------------------ 
% A_FOM_PROM_compariso_appendix.m
%
% Description: Accuracy and computational speed comparison between the FOM
% and the (P)ROM formulations. Additional results presented in the
% appendix of the paper.
%
% Last modified: 02/02/2026, Mathieu Dubied, ETH Zurich
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
kActu = 9.5*1e4;%0.8*1e5;
% propRigid = 0.35;
% muscleBoundaries = [0.95,0.35]
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
%
% Get maximum strain in the deformed configuration
out = Assembly.max_strain_in_structure(sol_FOM(:,stepMaxAmp));
fprintf("\nMax strain = %.3f in element %d (GP %d)\n", out.smax, out.element, out.gaussPoint);

%%
figure
elHighlightWeights = zeros(num_elements,1);
elHighlightWeights(out.element) = 1;
h = HighlightElement(nodes,elements,elHighlightWeights,"red");
%%
figure
elementPlot = elements(:,1:4);    
u = reshape(sol_FOM(:,stepMaxAmp), 3, []).';
h2 = PlotFieldonDeformedMesh_ext(nodes, elementPlot, u, 'lineWidth',0.2);

%% APPENDIX C.1 - ROB SELECTION - ITERATIVE TESTING OF MULTIPLE CASES _____

% Vector of element counts
elements_vec = [4270,8086, 16009]; % Number of elements for each input file
kActu_values = {[3.6, 3.0]*1e5 ...
                [3.1, 2.6]*1e5, ...
                [2.9, 2.5]*1e5
                };
            
% Set simulation parameters
h = 0.02;
tmax = 2.0;

% Set number of VMs and if using MDs
n_VMs = 2;      % 1 or 2
MDs_flag = 1;   % 1 (0) for (no) MDs in ROB   

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
    results_table_filename = sprintf('Results/Data/Appendix/A_results_app_%del_%dVMs_%dMDs.csv', num_elements, n_VMs, MDs_flag);
    csvwrite(results_table_filename,results_matrix);
end

%% APPENDIX C.1 - ROB SELECTION - RESULTS _________________________________
nElementsForResult = [4270, 8086, 16009]; % Number of elements for each input file

% Read files
for i=1:length(nElementsForResult)
    nElements = nElementsForResult(i);
    filename_main = sprintf('Results/Data/A_results_%del.csv', nElements);
    filename_VM1_MD0 = sprintf('Results/Data/Appendix/A_results_app_%del_%dVMs_%dMDs.csv', nElements, 1, 0);
    filename_VM2_MD1 = sprintf('Results/Data/Appendix/A_results_app_%del_%dVMs_%dMDs.csv', nElements, 2, 1);
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
actuValues = [3.6, 3.0, 3.1, 2.6, 2.9, 2.5 ]*1e5;
res_main = res_main(ismember(res_main(:,2), actuValues),[1,2,3,4,5]);
res_VM1_MD0 = res_VM1_MD0(:,[1,2,3]);
res_VM2_MD1 = res_VM2_MD1(:,[1,2,3]);

key_main = res_main(:,1:2);
[~, idx1] = ismember(key_main, res_VM1_MD0(:,1:2), 'rows');
res_VM1_MD0 = res_VM1_MD0(idx1, :);
[~, idx2] = ismember(key_main, res_VM2_MD1(:,1:2), 'rows');
res_VM2_MD1 = res_VM2_MD1(idx2, :);

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

%%  APPENDIX C.1 - ROB SELECTION - FIGURE 2ND VIBRATION MODE ______________                                                  
num_elements_fig = 8086;
filename = strcat('InputFiles/3d_rectangle_', num2str(num_elements_fig), 'el');
[Mesh_ROM, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);

% set boundary conditions
for l=1:length(nsetBC)
    Mesh_ROM.set_essential_boundary_condition([nsetBC{l}],1:3,0)  % all DOFs constrained to get VMs. Rigid body modes are added in build_ROM
end  

% FIGURE IN APPENDIX of the paper: 2nd VM
f = fig_VM(Mesh_ROM, nodes, elements,muscleBoundaries, esetBC, 2);
fig_filename = sprintf('Setup/Figures/Appendix/VM_2_%del.pdf', Mesh_ROM.nElements);
exportgraphics(f, fig_filename, 'Resolution', 1400);


