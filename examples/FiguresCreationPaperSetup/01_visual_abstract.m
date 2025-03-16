% -------------------------------------------------------------------------
% Creation of Figure 1, Visual Abstract
%
% Last modified: 06/02/2025, Mathieu Dubied, ETH Zurich
% -------------------------------------------------------------------------
clear; 
close all; 
clc
if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end

set(groot,'defaulttextinterpreter','latex');

%% PREPARE MODEL   

% load material parameters
load(fullfile(cd, '..', 'parameters.mat'));

% specify and create FE mesh
filename = 'InputFiles/3d_rectangle_8086el'; 
[MeshNominal, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
[Lx, Ly, Lz] = mesh_dimensions_3D(nodes);
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];
elementPlot = elements(:,1:4);
n_elements = size(elements,1);
lineWidthPatchForPlot = 0.2;


%% SHAPE VARIATIONS _______________________________________________________
% gather all shape variations
[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish, xz_concaveTail] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);

% create shape variation basis U
U = [z_tail,z_head,y_linLongTail,y_head,y_ellipseFish,...
    z_smallFish, z_notch, xz_concaveTail];


%% CREATE FIGURE __________________________________________________________

% Initialize figure
f_abstract = figure('units','centimeters','position',[3 3 9.5 9]);

% Arrow from top to bottom
annotation('rectangle', [0.14, 0.3, 0.02, 0.45], 'FaceColor', 'w', 'Color', 'None', 'LineWidth', 0.5, 'FaceAlpha', 0.4);
annotation('arrow', [0.15, 0.15], [0.75, 0.3], 'LineWidth', 2);

% 1. Nominal shape
xi = zeros(8,1);
v = reshape(U*xi, 3, []).';
annotation('textbox', [0.3, 0.9, 0.6, 0.05], 'String', '1. Load nominal shape', ...
               'EdgeColor', 'none', 'FontSize', 12, 'HorizontalAlignment', 'left', 'Interpreter','latex');
axes('Position', [0.02,0.73, 0.26, 0.26]);  
hf = PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', 1, 'lineWidth', lineWidthPatchForPlot);
hold off


% 2. Shape variations
annotation('textbox', [0.3, 0.8, 0.55, 0.05], 'String', '2. Define shape variations', ...
               'EdgeColor', 'none', 'FontSize', 12,  'Interpreter','latex', ...
               'BackgroundColor','White', 'VerticalAlignment', 'middle');
ax = axes('Position', [0.3,0.57, 0.2, 0.2]);  
v = reshape(U(:,7)*0.9, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', 1, 'lineWidth', lineWidthPatchForPlot);
axes('Position', [0.425,0.57, 0.2, 0.2]);
v = reshape(U(:,3)*0.9, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', 1, 'lineWidth', lineWidthPatchForPlot);
axes('Position', [0.55,0.57, 0.2, 0.2]);
v = reshape(U(:,8)*0.5, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', 1, 'lineWidth', lineWidthPatchForPlot);
axes('Position', [0.675,0.57, 0.2, 0.2]);
v = reshape(U(:,5)*0.9, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', 1, 'lineWidth', lineWidthPatchForPlot);
annotation('textbox', [0.79,0.6, 0.6, 0.05], 'String', '...', ...
               'EdgeColor', 'none', 'FontSize', 16, 'Interpreter','latex');
           
axes('Position', [0.33,0.55, 0.62, 0.275]); 
rectangle('Position',[0,0.0,1,1])
xlim([0 1])
ylim([0 1])
axis off


% 3. Optimisation pipeline
annotation('textbox', [0.3, 0.45, 0.44, 0.05], 'String', '3. Run optimisation', ...
               'EdgeColor', 'none', 'FontSize', 12, 'Interpreter','latex', ...
               'BackgroundColor','White', 'VerticalAlignment', 'middle');
annotation('textbox', [0.35, 0.4, 0.6, 0.05], 'String', '- Create PROM', ...
               'EdgeColor', 'none', 'FontSize', 12, 'Interpreter','latex');
annotation('textbox', [0.35, 0.35, 0.6, 0.05], 'String', '- Solve EoMs and sensitivity', ...
               'EdgeColor', 'none', 'FontSize', 12, 'Interpreter','latex');
annotation('textbox', [0.35, 0.3, 0.6, 0.05], 'String', '- Take gradient step', ...
               'EdgeColor', 'none', 'FontSize', 12, 'Interpreter','latex');

% box
% annotation('textbox', [0.325, 0.22, 0.65, 0.20], 'String','');
axes('Position', [0.33,0.275, 0.62, 0.2]); 
rectangle('Position',[0,0.0,1,1])
xlim([0 1])
ylim([0 1])
axis off



% 4. Optimal shape
annotation('textbox', [0.3, 0.15, 0.6, 0.05], 'String', '4. Obtain optimal shape', ...
               'EdgeColor', 'none', 'FontSize', 12, 'Interpreter','latex');
axes('Position', [0.,0.0, 0.28, 0.28]); 
xi =  [-0.3584,0.3080,0.3893,0.2906,0.4838,0.2701,0.2271,0.297].';
v = reshape(U*xi, 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', 1, 'lineWidth', lineWidthPatchForPlot);
plotcube(L, O, .05, [0 0 0]);


% 5. In-between shape
axes('Position', [0.05,0.6, 0.15, 0.15]);
v = reshape(U*xi*0.25, 3, []).'; % 0.25 of the optimal shape
PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', 1, 'lineWidth', lineWidthPatchForPlot);

axes('Position', [0.05,0.45, 0.15, 0.15]);
v = reshape(U*xi*0.5, 3, []).'; % 0.5 of the optimal shape
PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', 1, 'lineWidth', lineWidthPatchForPlot);

axes('Position', [0.05,0.30, 0.15, 0.15]);
v = reshape(U*xi*0.75, 3, []).'; % 0.75 of the optimal shape
PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', 1, 'lineWidth', lineWidthPatchForPlot);



% Save figure
exportgraphics(f_abstract,'Figures/01_visual_abstract.pdf','Resolution',1200)
