% ------------------------------------------------------------------------ 
% B_appendix_finite_difference.m
%
% Description: Compare the gradient from the PROM to FOM finite difference
% 
% Last modified: 05/02/2026, Mathieu Dubied, ETH Zurich
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
n_elements = 1272;
kActu = 10*1e4;
filename = strcat('InputFiles/3d_rectangle_', num2str(n_elements), 'el');

[MeshNominal, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
[Lx, Ly, Lz] = mesh_dimensions_3D(nodes);
n_elements = size(elements,1);

% select shape variation (one parameter at the time)
[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish, xz_concaveTail, y_fin, ...
    y_bumpBack, y_bumpFront] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);
U = [y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish, xz_concaveTail, y_fin...
    y_bumpBack, y_bumpFront];

paramIdx = 1;
U = U(:,paramIdx:2);
xi = zeros(size(U,2),1);

%% PROM

% build PROM
fprintf('____________________\n')
fprintf('Building PROM ... \n')
for l=1:length(nsetBC)
    MeshNominal.set_essential_boundary_condition([nsetBC{l}],1:3,0)   
end
[V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
build_PROM_3D(MeshNominal,nodes,elements,muscleBoundaries,U,USEJULIA,VOLUME,FORMULATION);

% Solve EoMs and get sensitivies
h = 0.02;
tmax = 2.0;
fprintf('____________________\n')
fprintf('Solving EoMs...\n') 
TI_NL_PROM = solve_EoMs_and_sensitivities(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,kActu,h,tmax);          
  
eta = TI_NL_PROM.Solution.q;
S = TI_NL_PROM.Solution.s;
eta_0k = TI_NL_PROM.Solution.q;
eta_k = eta;

% Get gradient
fprintf('____________________\n')
fprintf('Computing gradient...\n')   
wSize = 5;
nablaLr = gradient_cost_function_wo_constraints_TET4(xi,eta_k,S,V,wSize);

%% FOM
% Use Central difference
% parameter 1

% create deformed mesh 1
% solve FOM

% create deformed mesh 2
% solve FOM

% compute gradient based on the two difference













