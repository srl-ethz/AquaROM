% ------------------------------------------------------------------------ 
% A_appendix_simple_simulations.m
%
% Description: Forward simulate some of the models
% 
% Last modified: 10/02/2026, Mathieu Dubied, ETH Zurich
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
tmax = 10.0;
expName = "SO1";

% load mesh and optimization results
filename = strcat('InputFiles/3d_rectangle_', num2str(n_elements), 'el');
results_filename = sprintf('Results/Data/B_shape_optimization/%s_8086_el_kActu_%.3f.mat',expName, kActu);
load(results_filename)

%% PREPARE NOMINAL MODEL FOR REPRODUCIBILITY ______________________________
% nominal model
[MeshNominal, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
[Lx, Ly, Lz] = mesh_dimensions_3D(nodes);

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


%% RUN MODELS: NOMINAL ____________________________________________________

xi = [0;0;0];
U = [z_tail,z_head,y_thinFish]; 

[eta_nom, S_nom, V_nom]= get_results_PROM(nodes,elements,muscleBoundaries, ...
                                        myElementConstructor, nsetBC, ...                            
                                        U, xi, 'nominal', ...
                                        kActu, tmax, h, ...
                                        USEJULIA,VOLUME,FORMULATION, ...
                                        dorsalNodesStructFromUser, ...
                                        volVector) ;
                                    
[uHead_nom, uTail_nom] = compute_head_tail_motion(eta_nom, V_nom, headNode, tailNode); 

%%
timePlot = linspace(0,tmax-h,tmax/h+2);
f = fig_head_tail_motion_simple(timePlot, uHead_nom, uTail_nom, Lx, Ly, [3 3 18 5.0]);

     

%% RUN MODELS: SO1 ________________________________________________________
                                    %%_
% PROM used in optimization: xi=xi_0
[eta_A, S_A, V_A] = get_results_PROM(nodes, elements,muscleBoundaries, ...
                                    myElementConstructor, nsetBC, ...                            
                                    U, xi_A, 'A)', ...
                                    kActu, tmax, h, ...
                                    USEJULIA,VOLUME,FORMULATION, ...
                                    dorsalNodesStructFromUser, ...
                                    volVector);
[uHead_A, uTail_A] = compute_head_tail_motion(eta_A, V_A, headNode, tailNode); 

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

%%
%% ANIMATION ______________________________________________________________
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

nStep = size(eta_nom,2);
for t = 1:nStep              
    u(:, t) = V_nom*eta_nom(:,t);
end

actuationValues = zeros(size(u,2),1);
for t=1:size(u,2)
    actuationValues(t) = 0;
end

actuationValues2 = zeros(size(u,2),1);
for t=1:size(u,2)
    actuationValues2(t) = 0;
end
sol = u(:,1:end);
AnimateFieldonDeformedMeshActuation2Muscles(nodes, elementPlot,topMuscle,actuationValues,...
    bottomMuscle,actuationValues2,sol, ...
    'factor',1,'index',1:3,'filename','result_video','framerate',1/h)



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
function out = wrap_results(uHead, uTail, eta, V)
    out = struct();
    
    out.uHead = uHead;
    out.uTail = uTail;  
    out.eta   = eta;
    out.V     = V;

end

% compute grad etrelative error per element of the gradient
function [grad, grad_approx, rel_err] = compute_gradient_rel_error(xi, out, idx)
    if idx == 1
        grad_approx = gradient_cost_function_wo_constraints_TET4(xi, out.eta_approx_1,out.S_0,out.V0,5);
        grad= gradient_cost_function_wo_constraints_TET4(xi,out.eta_1,out.S_1,out.V1,5);
        rel_err = (grad_approx - grad) ./ abs(grad);
    elseif idx == 2
        grad_approx = gradient_cost_function_wo_constraints_TET4(xi, out.eta_approx_2,out.S_0,out.V0,5);
        grad = gradient_cost_function_wo_constraints_TET4(xi, out.eta_2,out.S_2,out.V2,5);
        rel_err = (grad_approx - grad) ./ abs(grad);
    end
    
end


                   
