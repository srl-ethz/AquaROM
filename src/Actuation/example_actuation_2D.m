% ------------------------------------------------------------------------ 
% example_actuation_2D.m
%
% Description: Sanity check for the muscle force formulation for 2D
% examples
%
% Last modified: 17/12/2024, Mathieu Dubied, ETH Zurich
% ------------------------------------------------------------------------
clear; 
close all; 
clc
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');

whichModel = 'CUSTOM'; % or 'ABAQUS' CUSTOM
elementType = 'QUAD4'; % or 'TRI3' QUAD4
filename = 'InputFiles/beam2D_272el'; % or 'beam2D_272el' 144

%% PREPARE MODEL __________________________________________________________     

% MATERIAL ________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 
thickness = .1;     % [m] beam's out-of-plane thickness

Material = KirchoffMaterial();
set(Material,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
Material.PLANE_STRESS = true;	% set "false" for plane_strain

% MESH_____________________________________________________________________
switch elementType
    case 'QUAD4'
        myElementConstructor = @()Quad4Element(thickness, Material);
    case 'TRI3'
        myElementConstructor = @()Tri3Element(thickness, Material);
end
Lx = 0.2;
Ly = 0.1;
nx = 16;
ny = 8;

% create mesh
switch upper( whichModel )
    case 'CUSTOM'
        [nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny,elementType);
    case 'ABAQUS'
        [nodes, elements, ~, ~] = mesh_ABAQUSread(filename);
end

Mesh = Mesh(nodes);
Mesh.create_elements_table(elements,myElementConstructor);

% boundary conditions (nset{3} contains the nodes of the right edge)
nel = size(elements,1);
Mesh.set_essential_boundary_condition(nset{3},1:2,0) 

% create assembly
Assembly = Assembly(Mesh);

% plot mesh
figure('units','centimeters','position',[1 11 16 7],'Name','Mesh')
PlotMesh(nodes, elements, 0);   
axis on

%% DEFINE MUSCLE ELEMENTS _________________________________________________

% Ok for QUAD4 elements
actuationDirection = [1;0;0]; %[1;0]-->[1;0;0] (Voigt notation)

% top muscle
topMuscle = zeros(nel,1);
for el=1:nel
    switch elementType
        case 'QUAD4'
            [elCenterX,elCenterY] = elementCenter_QUAD4(nodes,elements,el);
        case 'TRI3'
            [elCenterX,elCenterY] = elementCenter_TRI3(nodes,elements,el);
    end
    if elCenterY > Ly/2 &&  elCenterX > 0.5*Lx && elCenterX > 0
        topMuscle(el) = 1;
    end    
end

actuTop = compute_actuation_tensors_FOM(Assembly,topMuscle,actuationDirection);

% bottom muscle
bottomMuscle = zeros(nel,1);
for el=1:nel
    switch elementType
        case 'QUAD4'
            [elCenterX,elCenterY] = elementCenter_QUAD4(nodes,elements,el);
        case 'TRI3'
            [elCenterX,elCenterY] = elementCenter_TRI3(nodes,elements,el);
    end
    if elCenterY < Ly/2 &&  elCenterX > 0.5*Lx && elCenterX > 0
        bottomMuscle(el) = 1;
    end    
end
actuBottom = compute_actuation_tensors_FOM(Assembly, bottomMuscle, actuationDirection);

% plot muscles
figure('units','centimeters','position',[1 2 16 7],'Name','Muscles')
PlotMeshWith2Sets(nodes, elements, topMuscle, 'r', bottomMuscle, 'b', 'show',0);   
axis off

%% PLOT FORCES ON MESH ____________________________________________________
% define force
B1T = actuTop.B1;
B1B = actuBottom.B1;
B2T = actuTop.B2;
B2B = actuBottom.B2; 
kActu = 10;
a = 1;
fActu = @(q)  kActu/2*(a*(B1T+B2T*q) + (-a)*(B1B+B2B*q));

% UNDEFORMED MESH _________________________________________________________
% actuation force for u0
u0 = zeros(Mesh.nDOFs, 1);
fActu_0 = fActu(u0);
fActu_0 = reshape(fActu_0, 2, []).';

% plot
figure('units','centimeters','position',[17 11 16 7],'Name','Actuation forces (undeformed mesh)')
PlotMeshWith2Sets(nodes, elements, topMuscle, 'r', bottomMuscle, 'b', 'show',0);   
quiver(nodes(:,1),nodes(:,2),fActu_0(:,1),fActu_0(:,2),'k','linewidth',1)

%%
% DEFORMED MESH ___________________________________________________________
% get first vibration mode as a deformed state
Mn = Assembly.mass_matrix();
[Kn,~] = Assembly.tangent_stiffness_and_force(u0);
n_VMs = 1;
Kc = Assembly.constrain_matrix(Kn);
Mc = Assembly.constrain_matrix(Mn);
[VM,~] = eigs(Kc, Mc, n_VMs, 'SM');
VM = Assembly.unconstrain_vector(VM);
u1 = 0.05*VM;
nodesDeformed = nodes + reshape(u1, 2, []).';

% get corresponding actuation force for u1
fActu_1 = fActu(u1);
fActu_1 = reshape(fActu_1, 2, []).';

% plot
figure('units','centimeters','position',[17 2 16 7],'Name','Actuation forces (deformed mesh)')
PlotMeshWith2Sets(nodesDeformed, elements, topMuscle, 'r', bottomMuscle, 'b', 'show',0);
quiver(nodesDeformed(:,1),nodesDeformed(:,2),fActu_1(:,1),fActu_1(:,2),'k','linewidth',1)


%% SIMULATION _____________________________________________________________


%% HELPER FUNCTIONS _______________________________________________________
function [elCenterX,elCenterY] = elementCenter_QUAD4(nodes,elements,idx)
    elCenterX = (nodes(elements(idx,1),1)+nodes(elements(idx,2),1)+nodes(elements(idx,3),1)+nodes(elements(idx,4),1))/4;
    elCenterY = (nodes(elements(idx,1),2)+nodes(elements(idx,2),2)+nodes(elements(idx,3),2)+nodes(elements(idx,4),2))/4; 
end

function [elCenterX,elCenterY] = elementCenter_TRI3(nodes,elements,idx)
    elCenterX = (nodes(elements(idx,1),1)+nodes(elements(idx,2),1)+nodes(elements(idx,3),1))/3;
    elCenterY = (nodes(elements(idx,1),2)+nodes(elements(idx,2),2)+nodes(elements(idx,3),2))/3; 
end

