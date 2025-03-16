% -------------------------------------------------------------------------
% Vibration modes analysis
%
% Last updated: 03/03/2024, Mathieu Dubied, ETH Zurich
% -------------------------------------------------------------------------
clear; 
close all; 
clc

elementType = 'TET4';
fish = 0;

showMesh = 0;

%% PREPARE MODEL __________________________________________________________                                                   

% DATA ____________________________________________________________________
E       = 260000;      % Young's modulus [Pa]
rho     = 1070;         % density [kg/m^3]
nu      = 0.4;        % Poisson's ratio 

% material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	    % set "false" for plane_strain
switch elementType
    case 'TET4'
        myElementConstructor = @()Tet4Element(myMaterial);
        filename = '3d_rectangle_660el';
    case 'TET10'
        myElementConstructor = @()Tet10Element(myMaterial);
        filename = '3d_rectangle_660el_quad_tet';
end

% MESH ____________________________________________________________________

if fish
    filename = '3d_fish_for_mike';
end

[nodes, elements, ~, elset] = mesh_ABAQUSread(filename);
nel = size(elements,1);

if ~fish
    nodes = 0.01*nodes ;    % convert from cm to m
end

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of airfoil
Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of airfoil
Lz = abs(max(nodes(:,3))-min(nodes(:,3)));  % vertical length of airfoil

% plot nominal mesh
if showMesh
    elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
    figure('units','normalized','position',[.2 .1 .6 .8])
    PlotMeshAxis(nodes, elementPlot, 0);
    hold off
end

% boundary conditions of nominal mesh - Set 1
nset = {};
for el=1:nel   
    for n=1:size(elements,2)
        if  ~any(cat(2, nset{:}) == elements(el,n))
            nset{end+1} = elements(el,n);
        end
    end   
end

for l=1:length(nset)
    MeshNominal.set_essential_boundary_condition([nset{l}],3,0)
end   

% boundary condition of nominal mesh - Set 2
% desiredFixedPoint = - 0.45*Lx;
% fixedPoint = find_fixed_point(nodes,desiredFixedPoint);
% nel = size(elements,1);
% nset = {};
% marginFixedPoint = 0.02;
% 
% for el=1:nel   
%     for n=1:size(elements,2)
%         if  nodes(elements(el,n),1)> fixedPoint-marginFixedPoint && ...
%                 nodes(elements(el,n),1)< fixedPoint+marginFixedPoint && ...
%                 ~any(cat(2, nset{:}) == elements(el,n))
%             nset{end+1} = elements(el,n);
%         end
%     end   
% end
% 
% for l=1:length(nset)
%     MeshNominal.set_essential_boundary_condition([nset{l}],1:2,0)
% %     MeshNominal.set_essential_boundary_condition([nset{l}],1:3,0)
% end 

% ASSEMBLY _______________________________________________________________
NominalAssembly = Assembly(MeshNominal);
Mn = NominalAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( MeshNominal.nDOFs, 1);
[Kn,~] = NominalAssembly.tangent_stiffness_and_force(u0);

%% OBTAIN VIBRATION (AND RIGID BODY) MODES ________________________________

% solve eigenvalue problem
nVMs = 8;
Kc = NominalAssembly.constrain_matrix(Kn);
Mc = NominalAssembly.constrain_matrix(Mn);
[VMn,om] = eigs(Kc, Mc, nVMs, 'smallestabs','tol',1e-24);
[f0n,ind] = sort(sqrt(diag(om))/2/pi);
VMn = VMn(:,ind);
for ii = 1:nVMs
    VMn(:,ii) = VMn(:,ii)/max(sqrt(sum(VMn(:,ii).^2,2)));
end
VMn = NominalAssembly.unconstrain_vector(VMn);
if ~isreal(VMn)
    fprintf('Modes contain non-real parts \n')
end

% VMn  = orth(VMn)*10;

% plot the vibration modes
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
f1 = figure('units','centimeters','position',[1 2 28 15]);
nRows = int16(nVMs/4);
t = tiledlayout(nRows,4);

elementPlot = elements(:,1:4); 
for mod = 1:nVMs
    nexttile
    v = reshape(VMn(:,mod), 3, []).';
    PlotMesh(nodes, elementPlot, 0);
    PlotFieldonDeformedMesh(nodes, elementPlot, real(v), 'factor', max(nodes(:,2)));

    currentTitle = ['$$\Phi_{',num2str(mod),'}, f=',num2str(real(f0n(mod)),'%.1f'),'$$ Hz'];
    title(currentTitle,'Interpreter','latex')
end
figTitle = [elementType, ', $$\nu=',num2str(nu,'%.3f '),'$$'];
title(t,figTitle,' ','Interpreter','latex')

if fish
    start = 'Fish_';
else
    start = 'Rect_';
end

documentTitle = [start,elementType, '_nu_',num2str(nu,'%.3f '),'_constrained','.svg'];
print(f1,documentTitle,'-dsvg','-r800');























