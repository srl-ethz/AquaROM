% -------------------------------------------------------------------------
% Creation of a figure for the different swimming distances
%
% Last modified: 28/01/2024, Mathieu Dubied, ETH Zurich
% -------------------------------------------------------------------------
clear; 
close all; 
clc

elementType = 'TET4';

FORMULATION = 'N1t'; % N1/N1t/N0
VOLUME = 1;         % integration over defected (1) or nominal volume (0)

USEJULIA = 1;

%% PREPARE MODEL                                                    
% DATA ____________________________________________________________________
E       = 2600000;      % Young's modulus [Pa]
rho     = 1070;         % density [kg/m^3]
nu      = 0.499;        % Poisson's ratio 

% material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	    % set "false" for plane_strain
myElementConstructor = @()Tet4Element(myMaterial);

% MESH ____________________________________________________________________

% nominal mesh
filename = '3d_rectangle_660el';%'fish3_664el';
[nodes, elements, ~, elset] = mesh_ABAQUSread(filename);

nodes = nodes*0.01;
nodes(:,2)=0.8*nodes(:,2);

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of airfoil
Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of airfoil
Lz = abs(max(nodes(:,3))-min(nodes(:,3)));  % vertical length of airfoil

% plot nominal mesh
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMeshAxis(nodes, elementPlot, 0);
hold off

% boundary conditions of nominal mesh
nel = size(elements,1);
nset = {};
for el=1:nel   
    for n=1:size(elements,2)
        if  nodes(elements(el,n),1)>-Lx*0.1 && ~any(cat(2, nset{:}) == elements(el,n))
            nset{end+1} = elements(el,n);
        end
    end   
end

%% SHAPE VARIATIONS _______________________________________________________
% gather all shape variations
[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);

% create shape variation basis U
U = [z_tail,y_head];
U = [z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_head,y_linLongTail,y_ellipseFish];

% set xi
% xi = [0.16;-0.4831;0.3;-0.2801;0.2;0.2;0.475;0.45];
xi = zeros(8,1);

% plot
f0 = figure('units','centimeters','position',[3 3 10 7],'name','Shape-varied mesh');
elementPlot = elements(:,1:4); hold on 
v1 = reshape(U*xi, 3, []).';
S = 1;
hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];
plotcube(L,O,.05,[0 0 0]);
axis equal; grid on; box on; 

%% NOMINAL: BUILD ROM AND SOLVE EOMS ______________________________________
% Parameters
h = 0.005;
tmax = 2.0;
xi = zeros(8,1);
pActu = [0.2;2.0;0.0];

% Mesh        
df = U*xi;                    % displacement fields introduced by defects
ddf = [df(1:3:end) df(2:3:end) df(3:3:end)]; 
nodes_sv = nodes + ddf;   % nominal + d ---> defected 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);
for l=1:length(nset)
    svMesh.set_essential_boundary_condition([nset{l}],1:3,0)   
end

% build ROM
fprintf('____________________\n')
fprintf('Building ROM ... \n')

[V,ROM_Assembly,tensors_ROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
    build_ROM_3D(svMesh,nodes_sv,elements,USEJULIA);      

% Solve EoMs
tic 
fprintf('____________________\n')
fprintf('Solving EoMs...\n') 
TI_NL_ROM = solve_EoMs_and_sensitivities_actu(V,ROM_Assembly, ...
    tensors_ROM,tailProperties,spineProperties,dragProperties, ...
    actuTop,actuBottom,h,tmax,pActu);
toc

% Retrieving solutions    
solNom = TI_NL_ROM.Solution.u;

%% ACTU ONLY: BUILD ROM AND SOLVE EOMS ____________________________________
% Parameters
h = 0.005;
tmax = 2.0;
xi = zeros(8,1);
pActu = [0.3027,2.58,1.0];

% Mesh        
df = U*xi;                    % displacement fields introduced by defects
ddf = [df(1:3:end) df(2:3:end) df(3:3:end)]; 
nodes_sv = nodes + ddf;   % nominal + d ---> defected 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);
for l=1:length(nset)
    svMesh.set_essential_boundary_condition([nset{l}],1:3,0)   
end

% build ROM
fprintf('____________________\n')
fprintf('Building ROM ... \n')

[V,ROM_Assembly,tensors_ROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
    build_ROM_3D(svMesh,nodes_sv,elements,USEJULIA);      

% Solve EoMs
tic 
fprintf('____________________\n')
fprintf('Solving EoMs...\n') 
TI_NL_ROM = solve_EoMs_and_sensitivities_actu(V,ROM_Assembly, ...
    tensors_ROM,tailProperties,spineProperties,dragProperties, ...
    actuTop,actuBottom,h,tmax,pActu);
toc

% Retrieving solutions    
solActu = TI_NL_ROM.Solution.u;


%% SHAPE ONLY: BUILD ROM AND SOLVE EOMS ___________________________________
% Parameters
h = 0.005;
tmax = 2.0;
xi = [0.1973;-0.4848;0.2889;-0.2650;0.1222;0.3886;0.3969;0.3889];
pActu = [0.2;2.0;0];

% Mesh        
df = U*xi;                    % displacement fields introduced by defects
ddf = [df(1:3:end) df(2:3:end) df(3:3:end)]; 
nodes_sv = nodes + ddf;   % nominal + d ---> defected 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);
for l=1:length(nset)
    svMesh.set_essential_boundary_condition([nset{l}],1:3,0)   
end

% build ROM
fprintf('____________________\n')
fprintf('Building ROM ... \n')

[V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
    build_PROM_3D(svMesh,nodes_sv,elements,U,USEJULIA,VOLUME,FORMULATION);     

% Solve EoMs
tic 
fprintf('____________________\n')
fprintf('Solving EoMs...\n') 
TI_NL_PROM = solve_EoMs_and_sensitivities_co(V,PROM_Assembly,...
    tensors_PROM,tailProperties,spineProperties,dragProperties,...
    actuTop,actuBottom,h,tmax,pActu); 
toc

% Retrieving solutions    
solShape = TI_NL_PROM.Solution.u;


%% CO-OPTIMASATION: BUILD ROM AND SOLVE EOMS ______________________________
% Parameters
h = 0.0025;
tmax = 2.0;
xi = [0.1973;-0.4848;0.2889;-0.2650;0.1222;0.3886;0.3969;0.3889];
pActu = [0.3027,2.58,1.0];

% Mesh        
df = U*xi;                    % displacement fields introduced by defects
ddf = [df(1:3:end) df(2:3:end) df(3:3:end)]; 
nodes_sv = nodes + ddf;   % nominal + d ---> defected 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);
for l=1:length(nset)
    svMesh.set_essential_boundary_condition([nset{l}],1:3,0)   
end

% build ROM
fprintf('____________________\n')
fprintf('Building ROM ... \n')

[V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
    build_PROM_3D(svMesh,nodes_sv,elements,U,USEJULIA,VOLUME,FORMULATION); 

% Solve EoMs
tic 
fprintf('____________________\n')
fprintf('Solving EoMs...\n') 
TI_NL_PROM = solve_EoMs_and_sensitivities_co(V,PROM_Assembly,...
    tensors_PROM,tailProperties,spineProperties,dragProperties,...
    actuTop,actuBottom,h,tmax,pActu); 
toc

% Retrieving solutions    
solCO = TI_NL_PROM.Solution.u;


%% CREATE FIGURE
h1 = 0.005;
h2 = 0.0025;
uTailNominal = zeros(1,tmax/h1);
uTailActu = zeros(1,tmax/h1);
uTailShape = zeros(1,tmax/h1);
uTailCO = zeros(1,tmax/h2);
timePlot1 = linspace(0,tmax-h1,tmax/h1);
timePlot2 = linspace(0,tmax-h2,tmax/h2);
x0Tail = min(nodes(:,1));

for a=1:tmax/h1
    uTailNominal(1,a) = solNom(tailProperties.tailNode*3-2,a);
    uTailActu(1,a) = solActu(tailProperties.tailNode*3-2,a);
    uTailShape(1,a) = solShape(tailProperties.tailNode*3-2,a);
end

for a=1:tmax/h2
    uTailCO(1,a) = solCO(tailProperties.tailNode*3-2,a);
end

f = figure('units','centimeters','position',[3 3 9 4]);
plot(timePlot1,x0Tail+uTailNominal(1,:),'DisplayName','Nominal','LineStyle','--','linewidth', 1.0)
hold on
plot(timePlot1,x0Tail+uTailActu(1,:),'DisplayName','AO1','LineStyle','-.', 'linewidth', 1.0)
plot(timePlot1,x0Tail+uTailShape(1,:),'DisplayName','SO5','LineStyle',':', 'linewidth', 1.0)
plot(timePlot2,x0Tail+uTailCO(1,:),'DisplayName','CO1', 'linewidth', 1.0)
xlabel('Time [s]')
ylabel('x-position tail node')
legend('Location','northwest')

exportgraphics(f,'comparison_swimming_distance_V1.pdf','Resolution',600)


