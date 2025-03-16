% ------------------------------------------------------------------------ 
% 3D actuation signal optimization of a fish.
% 
% Last modified: 10/11/2023, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc

elementType = 'TET4';

FORMULATION = 'N1t'; % N1/N1t/N0
VOLUME = 1;         % integration over defected (1) or nominal volume (0)

USEJULIA = 0;

%% PREPARE MODEL                                                    

% DATA ____________________________________________________________________
E       = 2600000;      % Young's modulus [Pa]
rho     = 1070;         % density [kg/m^3]
nu      = 0.499;        % Poisson's ratio 

% material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	    % set "false" for plane_strain

% element
switch elementType
    case 'TET4'
        myElementConstructor = @()Tet4Element(myMaterial);
end
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
% for l=1:length(nset)
%     MeshNominal.set_essential_boundary_condition([nset{l}],2,0)   %
%     fixed head
%     MeshNominal.set_essential_boundary_condition([nset{l}],2,0)     % head on "rail"
% end




%% OPTIMIZATION PARAMETERS
dSwim = [1;0;0]; %swimming direction
h = 0.005;
tmax = 1.0;

%% OPTIMISATION TEST 1 ____________________________________________________

% A = [1;-1];
% b = [0.8;1.2];
A = [1 0 0;
    -1 0 0;
    0 1 0;
    0 -1 0;
    0 0 1;
    0 0 -1];
b = [1.1;-0.8;0.5;0.1;0.2;0.1];

tStart = tic;
[xiStar,xiEvo,LrEvo, ~] = optimise_actuation_3D(myElementConstructor,nset, ...
    nodes,elements,dSwim,h,tmax,A,b, ...
    'maxIteration',18, ...
    'convCrit',0.001, ...
    'convCritCost',0.0001, ...
    'barrierParam',20, ...
    'gStepSize',0.02, ...
    'nResolve',8, ...
    'resolveThreshold',0.1);
topti = toc(tStart);
fprintf('Computation time: %.2fmin\n',topti/60)



%% PLOT COST FUNCTION OVER ITERATIONS _____________________________________

figure
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
plot(LrEvo)
grid on
ylabel('$$L_r$$','Interpreter','latex')
xlabel('Iterations')



%% PLOT PARAMETERS' EVOLUTION OVER ITERATIONS _____________________________
figure('Position',[100,100,600,200])
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
subplot(1,3,1)
plot(xiEvo(1,:));
ylabel('$$p_1$$','Interpreter','latex')
xlabel('Iterations')
grid on
subplot(1,3,2)
plot(xiEvo(2,:));
ylabel('$$p_2$$','Interpreter','latex')
xlabel('Iterations')
grid on
subplot(1,3,3)
plot(xiEvo(3,:));
ylabel('$$p_3$$','Interpreter','latex')
xlabel('Iterations')
grid on


%% COST COMPUTATION ON FINAL ROM __________________________________________
%xiStar4 = [0.2;0.1;0.5;0.2;0.1]
xiFinal = xiStar;

% tmax=4;
% h=0.0025

% shape-varied mesh 
df = U*xiFinal;                       % displacement field introduced by shape variations
dd = [df(1:2:end) df(2:2:end)];   % rearrange as two columns matrix
nodes_sv = nodes + dd;          % nominal + dd ---> shape-varied nodes 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);
svMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)

% (P)ROM creation
FORMULATION = 'N1';VOLUME = 1; USEJULIA = 0;FOURTHORDER = 0; ACTUATION = 1;
[V,PROM_Assembly,tensors_PROM,tensors_hydro_PROM, tensors_topMuscle_PROM, tensors_bottomMuscle_PROM] = ...
        build_PROM(svMesh,nodes,elements,U,FORMULATION,VOLUME,USEJULIA,FOURTHORDER,ACTUATION);
       
% solve EoMs
TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax, ...
     'ACTUATION', ACTUATION,'topMuscle',tensors_topMuscle_PROM,'bottomMuscle',tensors_bottomMuscle_PROM);
% TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax);
eta = TI_NL_PROM.Solution.q;
etad = TI_NL_PROM.Solution.qd;
N = size(eta,2);

% compute cost Lr without barrier functions (no constraints, to obtain the
% the cost stemming from the hydrodynamic force only)
dr = reduced_constant_vector(dSwim,V);
AFinal = [];     % no constraint
bFinal= [];      % no constraint
barrierParam = 10;
Lr = reduced_cost_function_w_constraints(N,tensors_hydro_PROM,eta,etad,xiFinal,dr,AFinal,bFinal,barrierParam);
%fprintf('The cost function w/o constraint is: %.4f\n',Lr)

%% PLOTs __________________________________________________________________

% find a specific result node and corresponding DOF
tailNodeDOFS = MeshNominal.get_DOF_from_location([Lx, 0]);
tailNodeDOF = tailNodeDOFS(1); % y-direction
% time axis
tplot=linspace(0,tmax,tmax/h+1);

% plot
figure('units','normalized','position',[.1 .1 .8 .6],'name','Vertical displacement of the tail node')
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

plot(tplot,TI_NL_PROM.Solution.u(tailNodeDOF,1:end-1)*100, "--")

ylabel('$$u_y \mbox{ [cm]}$$','Interpreter','latex')
xlabel('Time [s]')
set(gca,'FontName','ComputerModern');
grid on
%legend({'FOM','FOM-t','ROM','PROM'}, 'Location', 'eastoutside','Orientation','vertical')
hold off



%% ANIMATIONS _____________________________________________________________
