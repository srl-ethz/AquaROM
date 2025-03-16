% ------------------------------------------------------------------------ 
% 2D optimization of a fish.
% 
% Last modified: 19/10/2023, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc

whichModel = 'ABAQUS';
elementType = 'TRI3';

FORMULATION = 'N1'; % N1/N1t/N0
VOLUME = 1;         % integration over defected (1) or nominal volume (0)

USEJULIA = 0;

%% PREPARE MODEL                                                    

% DATA ____________________________________________________________________
E       = 2600000;      % Young's modulus [Pa]
rho     = 1070;         % density [kg/m^3]
nu      = 0.499;        % Poisson's ratio 
thickness = .1;         % [m] out-of-plane thickness

% material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	    % set "false" for plane_strain

% element
switch elementType
    case 'TRI3'
        myElementConstructor = @()Tri3Element(thickness, myMaterial);
end
% MESH ____________________________________________________________________

% nominal mesh
switch upper( whichModel )
    case 'ABAQUS'
        filename = '2d_rectangle_120el';
        [nodes, elements, ~, elset] = mesh_ABAQUSread(filename);
end

nodes = nodes*0.01;

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of airfoil
Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of airfoil

% plot nominal mesh
elementPlot = elements(:,1:3); 
figure('units','normalized','position',[.2 .1 .6 .4],'name','Nominal mesh with element and node indexes')
PlotMesh(nodes, elementPlot, 0);

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


%% SHAPE VARIATIONS _______________________________________________________

% % (1) thinner airfoil 
% nodes_projected = [nodes(:,1), nodes(:,2)*0];   % projection on x-axis
% yDif = nodes_projected(:,2) - nodes(:,2);       % y-difference projection vs nominal
% thinAirfoil = zeros(numel(nodes),1);            % create a single long vectors [x1 y1 x2 y2 ...]^T
% thinAirfoil(2:2:end) = yDif;                    % fill up all y-positions
% 
% % shape variations basis
% U = thinAirfoil;    % shape variations basis
% 


% shape variations 
[thinFish,shortFish,linearTail,longTail,shortTail,linearHead,longHead,shortHead] = shape_variations_2D(nodes,Lx,Ly);
U = thinFish;

% plot the two meshes
xiPlot = ones(size(U,2),1)*0.6;
f1 = figure('units','centimeters','position',[3 3 15 7],'name','Shape-varied mesh');
elementPlot = elements(:,1:3); hold on 
PlotMesh(nodes, elementPlot, 0); 
v1 = reshape(U*xiPlot, 2, []).';
S = 1;
hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
axis equal; grid on; box on; set(hf{1},'FaceAlpha',.7); drawnow
set(f1,'PaperUnits','centimeters');
set(f1,'PaperPositionMode','auto');
% set(f1,'PaperSize',[10 3.5]); % Canvas Size
set(f1,'Units','centimeters');
% 
%%
xiFinal = [0.5];
df = U*xiFinal;                       % displacement field introduced by shape variations
dd = [df(1:2:end) df(2:2:end)];   % rearrange as two columns matrix
nodes_sv = nodes + dd;          % nominal + dd ---> shape-varied nodes 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);

%%
% plot the two meshes
xiPlot = xiFinal;
f1 = figure('units','centimeters','position',[3 3 15 7],'name','Shape-varied mesh');
elementPlot = elements(:,1:3); hold on 
PlotMesh(nodes, elementPlot, 0); 
v1 = reshape(U*xiPlot, 2, []).';
S = 1;
hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
axis equal; grid on; box on; set(hf{1},'FaceAlpha',.7); drawnow
set(f1,'PaperUnits','centimeters');
set(f1,'PaperPositionMode','auto');
% set(f1,'PaperSize',[10 3.5]); % Canvas Size
set(f1,'Units','centimeters');
% 

%% OPTIMIZATION PARAMETERS
dSwim = [1;0]; %swimming direction
h = 0.01;
tmax = 1.5;

%% OPTIMISATION TEST 1 ____________________________________________________
% A = [1 0;
%     -1 0;
%     0 1;
%     0 -1];
% b = [0.2;0.2;0.2;0.2];
A = [1;
    -1];
b = [0.2;0.2];

tStart = tic;
[xiStar,xiEvo,LrEvo] = optimise_shape_2D(myElementConstructor,nset, ...
    nodes_sv,elements,U,dSwim,h,tmax,A,b,'maxIteration',15,'convCrit',0.002,'barrierParam',100,'gStepSize', 1e-60,'nRebuild',3);
topti = toc(tStart);
fprintf('Computation time: %.2fmin\n',topti/60)


%% VISUALIZATION __________________________________________________________

% shape-varied mesh 
df = U*xiStar;                       % displacement field introduced by shape variations
dd = [df(1:2:end) df(2:2:end)];   % rearrange as two columns matrix
nodes_sv = nodes + dd;          % nominal + dd ---> shape-varied nodes 
svMesh = Mesh(nodes_sv);
svMesh.create_elements_table(elements,myElementConstructor);
svMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)

% plot the two meshes
figure('units','normalized','position',[.2 .3 .6 .4],'name','Shape-varied mesh');
elementPlot = elements(:,1:3); hold on 
PlotMesh(nodes, elementPlot, 0); 
v1 = reshape(U*xiStar, 2, []).';
S = 1;
hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
axis equal; grid on; box on; set(hf{1},'FaceAlpha',.7); drawnow
%%
% plot the two meshes
xiPlot = xiStar;
f1 = figure('units','centimeters','position',[3 3 15 7],'name','Shape-varied mesh');
elementPlot = elements(:,1:3); hold on 
PlotMesh(nodes, elementPlot, 0); 
v1 = reshape(U*xiPlot, 2, []).';
S = 1;
hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
axis equal; grid on; box on; set(hf{1},'FaceAlpha',.7); drawnow
set(f1,'PaperUnits','centimeters');
set(f1,'PaperPositionMode','auto');
% set(f1,'PaperSize',[10 3.5]); % Canvas Size
set(f1,'Units','centimeters');
% 

%% PLOT COST FUNCTION OVER ITERATIONS _____________________________________

figure
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
plot(LrEvo)
grid on
ylabel('$$L_r$$','Interpreter','latex')
xlabel('Iterations')

%% PLOT XI PARAMATER OVER ITERATIONS ______________________________________

figure
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
plot3(xiEvo(1,:),xiEvo(2,:),linspace(1,size(xiEvo,2),size(xiEvo,2)));
grid on
xlabel('$$\xi_1$$','Interpreter','latex')
ylabel('$$\xi_2$$','Interpreter','latex')
zlabel('Iterations')

%%
figure
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
plot(xiEvo(2,:));
grid on
xlabel('$$\xi_1$$','Interpreter','latex')
ylabel('Iterations')




%% COST COMPUTATION ON FINAL ROM __________________________________________
xiStar = [0.2;0.2]
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

TI_NL_PROM.Solution.u = V * TI_NL_PROM.Solution.q; % get full order solution

% ROM TENSORS - ACTUATION FORCES __________________________________________

elementPlot = elements(:,1:3); 
nel = size(elements,1);
actuationDirection = [1;0;0];%[1;0]-->[1;0;0] (Voigt notation)

% top muscle
topMuscle = zeros(nel,1);
for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2))/3;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1))/3;
    if elementCenterY>0.00 &&  elementCenterX > Lx*0.25
        topMuscle(el) = 1;
    end    
end

% bottom muscle
bottomMuscle = zeros(nel,1);
for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2))/3;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1))/3;
    if elementCenterY<0.00 &&  elementCenterX > Lx*0.25
        bottomMuscle(el) = 1;
    end    
end


actuationValues = zeros(size(TI_NL_PROM.Solution.u,2),1);
for t=1:size(TI_NL_PROM.Solution.u,2)
    actuationValues(t) = 1+0.005*sin(t*h*2*pi);
end

actuationValues2 = zeros(size(TI_NL_PROM.Solution.u,2),1);
for t=1:size(TI_NL_PROM.Solution.u,2)
    actuationValues2(t) = 1-0.005*sin(t*h*2*pi);
end

AnimateFieldonDeformedMeshActuation2Muscles(nodes_sv, elementPlot,topMuscle,actuationValues,...
    bottomMuscle,actuationValues2,TI_NL_PROM.Solution.u, ...
    'factor',1,'index',1:2,'filename','result_video','framerate',1/h)

%% COST COMPUTATION ON INITIAL ROM ________________________________________
% % 88.57
% xiFinal = zeros(size(U,2),1);
% 
% % shape-varied mesh 
% df = U*xiFinal;                       % displacement field introduced by shape variations
% dd = [df(1:2:end) df(2:2:end)];   % rearrange as two columns matrix
% nodes_sv = nodes + dd;          % nominal + dd ---> shape-varied nodes 
% svMesh = Mesh(nodes_sv);
% svMesh.create_elements_table(elements,myElementConstructor);
% svMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)
% 
% % (P)ROM creation
% FORMULATION = 'N1';VOLUME = 1; USEJULIA = 0;FOURTHORDER = 0; ACTUATION = 1;
% [V,PROM_Assembly,tensors_PROM,tensors_hydro_PROM, tensors_topMuscle_PROM, tensors_bottomMuscle_PROM] = ...
%         build_PROM(svMesh,nodes,elements,U,FORMULATION,VOLUME,USEJULIA,FOURTHORDER,ACTUATION);
%         
% % solve EoMs
% TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax,...
%                     'ACTUATION', ACTUATION,'topMuscle',tensors_topMuscle_PROM,'bottomMuscle',tensors_bottomMuscle_PROM);
% 
% eta = TI_NL_PROM.Solution.q;
% etad = TI_NL_PROM.Solution.qd;
% N = size(eta,2);
% 
% % compute cost Lr without barrier functions (no constraints, to obtain the
% % the cost stemming from the hydrodynamic force only)
% dr = reduced_constant_vector(dSwim,V);
% AFinal = [];     % no constraint
% bFinal= [];      % no constraint
% barrierParam = 10;
% Lr = reduced_cost_function_w_constraints(N,tensors_hydro_PROM,eta,etad,xiFinal,dr,AFinal,bFinal,barrierParam);
% fprintf('The cost function w/o constraint is: %.4f\n',Lr)

% %% PLOTs __________________________________________________________________
% 
% % find a specific result node and corresponding DOF
% tailNodeDOFS = MeshNominal.get_DOF_from_location([Lx, 0]);
% tailNodeDOF = tailNodeDOFS(1); % y-direction
% % time axis
% tplot=linspace(0,tmax,tmax/h+1);
% 
% % plot
% figure('units','normalized','position',[.1 .1 .8 .6],'name','Vertical displacement of the tail node')
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% set(groot,'defaultAxesTickLabelInterpreter','latex'); 
% 
% plot(tplot,TI_NL_PROM.Solution.u(tailNodeDOF,1:end-1)*100, "--")
% 
% ylabel('$$u_y \mbox{ [cm]}$$','Interpreter','latex')
% xlabel('Time [s]')
% set(gca,'FontName','ComputerModern');
% grid on
% %legend({'FOM','FOM-t','ROM','PROM'}, 'Location', 'eastoutside','Orientation','vertical')
% hold off
