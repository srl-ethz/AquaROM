% ------------------------------------------------------------------------ 
% Actuation signal optimisation of a fish.
% 
% Last modified: 26/01/2024, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc

elementType = 'TET4';

FORMULATION = 'N1t'; % N1/N1t/N0
VOLUME = 1;         % integration over defected (1) or nominal volume (0)
USEJULIA = 1;

%% PREPARE MODEL                                                    

% DATA ____________________________________________________________________
E       = 260000;      % Young's modulus [Pa]
rho     = 1070;         % density [kg/m^3]
nu      = 0.4;        % Poisson's ratio 

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

nel = size(elements,1);
fixedPortion = 0.58;
nset = {};

fixedElements = zeros(nel,1);
for el=1:nel   
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
    if elementCenterX >= -Lx*fixedPortion
        for n=1:size(elements,2) 
            if  ~any(cat(2, nset{:}) == elements(el,n))
                nset{end+1} = elements(el,n); 
            end
        
        end   
    end
    
end


%% OPTIMIZATION PARAMETERS
dSwim = [1;0;0]; %swimming direction
h = 0.005;
tmax = 2.0;

%% OPTIMISATION ___________________________________________________________
% The actuation force needs to be chosen in the following script:
%   solve_EoMs_and_sensitivities_actu
% The initial conditions should be changed in the following script:
%   optimise_actuation_3D
%% TESTING 

AActu = [1 ;
        -1];

% bActu = [0.31;-0.1;2.6;-1.4];
bActu = [2.5;0-1.5];
% start with [0.2;2;0] % amplitude, frequency, quantity of easing
barrierParam = [1,1];
gradientWeights = [1];


tStart = tic;
[pStar,pEvo,LEvo, LwoBEvo] = optimise_actuation_3D(myElementConstructor,nset, ...
    nodes,elements,dSwim,h,tmax,AActu,bActu, ...
    'maxIteration',30, ...
    'convCrit',0.0002, ...
    'convCritCost',0.002, ...
    'barrierParam',barrierParam, ...
    'gradientWeights',gradientWeights, ...
    'gStepSize',0.01, ...
    'nResolve',6, ...
    'resolveThreshold',0.03, ...
    'USEJULIA',1);
topti = toc(tStart);
fprintf('Computation time: %.2fmin\n',topti/60)

%% TESTING 

AActu = [1 0 ;
        -1 0 ;
        0 1 ;
        0 -1];

% bActu = [0.31;-0.1;2.6;-1.4];
bActu = [0.1;0.1;0.6;0.6];
% start with [0.2;2;0] % amplitude, frequency, quantity of easing
barrierParam = [30,30,300,300];
gradientWeights = [1,1];


tStart = tic;
[pStar,pEvo,LEvo, LwoBEvo] = optimise_actuation_3D(myElementConstructor,nset, ...
    nodes,elements,dSwim,h,tmax,AActu,bActu, ...
    'maxIteration',30, ...
    'convCrit',0.02, ...
    'convCritCost',0.002, ...
    'barrierParam',barrierParam, ...
    'gradientWeights',gradientWeights, ...
    'gStepSize',0.0001, ...
    'nResolve',12, ...
    'resolveThreshold',0.1);
topti = toc(tStart);
fprintf('Computation time: %.2fmin\n',topti/60)


%% AO1
% actuation_force_6, ran with k=300. Results: [0.3027,2.58.43,1.0161]
AActu = [1 0 0;
        -1 0 0;
        0 1 0;
        0 -1 0;
        0 0 1;
        0 0 -1];
bActu = [0.31;-0.1;2.6;-1.4;1.05;0.05];
% start with [0.2;2;0] % amplitude, frequency, quantity of easing
barrierParam = [30,30,300,300,300,300];
gradientWeights = [1,1,1];


tStart = tic;
[pStar,pEvo,LEvo, LwoBEvo] = optimise_actuation_3D(myElementConstructor,nset, ...
    nodes,elements,dSwim,h,tmax,AActu,bActu, ...
    'maxIteration',30, ...
    'convCrit',0.02, ...
    'convCritCost',0.002, ...
    'barrierParam',barrierParam, ...
    'gradientWeights',gradientWeights, ...
    'gStepSize',0.01, ...
    'nResolve',12, ...
    'resolveThreshold',0.1);
topti = toc(tStart);
fprintf('Computation time: %.2fmin\n',topti/60)

%% AO 2
% actuation_force_4, ran with k=300
A = [1 0 0 0;
    -1 0 0 0;
    0 1 0 0;
    0 -1 0 0;
    0 0 1 0;
    0 0 -1 0;
    0 0 0 1;
    0 0 0 -1];
b = [1.15;-0.75;0.4;0.4;0.2;0.2;1.3;-0.7];

barrierParam = [3,3,3,3,3,3,30,30];
gradientWeights = [2,2,2,0.5];

tStart = tic;
[pStar,pEvo,LEvo,LwoBEvo] = optimise_actuation_3D(myElementConstructor,nset, ...
    nodes,elements,dSwim,h,tmax,A,b, ...
    'maxIteration',50, ...
    'convCrit',0.002, ...
    'convCritCost',0.002, ...
    'barrierParam',barrierParam, ...
    'gradientWeights',gradientWeights, ...
    'gStepSize',0.003, ...
    'nResolve',12, ...
    'resolveThreshold',0.25);
topti = toc(tStart);
fprintf('Computation time: %.2fmin\n',topti/60)

%% PLOT COST FUNCTION OVER ITERATIONS _____________________________________

fCost = figure;
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
plot(LwoBEvo)
grid on
ylabel('$$L$$','Interpreter','latex')
xlabel('Iterations')

exportgraphics(fCost,'AO1_cost_V1.pdf','Resolution',600)

%% PLOT ACTUATION SIGNAL _____________________________________

f1 = figure('units','centimeters','position',[3 3 9 4]);
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
timePlot = linspace(0,tmax-h,tmax/h);
p0 = [0.2;2;0];
%p0 = [1;0;0;1];
k=300;
actu0 = zeros(1,length(timePlot));
actuStar = zeros(1,length(timePlot));
for i = 1:length(timePlot)
    actu0(i) = actuation_signal_6(k,timePlot(i),p0);
    actuStar(i) = actuation_signal_6(k,timePlot(i),pStar);
end
plot(timePlot,actu0,'--')
hold on
plot(timePlot,actuStar)
grid on
ylabel('Actuation signal','Interpreter','latex')
xlabel('Time [s]')
legend('Initial','Optimised','Interpreter','latex','Location','northoutside','Orientation','horizontal')
hold off
exportgraphics(f1,'AO1_signal_V1.pdf','Resolution',600)

%% PLOT PARAMETERS' EVOLUTION OVER ITERATIONS _____________________________
figure('Position',[100,100,600,200])
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
subplot(1,3,1)
plot(pEvo(1,:));
ylabel('$$p_1$$','Interpreter','latex')
xlabel('Iterations')
grid on
subplot(1,3,2)
plot(pEvo(2,:));
ylabel('$$p_2$$','Interpreter','latex')
xlabel('Iterations')
grid on
subplot(1,3,3)
plot(pEvo(3,:));
ylabel('$$p_3$$','Interpreter','latex')
xlabel('Iterations')
grid on


