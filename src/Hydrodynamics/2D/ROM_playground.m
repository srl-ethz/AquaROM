% ------------------------------------------------------------------------ 
% Playground to test fish locomotion forces in ROM
% Used element type: TRI3.
% 
% Last modified: 18/10/2023, Mathieu Dubied, ETH Zurich
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

% PROM parameters
xi1 = 0.2;

% MESH ____________________________________________________________________

% nominal mesh
switch upper( whichModel )
    case 'ABAQUS'
        filename = 'naca0012_76el_2';
        [nodes, elements, ~, elset] = mesh_ABAQUSread(filename);
end

nodes = nodes;

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of airfoil
Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of airfoil

%%
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
for l=1:length(nset)
    MeshNominal.set_essential_boundary_condition([nset{l}],1:2,0)   %
    % fixed head
    % MeshNominal.set_essential_boundary_condition([nset{l}],2,0)     % head on "rail"
end

%% VISUALISATION __________________________________________________________
% plot nominal mesh
elementPlot = elements(:,1:3); 
figure('units','normalized','position',[.2 .1 .6 .4],'name','Nominal mesh')
PlotMesh(nodes, elementPlot, 1);


%% ROM CONSTRUCTION _______________________________________________________
ACTUATION =1;
mTilde = 10;
[V,ROM_Assembly,tensors_ROM,tailProp,spineProp,dragProp,actuTop,actuBottom] = ...
    build_ROM_non_optimisation(MeshNominal,nodes,elements,mTilde,USEJULIA,ACTUATION);

%% TIME INTEGRATION _______________________________________________________

% time step for integration
h = 0.01;
tmax = 2.0; 

% initial condition: equilibrium
fprintf('solver')
tic
q0 = zeros(size(V,2),1);        % q, qd are reduced order DOFs
qd0 = zeros(size(V,2),1);
% qd0(1:2:end) = -0.2;
qdd0 = zeros(size(V,2),1);

% actuation force
B1T = actuTop.B1;
B1B = actuBottom.B1;
B2T = actuTop.B2;
B2B = actuBottom.B2;
k=1;

actuSignalT = @(t) k/2*(1-(1+0.02*sin(t*2*pi)));    % to change below as well if needed
actuSignalB = @(t) k/2*(1-(1-0.02*sin(t*2*pi)));

fActu = @(t,q)  k/2*(1-(1+0.02*sin(t*2*pi)))*(B1T+B2T*q) + ...
                k/2*(1-(1-0.02*sin(t*2*pi)))*(B1B+B2B*q);

% tail pressure force 
A = tailProp.A;
B = tailProp.B;
R = [0 -1 0 0 0 0;
     1 0 0 0 0 0;
     0 0 0 -1 0 0;
     0 0 1 0 0 0;
     0 0 0 0 0 1;
     0 0 0 0 -1 0];     % 90 degrees rotation counterclock-wise
wTail = tailProp.w;
VTail = tailProp.V;
nodesInPos = V.'*reshape(nodes.',[],1);     % initial node position expressed in the ROM
x0 = reshape(nodes.',[],1);  
                   
fTail = @(q,qd)  0.5*mTilde*2*wTail^3*VTail.'*(dot(A*VTail*qd,R*B*(x0(tailProp.iDOFs)+VTail*q))).^2* ...
                    B*(x0(tailProp.iDOFs)+VTail*q);

% Constant force test
% forceTest = zeros(size(reshape(nodes.',[],1)));
% forceTest(1:2:end)=2;
% fTail = @(q,qd)  V.'*forceTest;

% spine change in momentum
Txx = spineProp.tensors.Txx;
TxV = spineProp.tensors.TxV;
TVx = spineProp.tensors.TVx;
TVV = spineProp.tensors.TVV;

fSpine = @(q,qd,qdd) double(Txx)*qdd ...
    + double(ttv(ttv(TxV,q,3),qdd,2) ...
    + ttv(ttv(TVx,q,3),qdd,2) ...
    + ttv(ttv(ttv(TVV,q,4),q,3),qdd,2) ...
    + ttv(ttv(TVx,qd,3),qd,2) ...
    + ttv(ttv(ttv(TVV,q,4),qd,3),qd,2) ...
    + ttv(ttv(TxV,qd,3),qd,2) ...
    + ttv(ttv(ttv(TVV,qd,4),q,3),qd,2));

% drag force
T3 = dragProp.tensors.Tr3;
fDrag = @(qd)  0.5*double(ttv(ttv(T3,qd,3),qd,2));

% instantiate object for nonlinear time integration
TI_NL_ROM = ImplicitNewmark('timestep',h,'alpha',0.005,'MaxNRit',60,'RelTol',1e-6);

% modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_actu_hydro(q,qd, ...
    qdd,t,ROM_Assembly,tensors_ROM,fActu,fTail,fSpine,fDrag,actuTop,actuBottom,actuSignalT,actuSignalB,tailProp,spineProp,dragProp,R,x0);

% nonlinear Time Integration
TI_NL_ROM.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_ROM.Solution.u = V * TI_NL_ROM.Solution.q; % get full order solution

TI_NL_ROM.Solution.ud = V * TI_NL_ROM.Solution.qd; % get full order solution

TI_NL_ROM.Solution.udd = V * TI_NL_ROM.Solution.qdd; % get full order solution
toc

%% PLOT ___________________________________________________________________
tailNode = 1;

initialPosFOM = reshape(nodes.',[],1);
figure(Position=[300,150,300,300])
plot(initialPosFOM(tailNode*2)+TI_NL_ROM.Solution.u(tailNode*2,:))
title('tail y position')
figure(Position=[600,150,300,300])
plot(TI_NL_ROM.Solution.ud(tailNode*2,:))
title('tail y velocity')

%% PLOT ___________________________________________________________________

initialPosFOM = reshape(nodes.',[],1);
figure(Position=[300,150,300,300])
plot(initialPosFOM(tailNode*2-1)+TI_NL_ROM.Solution.u(tailNode*2-1,:))
title('tail x position')
figure(Position=[600,150,300,300])
plot(TI_NL_ROM.Solution.ud(tailNode*2-1,:))
title('tail x velocity')
%% lines at t=0.5, 1.0 etc
timeSnapshot = 1:50:501;
figure(Position=[600,150,300,300])
plot(timeSnapshot,TI_NL_ROM.Solution.ud(tailNode*2-1,timeSnapshot))
title('tail velocity')

%%
headNode = 2;

initialPosFOM = reshape(nodes.',[],1);
figure(Position=[300,150,300,300])
plot(initialPosFOM(headNode*2)+TI_NL_ROM.Solution.u(headNode*2,:))
title('tail position')
figure(Position=[600,150,300,300])
plot(TI_NL_ROM.Solution.ud(headNode*2,:))
title('tail velocity')

%% CHECK TAIL FORCE _______________________________________________________
upTo = 500;
force = zeros(6,upTo);
pos = zeros(6,upTo);
tvec = zeros(6,upTo);
for i=1:upTo
    force(:,i) = 0.5*mTilde*2*wTail^3*(dot(A*VTail*TI_NL_ROM.Solution.qd(:,i),R*B*(x0(tailProp.iDOFs)+VTail*TI_NL_ROM.Solution.q(:,i)))).^2*...
                    B*(x0(tailProp.iDOFs)+VTail*TI_NL_ROM.Solution.q(:,i));
    pos(:,i)=VTail*TI_NL_ROM.Solution.q(:,i);
    tvec(:,i) = B*(x0(tailProp.iDOFs)+VTail*TI_NL_ROM.Solution.q(:,i));
end


figure
plot(force(5,:))
hold on
plot(force(6,:))
plot(pos(6,:)*10)
legend('force x','force y', 'y position*10')
hold off

figure
plot(force(5,:))
hold on
plot(force(6,:))
plot(tvec(5,:))
plot(tvec(6,:))
legend('force x','force y', 't x','t_y')
hold off


%% CHECK SPINE FORCE AT TAIL ELEMENT ______________________________________
upTo = 300;
force = zeros(6,upTo);
pos = zeros(6,upTo);
tvec = zeros(6,upTo);

for i=1:upTo
    force(:,i) = -mTilde*wTail*( ...
        (A*VTail*TI_NL_ROM.Solution.qdd(:,i)).'*R*B*(x0(tailProp.iDOFs)+VTail*TI_NL_ROM.Solution.q(:,i))*R*B*(x0(tailProp.iDOFs)+VTail*TI_NL_ROM.Solution.q(:,i)) ...
        +(A*VTail*TI_NL_ROM.Solution.qd(:,i)).'*R*B*VTail*TI_NL_ROM.Solution.qd(:,i)*R*B*(x0(tailProp.iDOFs)+VTail*TI_NL_ROM.Solution.q(:,i))...
        +(A*VTail*TI_NL_ROM.Solution.qd(:,i)).'*R*B*(x0(tailProp.iDOFs)+VTail*TI_NL_ROM.Solution.q(:,i))*R*B*VTail*TI_NL_ROM.Solution.qd(:,i));
    pos(:,i)=VTail*TI_NL_ROM.Solution.q(:,i);
    % tvec(:,i) = B*(x0(tailProp.iDOFs)+VTail*TI_NL_ROM.Solution.q(:,i));
end


figure
plot(force(5,:))
hold on
plot(force(6,:)*0)
plot(pos(6,:)*10)
legend('force x','0-line', 'y position*10')
hold off

%%
figure
plot(force(5,:))
hold on
plot(force(6,:))
plot(tvec(5,:))
plot(tvec(6,:))
legend('force x','force y', 't x','t_y')
hold off

%% ANIMATIONS _____________________________________________________________

TI_NL_PROM.Solution.u = V * TI_NL_ROM.Solution.q; % get full order solution

elementPlot = elements(:,1:3); 
nel = size(elements,1);
actuationDirection = [1;0;0];%[1;0]-->[1;0;0] (Voigt notation)

% top muscle
topMuscle = zeros(nel,1);
for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2))/3;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1))/3;
    if elementCenterY>0.00 &&  elementCenterX < -Lx*0.25 && elementCenterX > -Lx*0.8
        topMuscle(el) = 1;
    end    
end

% bottom muscle
bottomMuscle = zeros(nel,1);
for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2))/3;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1))/3;
    if elementCenterY<0.00 &&  elementCenterX < -Lx*0.25 && elementCenterX > -Lx*0.8
        bottomMuscle(el) = 1;
    end    
end

actuationValues = zeros(size(TI_NL_PROM.Solution.u,2),1);
for t=1:size(TI_NL_PROM.Solution.u,2)
    actuationValues(t) = -1*(actuSignalT(t*h)*2/k-1);
end

actuationValues2 = zeros(size(TI_NL_PROM.Solution.u,2),1);
for t=1:size(TI_NL_PROM.Solution.u,2)
    actuationValues2(t) = -1*(actuSignalB(t*h)*2/k-1);
end

AnimateFieldonDeformedMeshActuation2Muscles(nodes, elementPlot,topMuscle,actuationValues,...
    bottomMuscle,actuationValues2,TI_NL_PROM.Solution.u, ...
    'factor',1,'index',1:2,'filename','result_video','framerate',1/h)


%% ANIMATION ______________________________________________________________
elementPlot = elements(:,1:3); % plot only corners (otherwise it's a mess)
AnimateFieldonDeformedMesh(nodes, elementPlot,TI_NL_ROM.Solution.u, ...
    'factor',1,'index',1:2,'filename','result_video','framerate',1/h)

