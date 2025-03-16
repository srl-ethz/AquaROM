% ------------------------------------------------------------------------ 
% Test of the FOM.
% 
% Last modified: 24/02/2024, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc

elementType = 'TET4';

FORMULATION = 'N1t'; % N1/N1t/N0
VOLUME = 1;         % integration over defected (1) or nominal volume (0)

USEJULIA = 1;

%% PREPARE MODEL __________________________________________________________                                                   

% DATA ____________________________________________________________________
E       = 260000;      % Young's modulus [Pa]
rho     = 1070;         % density [kg/m^3]
nu      = 0.4;        % Poisson's ratio 

% material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	    % set "false" for plane_strain
myElementConstructor = @()Tet4Element(myMaterial);

% MESH ____________________________________________________________________

filename = '3d_fish_for_mike';
[nodes, elements, ~, elset] = mesh_ABAQUSread(filename);
nel = size(elements,1);

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
desiredFixedPoint = - 0.3*Lx;
fixedPoint = find_fixed_point(nodes,desiredFixedPoint);
nel = size(elements,1);
nset = {};
marginFixedPoint = 0.02;
for el=1:nel   
    for n=1:size(elements,2)
        if  ~any(cat(2, nset{:}) == elements(el,n))
            nset{end+1} = elements(el,n);
        end
    end   
end
[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);

% for testing
U = [z_tail,y_head];

for l=1:length(nset)
    MeshNominal.set_essential_boundary_condition([nset{l}],3,0)
end   

desiredFixedPoint = - 0.45*Lx;
fixedPoint = find_fixed_point(nodes,desiredFixedPoint);
nel = size(elements,1);
nset2 = {};
marginFixedPoint = 0.02;
for el=1:nel   
    for n=1:size(elements,2)
        if  nodes(elements(el,n),1)> fixedPoint-marginFixedPoint && ...
                nodes(elements(el,n),1)< fixedPoint+marginFixedPoint && ...
                ~any(cat(2, nset2{:}) == elements(el,n))
            nset2{end+1} = elements(el,n);
        end
    end   
end

% for l=1:length(nset2)
%     MeshNominal.set_essential_boundary_condition([nset2{l}],2,0)
% %     MeshNominal.set_essential_boundary_condition([nset{l}],1:3,0)
% end   

%% SIMULATION PARAMETERS __________________________________________________
h = 0.005;
tmax = 1.2;

%% FOM ____________________________________________________________________
tStartFOM = tic;

% build PROM
fprintf('____________________\n')
fprintf('Building FOM ... \n')
[Assembly,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
build_FOM_3D(MeshNominal,nodes,elements);   
% %% EoMs: ACTUATION TEST ___________________________________________________
% % SETTING UP ______________________________________________________________
% % initial conditions
% nUncDOFs = size(Assembly.Mesh.EBC.unconstrainedDOFs,2);
% nDOFs = Assembly.Mesh.nDOFs;
% q0 = zeros(nUncDOFs,1);
% qd0 = zeros(nUncDOFs,1);
% qdd0 = zeros(nUncDOFs,1);
% 
% % forces
% B1T = actuTop.B1;
% B1B = actuBottom.B1;
% B2T = actuTop.B2;
% B2B = actuBottom.B2;
% k=8; 
% 
% actuSignalT = @(t) k/2*(-0.2*sin(t*2*pi));    % to change below as well if needed
% actuSignalB = @(t) k/2*(0.2*sin(t*2*pi));
% 
% fActu = @(t,q)  k/2*(-0.2*sin(t*2*pi))*(B1T+B2T*q) + ...
%                 k/2*(0.2*sin(t*2*pi))*(B1B+B2B*q);
%             
% % NONLINEAR TIME INTEGRATION ______________________________________________
%     
% % instantiate object for nonlinear time integration
% TI_NL_FOM = ImplicitNewmark('timestep',h,'alpha',0.005,'MaxNRit',60,'RelTol',1e-6);
% 
% % modal nonlinear Residual evaluation function handle
% Residual_NL = @(q,qd,qdd,t)residual_nonlinear_actu(q,qd,qdd, ...
%     t,Assembly,fActu,actuTop,actuBottom,actuSignalT,actuSignalB);
% 
% % time integration 
% TI_NL_FOM.Integrate(q0,qd0,qdd0,tmax,Residual_NL);
% TI_NL_FOM.Solution.u = Assembly.unconstrain_vector(TI_NL_FOM.Solution.q);
% 
% fprintf('Time needed to solve the problem using FOM: %.2fsec\n',toc(tStartFOM))
% timeFOM = toc(tStartFOM);

%% EoMs: HYDRODYNAMIC TEST ________________________________________________
% SETTING UP ______________________________________________________________
% initial conditions
if isempty(Assembly.Mesh.EBC)
     nUncDOFs = Assembly.Mesh.nDOFs;  % no constraints
     fprintf('Info: no constraints \n')
else
    nUncDOFs = size(Assembly.Mesh.EBC.unconstrainedDOFs,2);% some constraints 
end

nDOFs = Assembly.Mesh.nDOFs;
q0 = zeros(nUncDOFs,1);
qd0 = zeros(nUncDOFs,1);
qdd0 = zeros(nUncDOFs,1);

% forces
B1T = actuTop.B1;
B1B = actuBottom.B1;
B2T = actuTop.B2;
B2B = actuBottom.B2;
k=1.0; 

actuSignalT = @(t) k/2*(-0.2*sin(t*4*pi));    % to change below as well if needed
actuSignalB = @(t) k/2*(0.2*sin(t*4*pi));

fActu = @(t,q)  k/2*(-0.2*sin(t*4*pi))*(B1T+B2T*q) + ...
                k/2*(0.2*sin(t*4*pi))*(B1B+B2B*q);
            
fTail = @(q,qd) tail_force_FOM(q,qd,Assembly,elements,tailProperties);

fSpine = @(q,qd,qdd) spine_force_FOM(q,qd,qdd,Assembly,elements,spineProperties);

fDrag = @(qd) drag_force_FOM(qd,dragProperties);
       
% NONLINEAR TIME INTEGRATION ______________________________________________
    
% instantiate object for nonlinear time integration
TI_NL_FOM = ImplicitNewmark('timestep',h,'alpha',0.005,'MaxNRit',60,'RelTol',1e-6);

% modal nonlinear Residual evaluation function handle
Residual_NL = @(q,qd,qdd,t)residual_nonlinear_actu_hydro(q,qd,qdd, ...
    t,Assembly,fActu,fTail,fSpine,fDrag,actuTop,actuBottom,actuSignalT,actuSignalB);

% time integration 
TI_NL_FOM.Integrate(q0,qd0,qdd0,tmax,Residual_NL);
TI_NL_FOM.Solution.u = Assembly.unconstrain_vector(TI_NL_FOM.Solution.q);

fprintf('Time needed to solve the problem using FOM: %.2fsec\n',toc(tStartFOM))
timeFOM = toc(tStartFOM);


%% PLOT ___________________________________________________________________
% unconstrained solution
sol = Assembly.unconstrain_vector(TI_NL_FOM.Solution.q);
uTail = zeros(3,tmax/h);
timePlot = linspace(0,tmax-h,tmax/h);
x0Tail = min(nodes(:,1));
for a=1:tmax/h
    uTail(:,a) = sol(tailProperties.tailNode*3-2:tailProperties.tailNode*3,a);
end

figure
subplot(2,1,1);
plot(timePlot,x0Tail+uTail(1,:),'DisplayName','k=0')
hold on
xlabel('Time [s]')
ylabel('x-position tail node')
legend('Location','northwest')
subplot(2,1,2);
plot(timePlot,uTail(2,:),'DisplayName','k=0')
hold on
xlabel('Time [s]')
ylabel('y-position tail node')
legend('Location','southwest')
drawnow

%%reshape(PROM_Assembly.Mesh.nodes.',[],1);

%% ANIMATION ______________________________________________________________
elementPlot = elements(:,1:4); 
nel = size(elements,1);

% top muscle
topMuscle = zeros(nel,1);

for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
    if elementCenterY>0.00 &&  elementCenterX < -Lx*0.25 && elementCenterX > -Lx*0.85
        topMuscle(el) = 1;%0.25,0.85
    end    
end

% bottom muscle
bottomMuscle = zeros(nel,1);
for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
    if elementCenterY<0.00 &&  elementCenterX < -Lx*0.25 && elementCenterX > -Lx*0.85
        bottomMuscle(el) = 1;
    end    

end

actuationValues = zeros(size(TI_NL_FOM.Solution.u,2),1);
for t=1:size(TI_NL_FOM.Solution.u,2)
    actuationValues(t) = 0.2;
end

actuationValues2 = zeros(size(TI_NL_FOM.Solution.u,2),1);
for t=1:size(TI_NL_FOM.Solution.u,2)
    actuationValues2(t) = -0.2;
end
sol = TI_NL_FOM.Solution.u(:,1:2:end);
AnimateFieldonDeformedMeshActuation2Muscles(nodes, elementPlot,topMuscle,actuationValues,...
    bottomMuscle,actuationValues2,sol, ...
    'factor',1,'index',1:3,'filename','result_video','framerate',1/h*0.5)






