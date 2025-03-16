% ------------------------------------------------------------------------ 
% Actuates fish without fluid, and compare results with FEniCS results.
% 
% Last modified: 04/02/2024, Mathieu Dubied, ETH Zurich
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
nu      = 0.499;        % Poisson's ratio 

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
    MeshNominal.set_essential_boundary_condition([nset{l}],1:3,0)
end

% ASSEMBLY ________________________________________________________________
% nominal
NominalAssembly = Assembly(MeshNominal);
Mn = NominalAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( MeshNominal.nDOFs, 1);
[Kn,~] = NominalAssembly.tangent_stiffness_and_force(u0);
% store matrices
NominalAssembly.DATA.K = Kn;
NominalAssembly.DATA.M = Mn;
    
% DAMPING _________________________________________________________________
alfa = 0.912;
beta = 0.002;

Dn = alfa*Mn + beta*Kn; % Rayleigh damping 
NominalAssembly.DATA.D = Dn;
NominalAssembly.DATA.C = Dn;
Dc = NominalAssembly.constrain_matrix(Dn);

%% ROM ____________________________________________________________________

% EIGENMODES - VIBRATION MODES (VMs) ______________________________________
n_VMs = 6;  % number of vibration modes to include in the ROM

% eigenvalue problem
Kc = NominalAssembly.constrain_matrix(Kn);
Mc = NominalAssembly.constrain_matrix(Mn);
[VMn,om] = eigs(Kc, Mc, n_VMs, 'SM');
[f0n,ind] = sort(sqrt(diag(om))/2/pi);
VMn = VMn(:,ind);
for ii = 1:n_VMs
    VMn(:,ii) = VMn(:,ii)/max(sqrt(sum(VMn(:,ii).^2,2)));
end
VMn = NominalAssembly.unconstrain_vector(VMn);

% MODAL DERIVATIVES (MDs) _________________________________________________                   
[MD, MDname] = modal_derivatives(NominalAssembly, elements, VMn);

% ROM TENSORS (INTERNAL FORCES) ___________________________________________                       
% define reduced order basis
V = [VMn MD];     % reduced order basis (ROM-n)

% orthonormalize reduction basis
V = orth(V);	% ROM-n

tensors_ROMn = reduced_tensors_ROM(NominalAssembly, elements, V, USEJULIA);

% REDUCED ASSEMBLY ________________________________________________________

ROM_Assembly = ReducedAssembly(MeshNominal, V);
ROM_Assembly.DATA.M = ROM_Assembly.mass_matrix();   % reduced mass matrix 
ROM_Assembly.DATA.C = V.'*Dn*V;                   % reduced damping matrix 
ROM_Assembly.DATA.K = V.'*Kn*V;                   % reduced stiffness matrix 

%% PLOT VMs 
mod = 1;
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes, elementPlot, 0);
v1 = reshape(VMn(:,mod), 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', max(nodes(:,2)));
title(['\Phi_' num2str(mod)])

%%
% ACTUATION FORCES ____________________________________________________
actuationDirection = [1;0;0;0;0;0]; %[1;0;0]-->[1;0;0;0;0;0] (Voigt notation)

% left muscle
leftMuscle = zeros(nel,1);
for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
    if elementCenterY>0.00 &&  elementCenterX < -Lx*0.25 && elementCenterX > -Lx*0.8
        leftMuscle(el) = 1;
    end      
end
actuLeft = reduced_tensors_actuation_ROM(NominalAssembly, V, leftMuscle, actuationDirection);

% right muscle
rightMuscle = zeros(nel,1);
for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
    if elementCenterY<0.00 &&  elementCenterX < -Lx*0.25 && elementCenterX > -Lx*0.8
        rightMuscle(el) = 1;
    end    
end
actuRight = reduced_tensors_actuation_ROM(NominalAssembly, V, rightMuscle, actuationDirection);

%% PLOT ACTUATION ELEMENT IN RED __________________________________________
disp=zeros(size(nodes,1),3);
figure('units','normalized','position',[.1 .2 .4 .4])
PlotFieldonDeformedMeshActuation(nodes,elements,leftMuscle,0.8,disp,'factor',1,'color', 'k') ;
title('Left muscle')
figure('units','normalized','position',[.5 .2 .4 .4])
title('Right muscle')
PlotFieldonDeformedMeshActuation(nodes,elements,rightMuscle,0.8,disp,'factor',1,'color', 'k') ;

%% TIME INTEGRATION _______________________________________________________
% time step for integration
h = 0.005;
% initial condition: equilibrium
q0 = zeros(size(V,2),1);
qd0 = zeros(size(V,2),1);
qdd0 = zeros(size(V,2),1);

% actuation forces
B1L = actuLeft.B1;
B1R = actuRight.B1;
B2L = actuLeft.B2;
B2R = actuRight.B2;
k = 30; 

actuSignalL = @(t) k/2*(-0.2*sin(t*2*pi));    
actuSignalR = @(t) k/2*(0.2*sin(t*2*pi));

fActu = @(t,q)  k/2*(-0.2*sin(t*2*pi))*(B1L+B2L*q) + ...
                k/2*(0.2*sin(t*2*pi))*(B1R+B2R*q);

% instantiate object for nonlinear time integration
TI_NL_ROM = ImplicitNewmark('timestep',h,'alpha',0.005);

% modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_actuation(q, ...
    qd,qdd,t,ROM_Assembly,fActu,actuLeft,actuRight,actuSignalL,actuSignalR);

% nonlinear Time Integration
tmax = 2; 
tic
TI_NL_ROM.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
toc

%% PLOT ___________________________________________________________________
[~, spineElements, ~, nodeIdxPosInElements] = find_spine_TET4(elements,nodes);
[tailNode, ~, ~] = find_tail(elements,nodes,spineElements,nodeIdxPosInElements);
        
uTail = zeros(3,tmax/h);
timePlot = linspace(0,tmax-h,tmax/h);
x0Tail = min(nodes(:,1));
for a=1:tmax/h
    uTail(:,a) = V(tailNode*3-2:tailNode*3,:)*TI_NL_ROM.Solution.q(:,a);
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


