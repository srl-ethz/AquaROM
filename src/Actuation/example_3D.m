% ------------------------------------------------------------------------ 
% Script to test the implementation of actuation forces on 3D structures
% with TET4 elements, using a ROM formulation
% 
% Last modified: 12/03/2023, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc

FORMULATION = 'N1'; % N1/N1t/N0
VOLUME = 1;         % integration over defected (1) or nominal volume (0)

USEJULIA = 0;

%% PREPARE MODEL __________________________________________________________                                                  

% DATA ____________________________________________________________________
E       = 0.6*263824;     % Young's modulus [Pa]
rho     = 1070;     % density [kg/m^3]
nu      = 0.499;     % Poisson's ratio 
thickness = .1;     % [m] beam's out-of-plane thickness

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain
% Element
myElementConstructor = @()Tet4Element(myMaterial);

% MESH_____________________________________________________________________
filename = 'naca0012_TET4_350';%'testPartTET4D4';
[nodes, elements, nset, elset] = mesh_ABAQUSread(filename);

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

%% PLOT MESH WITH NODES AND ELEMENTS ______________________________________

elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMeshAxis(nodes, elementPlot, 0);
hold off

%% BOUNDARY CONDITIONS ____________________________________________________

nel = size(elements,1);
nset = {};
Lx = max(nodes(:,1)) - min(nodes(:,1));
Ly = max(nodes(:,2)) - min(nodes(:,2));

for el=1:nel   
    nodePosX = 0;
    for n=1:size(elements,2)
        if  nodes(elements(el,n),1)<Lx/5 && ~any(cat(2, nset{:}) == elements(el,n))
              nset{end+1} = elements(el,n);         
        end
    end   
end

for l=1:length(nset)
    MeshNominal.set_essential_boundary_condition([nset{l}],1:3,0)
end

%% ASSEMBLY _______________________________________________________________

% nominal
NominalAssembly = Assembly(MeshNominal);
Mn = NominalAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( MeshNominal.nDOFs, 1);
[Kn,~] = NominalAssembly.tangent_stiffness_and_force(u0);
    % store matrices
    NominalAssembly.DATA.K = Kn;
    NominalAssembly.DATA.M = Mn;

%% DAMPING ________________________________________________________________

alfa = 3.1;
beta = 6.3*1e-6;

Dn = alfa*Mn + beta*Kn; % Rayleigh damping 
NominalAssembly.DATA.D = Dn;
Dc = NominalAssembly.constrain_matrix(Dn);

%% EIGENMODES - VIBRATION MODES (VMs) _____________________________________

n_VMs = 3;

% eigenvalue problem
Kc = NominalAssembly.constrain_matrix(Kn);
Mc = NominalAssembly.constrain_matrix(Mn);
[VMn,om] = eigs(Kc, Mc, n_VMs, 'SM');
[f0n,ind] = sort(sqrt(diag(om))/2/pi)
VMn = VMn(:,ind)
for ii = 1:n_VMs
    VMn(:,ii) = VMn(:,ii)/max(sqrt(sum(VMn(:,ii).^2,2)));
end
VMn = NominalAssembly.unconstrain_vector(VMn);

% plot
mod = 1;
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes, elementPlot, 0);
v1 = reshape(VMn(:,mod), 3, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', max(nodes(:,2)));
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0n(mod),3) ' Hz']);

%% MODAL DERIVATIVES (MDs) ________________________________________________                    

[MDn, MDname] = modal_derivatives(NominalAssembly, elements, VMn);

%% ROM TENSORS (INTERNAL FORCES) __________________________________________                      

% define reduced order basis
Vn = [VMn MDn];     % reduced order basis (ROM-n)

% orthonormalize reduction basis
Vn = orth(Vn);	% ROM-n

% ROM-n: standard reduced order model (no defects in the mesh)
tensors_ROMn = reduced_tensors_ROM(NominalAssembly, elements, Vn, USEJULIA);

%% REDUCED ASSEMBLIES _____________________________________________________

% ROM-n ___________________________________________________________________
ROMn_Assembly = ReducedAssembly(MeshNominal, Vn);
ROMn_Assembly.DATA.M = ROMn_Assembly.mass_matrix();  % reduced mass matrix (ROM-n)
ROMn_Assembly.DATA.C = Vn.'*Dn*Vn;    % reduced damping matrix (ROM-n), using C as needed by the residual function of Newmark integration
ROMn_Assembly.DATA.K = Vn.'*Kn*Vn;    % reduced stiffness matrix (ROM-n)

%% ROM TENSORS - ACTUATION FORCES _________________________________________

[skin,allfaces,skinElements, skinElementFaces] = getSkin3D(elements);

nel = size(elements,1);
actuationElements = zeros(nel,1);
for el=1:nel
    for n=1:size(elements,2)
        if nodes(elements(el,n),1)>Lx/5 && nodes(elements(el,n),1)<Lx*0.5 && nodes(elements(el,n),2)>0.2*Ly
            actuationElements(el) = 1;   
        end
    end
end

actuationDirection = [0;1;0;0;0;0]; % [0;1;0]-->[x^2,y^2,z^2,xy,xz,yz]

% ROM-n
tensors_actuation_ROMn = reduced_tensors_actuation_ROM(NominalAssembly, Vn, actuationElements, actuationDirection);

%% PLOT ACTUATION ELEMENT IN RED __________________________________________

disp=zeros(size(nodes,1),3);
PlotFieldonDeformedMeshActuation(nodes,elements,actuationElements,0.8,disp,'factor',1,'color', 'k','cameraPos',[-2,2,5]) ;

%% TIME INTEGRATION _______________________________________________________

% ROM-n __________________________________________________________________
% time step for integration
h = 0.05;
% initial condition: equilibrium
q0 = zeros(size(Vn,2),1);
qd0 = zeros(size(Vn,2),1);
qdd0 = zeros(size(Vn,2),1);

% Actuation forces
B1 = tensors_actuation_ROMn.B1;
B2 = tensors_actuation_ROMn.B2;
F_ext = @(t,q) 1/2*(1-(1+0.01*sin(t*2*pi/10)))*(B1+B2*q); % q is the reduced order DOFs

% instantiate object for nonlinear time integration
TI_NL_ROMn = ImplicitNewmark('timestep',h,'alpha',0.005);

% modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_actuation(q,qd,qdd,t,ROMn_Assembly,F_ext);

% nonlinear Time Integration
tmax = 10; 
TI_NL_ROMn.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_ROMn.Solution.u = Vn * TI_NL_ROMn.Solution.q; % get full order solution

%% ANIMATION ______________________________________________________________

actuationValues = zeros(size(TI_NL_ROMn.Solution.u,2),1);
for t=1:size(TI_NL_ROMn.Solution.u,2)
    actuationValues(t) = 1+0.01*sin(t*h*2*pi/10);
end

AnimateFieldonDeformedMeshActuation(nodes, elementPlot,actuationElements,actuationValues,TI_NL_ROMn.Solution.u,'factor',1,'index',1:3,'cameraPos',[-2,3,5],'filename','result_video','framerate',1/h)

%% PLOT ACTUATION SPECIFIC FRAME __________________________________________

frameNumber = 7.5/h;
U = reshape(TI_NL_ROMn.Solution.u(:,frameNumber),3,[]).';
disp = U(:,1:3);
normalizedActuationValues = normalize(actuationValues,'range',[-1 1]);
PlotFieldonDeformedMeshActuation(nodes,elements,actuationElements,normalizedActuationValues(frameNumber),disp,'factor',1,'color', 'k','cameraPos',[-2,3,5]) ;


