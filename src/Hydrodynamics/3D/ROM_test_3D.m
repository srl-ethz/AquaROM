% ------------------------------------------------------------------------ 
% Script to test the implementation of hydrodynamic forces on 3D structures
% with TET4 elements, using a PROM formulation
% 
% Last modified: 13/03/2023, Mathieu Dubied, ETH Zurich
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
filename = 'testPartTET4D4';%'fish1_728el_rot';%'testPartTET4D4';
[nodes, elements, nset, elset] = mesh_ABAQUSread(filename);
% nodes = nodes*0.01;

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
[f0n,ind] = sort(sqrt(diag(om))/2/pi);
VMn = VMn(:,ind);
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

%% ROM TENSORS - HYDRODYNAMIC FORCES ______________________________________
[skin,allfaces,skinElements, skinElementFaces] = getSkin3D(elements);
vwater = [3;6;0.0];
rho = 1;

% ROM-n
tensors_hydro_ROMn = reduced_tensors_hydro_ROM(NominalAssembly, elements, Vn, skinElements, skinElementFaces, vwater, rho);
% ROM-d
%tensors_hydro_ROMd = reduced_tensors_hydro_ROM(DefectedAssembly, elements, Vd, skinElements, skinElementFaces, vwater, rho);
% PROM
%tensors_hydro_PROM = reduced_tensors_hydro_PROM(NominalAssembly, elements, V, U, FOURTHORDER, skinElements, skinElementFaces, vwater, rho);

%% TIME INTEGRATION

% Parameters' initialization for all models _______________________________
% time step for integration
h = 0.05;

%% ROM-n __________________________________________________________________
% initial condition: equilibrium
q0 = zeros(size(Vn,2),1);
qd0 = zeros(size(Vn,2),1);
qdd0 = zeros(size(Vn,2),1);

% hydrodynamic forces
F_ext = @(t,q,qd) (double(tensors_hydro_ROMn.Tr1) + ...
     double(tensors_hydro_ROMn.Tru2*q) + double(tensors_hydro_ROMn.Trudot2*qd)); %+ ...
%     double(ttv(ttv(tensors_hydro_ROMn.Truu3,q,3), q,2)) + ...
%     double(ttv(ttv(tensors_hydro_ROMn.Truudot3,qd,3), q,2)) + ...
%     double(ttv(ttv(tensors_hydro_ROMn.Trudotudot3,qd,3), qd,2))); % q, qd are reduced order DOFs

% instantiate object for nonlinear time integration
TI_NL_ROMn = ImplicitNewmark('timestep',h,'alpha',0.005);

% modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hydro(q,qd,qdd,t,ROMn_Assembly,F_ext);

% nonlinear Time Integration
tmax = 4.0; 
TI_NL_ROMn.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_ROMn.Solution.u = Vn * TI_NL_ROMn.Solution.q; % get full order solution

%% ANIMATION
actuationValues = zeros(size(TI_NL_ROMn.Solution.u,2),1);
actuationElements = zeros(size(elements,1));

AnimateFieldonDeformedMeshActuation(nodes, elementPlot,actuationElements,actuationValues,TI_NL_ROMn.Solution.u,'factor',1,'index',1:3,'cameraPos',[-2,3,5],'filename','result_video','framerate',1/h)

















%% PREPARE MODEL                                                    

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

% PROM parameters
xi1 = 0.2;

% MESH_____________________________________________________________________
filename = 'naca0012_TET4_350';
[nodes, elements, ~, elset] = mesh_ABAQUSread(filename);

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);




%%
% boundary conditions of nominal mesh: front and back nodes fixed
frontNode = find_node_2D(0,0,nodes);
Lx = 0.15;
backNode = find_node_2D(Lx,0,nodes);
frontNode2 = find_node_2D(0.005,0,nodes);
%nset = {frontNode,backNode};
nset = {frontNode, frontNode2};
MeshNominal.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)

% defect shapes
% (1) thinning airfoil 
Lx=0.15;
nodes_translated = [nodes(:,1), nodes(:,2)*0];
vdA = nodes_translated(:,2) - nodes(:,2);
thinAirfoil = zeros(numel(nodes),1);
thinAirfoil(2:2:end) = vdA;

% defected mesh
U = thinAirfoil;   % defect basis
xi = xi1;           % parameter vector
m = length(xi);     % number of parameters

% update defected mesh nodes
d = U*xi;                       % displacement fields introduced by defects
dd = [d(1:2:end) d(2:2:end)]; 
nodes_defected = nodes + dd;    % nominal + d ---> defected 
DefectedMesh = Mesh(nodes_defected);
DefectedMesh.create_elements_table(elements,myElementConstructor);
DefectedMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)

figure('units','normalized','position',[.2 .3 .6 .4])
elementPlot = elements(:,1:3); hold on % plot only corners (otherwise it's a mess)
PlotMesh(nodes_defected, elementPlot, 0); 
PlotMesh(nodes,elementPlot,0);
v1 = reshape(U*xi, 2, []).';
S = 1;
hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
title(sprintf('Defect, \\xi=[%.1f], S=%.1f\\times',...
    xi1, S))
axis equal; grid on; box on; set(hf{1},'FaceAlpha',.7); drawnow

%% ASSEMBLY ________________________________________________________________

% nominal
NominalAssembly = Assembly(MeshNominal);
Mn = NominalAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( MeshNominal.nDOFs, 1);
[Kn,~] = NominalAssembly.tangent_stiffness_and_force(u0);
    % store matrices
    NominalAssembly.DATA.K = Kn;
    NominalAssembly.DATA.M = Mn;

 %%
[skin,allfaces,skinElements,skinElementFaces] = getSkin3D(elements);
F = NominalAssembly.skin_force('force_length_prop_skin_normal', 'weights', skinElements, skinElementFaces);
PlotMeshandForce(nodes,elements,0,F);
%%
% defected
DefectedAssembly = Assembly(DefectedMesh);
Md = DefectedAssembly.mass_matrix();
[Kd,~] = DefectedAssembly.tangent_stiffness_and_force(u0);
    % store matrices
    DefectedAssembly.DATA.K = Kd;
    DefectedAssembly.DATA.M = Md;


%% PLOT MESH WITH NODES AND ELEMENTS
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes, elementPlot, 1);

%% DAMPING 
alfa = 3.1;
beta = 6.3*1e-6;
% nominal
Dn = alfa*Mn + beta*Kn; % Rayleigh damping 
NominalAssembly.DATA.D = Dn;
Dc = NominalAssembly.constrain_matrix(Dn);
% defected
Dd = alfa*Md + beta*Kd; % Rayleigh damping 
DefectedAssembly.DATA.D = Dd;
Ddc = DefectedAssembly.constrain_matrix(Dd);

%% EIGENMODES - VIBRATION MODES (VMs)

n_VMs = 2;

% NOMINAL _________________________________________________________________
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

% plot
mod = 1;
elementPlot = elements(:,1:3); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes, elementPlot, 0);
v1 = reshape(VMn(:,mod), 2, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', max(nodes(:,2)));
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0n(mod),3) ' Hz']);

% DEFECTED ________________________________________________________________
% eigentvalue problem
Kdc = DefectedAssembly.constrain_matrix(Kd);
Mdc = DefectedAssembly.constrain_matrix(Md);
[VMd,om] = eigs(Kdc, Mdc, n_VMs, 'SM');
[f0d,ind] = sort(sqrt(diag(om))/2/pi);
VMd = VMd(:,ind);
for ii = 1:n_VMs
    VMd(:,ii) = VMd(:,ii)/max(sqrt(sum(VMd(:,ii).^2,2)));
end
VMd = DefectedAssembly.unconstrain_vector(VMd);

% plot
mod = 1;
elementPlot = elements(:,1:3); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes, elementPlot, 0);
v1 = reshape(VMd(:,mod), 2, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', max(nodes(:,2)));
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0d(mod),3) ' Hz'])

%% MODAL DERIVATIVES (MDs)                     

% nominal
[MDn, MDname] = modal_derivatives(NominalAssembly, elements, VMn);
% defected
MDd = modal_derivatives(DefectedAssembly, elements, VMd);
 
% defect sensitivities
[DS, names] = defect_sensitivities(NominalAssembly, elements, VMn, U, ...
    FORMULATION);

%% ROM TENSORS (INTERNAL FORCES)                         
% define reduced order basis
Vn = [VMn MDn];     % reduced order basis (ROM-n)
Vd = [VMd MDd];   	% reduced order basis (ROM-d)
V  = [VMn MDn DS]; 	% reduced order basis (PROM) 

% orthonormalize reduction basis
Vn = orth(Vn);	% ROM-n
Vd = orth(Vd);	% ROM-d
V  = orth(V);	% PROM

% ROM-n: standard reduced order model (no defects in the mesh)
tensors_ROMn = reduced_tensors_ROM(NominalAssembly, elements, Vn, USEJULIA);

% ROM-d: standard reduced order model (defects in the mesh)
tensors_ROMd = reduced_tensors_ROM(DefectedAssembly, elements, Vd, USEJULIA);
tensors_ROMd.xi = xi; % save for which xi ROMd is computed

% PROM: parametric formulation for defects
tensors_PROM = reduced_tensors_DpROM(NominalAssembly, elements, ...
    V, U, FORMULATION, VOLUME, USEJULIA); %compute tensors

% evaluate the defected tensors at xi
[Q2, Q3, Q4, Q3t, Q4t, M] = DefectedTensors(tensors_PROM, xi);

%% REDUCED ASSEMBLIES

% ROM-n ___________________________________________________________________
ROMn_Assembly = ReducedAssembly(MeshNominal, Vn);
ROMn_Assembly.DATA.M = ROMn_Assembly.mass_matrix();  % reduced mass matrix (ROM-n)
ROMn_Assembly.DATA.C = Vn.'*Dn*Vn;    % reduced damping matrix (ROM-n), using C as needed by the residual function of Newmark integration
ROMn_Assembly.DATA.K = Vn.'*Kn*Vn;    % reduced stiffness matrix (ROM-n)

% ROM-d ___________________________________________________________________
ROMd_Assembly = ReducedAssembly(DefectedMesh, Vd);
ROMd_Assembly.DATA.M = ROMd_Assembly.mass_matrix();  % reduced mass matrix (ROM-d)
ROMd_Assembly.DATA.C = Vd.'*Dd*Vd;    % reduced damping matrix (ROM-d), using C as needed by the residual function of Newmark integration
ROMd_Assembly.DATA.K = Vd.'*Kd*Vd;    % reduced stiffness matrix (ROM-d)

% PROM ____________________________________________________________________
PROM_Assembly = ReducedAssembly(MeshNominal, V);
PROM_Assembly.DATA.M = PROM_Assembly.mass_matrix();  % reduced mass matrix (PROM)
PROM_Assembly.DATA.C = V.'*Dn*V;    % reduced damping matrix (PROM), using C as needed by the residual function of Newmark integration
PROM_Assembly.DATA.K = V.'*Kn*V;    % reduced stiffness matrix (PROM)


%% ROM TENSORS - ACTUATION FORCES
[skin,allfaces,skinElements, skinElementFaces] = getSkin2D(elements);
actuationDirection = [1;0];

% ROM-n
tensors_hydro_ROMn = reduced_tensors_hydro_ROM(NominalAssembly, Vn, skinElements, actuationDirection);
% ROM-d
%tensors_hydro_ROMd = reduced_tensors_hydro_ROM(DefectedAssembly, elements, Vd, skinElements, skinElementFaces, vwater, rho);
% PROM
%tensors_hydro_PROM = reduced_tensors_hydro_PROM(NominalAssembly, elements, V, U, FOURTHORDER, skinElements, skinElementFaces, vwater, rho);

%% TIME INTEGRATION

% Parameters' initialization for all models _______________________________
% time step for integration
h = 0.05;

%% ROM-n __________________________________________________________________
% initial condition: equilibrium
q0 = zeros(size(Vn,2),1);
qd0 = zeros(size(Vn,2),1);
qdd0 = zeros(size(Vn,2),1);

% hydrodynamic forces
F_ext = @(t,q,qd) (double(tensors_hydro_ROMn.Tr1) + ...
    double(tensors_hydro_ROMn.Tru2*q) + double(tensors_hydro_ROMn.Trudot2*qd) + ...
    double(ttv(ttv(tensors_hydro_ROMn.Truu3,q,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_ROMn.Truudot3,qd,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_ROMn.Trudotudot3,qd,3), qd,2))); % q, qd are reduced order DOFs

% instantiate object for nonlinear time integration
TI_NL_ROMn = ImplicitNewmark('timestep',h,'alpha',0.005);

% modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hydro(q,qd,qdd,t,ROMn_Assembly,F_ext);

% nonlinear Time Integration
tmax = 4.0; 
TI_NL_ROMn.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_ROMn.Solution.u = Vn * TI_NL_ROMn.Solution.q; % get full order solution

%% ROM-d __________________________________________________________________
% initial condition: equilibrium
q0 = zeros(size(Vd,2),1);
qd0 = zeros(size(Vd,2),1);
qdd0 = zeros(size(Vd,2),1);

% hydrodynamic forces
F_ext = @(t,q,qd) (double(tensors_hydro_ROMd.Tr1) + ...
    double(tensors_hydro_ROMd.Tru2*q) + double(tensors_hydro_ROMd.Trudot2*qd) + ...
    double(ttv(ttv(tensors_hydro_ROMd.Truu3,q,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_ROMd.Truudot3,qd,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_ROMd.Trudotudot3,qd,3), qd,2))); % q, qd are reduced order DOFs

% instantiate object for nonlinear time integration (TI)
TI_NL_ROMd = ImplicitNewmark('timestep',h,'alpha',0.005);

% nonlinear residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hydro(q,qd,qdd,t,ROMd_Assembly,F_ext);

% nonlinear time integration (TI)
TI_NL_ROMd.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_ROMd.Solution.u = Vd * TI_NL_ROMd.Solution.q; % get full order solution

%% PROM - nominal (defect=0) ______________________________________________
% initial condition: equilibrium
q0 = zeros(size(V,2),1);
qd0 = zeros(size(V,2),1);
qdd0 = zeros(size(V,2),1);

% hydrodynamic forces
F_ext = @(t,q,qd) (double(tensors_hydro_PROM.Tr1) + ...
    double(tensors_hydro_PROM.Tru2*q) + double(tensors_hydro_PROM.Trudot2*qd) + ...
    double(ttv(ttv(tensors_hydro_PROM.Truu3,q,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_PROM.Truudot3,qd,3), q,2)) + ...
    double(ttv(ttv(tensors_hydro_PROM.Trudotudot3,qd,3), qd,2))); % q, qd are reduced order DOFs

% instantiate object for nonlinear time integration
TI_NL_PROM = ImplicitNewmark('timestep',h,'alpha',0.005);

% modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hydro(q,qd,qdd,t,PROM_Assembly,F_ext);

% nonlinear Time Integration
tmax = 4.0; 
TI_NL_PROM.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_PROM.Solution.u = V * TI_NL_PROM.Solution.q; % get full order solution

%% PROM-d (defect=xi) _____________________________________________________
% % initial condition: equilibrium
% q0 = zeros(size(V,2),1);
% qd0 = zeros(size(V,2),1);
% qdd0 = zeros(size(V,2),1);
% 
% % hydrodynamic forces
% F_ext = @(t,q,qd) (double(tensors_hydro_PROM.Tr1) + double(tensors_hydro_PROM.Tr2)*xi + ...
%     double(tensors_hydro_PROM.Tru2*q) + double(ttv(ttv(tensors_hydro_PROM.Tru3,xi,3), q,2)) + ...
%     double(tensors_hydro_PROM.Trudot2*qd) + double(ttv(ttv(tensors_hydro_PROM.Trudot3,xi,3), qd,2)) + ...
%     double(ttv(ttv(tensors_hydro_PROM.Truu3,q,3), q,2)) + ...
%     double(ttv(ttv(tensors_hydro_PROM.Truudot3,qd,3), q,2)) + ...
%     double(ttv(ttv(tensors_hydro_PROM.Trudotudot3,qd,3), qd,2))); % q, qd are reduced order DOFs
% 
% % instantiate object for nonlinear time integration
% TI_NL_PROMd = ImplicitNewmark('timestep',h,'alpha',0.005);
% 
% % modal nonlinear Residual evaluation function handle
% Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hydro(q,qd,qdd,t,PROM_Assembly,F_ext);
% 
% % nonlinear Time Integration
% tmax = 4.0; 
% TI_NL_PROMd.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
% TI_NL_PROMd.Solution.u = V * TI_NL_PROMd.Solution.q; % get full order solution

%% SENSITIVITY ANALYSIS ___________________________________________________
% initial conditions
s0 = zeros(size(V,2),1);
sd0 = zeros(size(V,2),1);
sdd0 = zeros(size(V,2),1);

% evaluate partial derivatives along solution
pd_fext_PROM = @(q,qd)DpROM_hydro_derivatives(q,qd,xi,tensors_hydro_PROM,FOURTHORDER);
pd_fint_PROM = @(q)DpROM_derivatives(q,tensors_PROM); 

% instantiate object for time integration
TI_sens = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true,'sens',true);

% residual function handle
Residual_sens = @(s,sd,sdd,t,it)residual_linear_sens(s,sd,sdd,t,PROM_Assembly,TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd, TI_NL_PROM.Solution.qdd,pd_fext_PROM,pd_fint_PROM,h);

% time integration (TI)
TI_sens.Integrate(q0,qd0,qdd0,tmax,Residual_sens);
TI_sens.Solution.s = V * TI_sens.Solution.q; % get full order solution


%% VISUALIZE ______________________________________________________________

%% 1-DOF PLOT _____________________________________________________________

% find a specific result node and corresponding DOF
rNodeDOFS = MeshNominal.get_DOF_from_location([0.9*Lx, 0]);
rNodeDOF = rNodeDOFS(2); % y-direction

% plot
figure 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
tplot=linspace(0,tmax,tmax/h);
plot(tplot,TI_NL_ROMn.Solution.u(rNodeDOF,1:end-2)*100)
hold on
plot(tplot,TI_NL_ROMd.Solution.u(rNodeDOF,1:end-2)*100, "-.")
plot(tplot,TI_NL_PROM.Solution.u(rNodeDOF,1:end-2)*100, "--")
%plot(tplot,TI_NL_PROMd.Solution.u(rNodeDOF,1:end-2)*100, ":")
%plot(tplot,(TI_NL_PROM.Solution.u(rNodeDOF,1:end-2)+TI_sens.Solution.s(rNodeDOF,1:end-2)*xi)*100,":")
approx = V*(TI_NL_PROM.Solution.q(:,1:end-2)+TI_sens.Solution.q(:,1:end-2)*xi)*100;
plot(tplot,approx(rNodeDOF,:))

title('Vertical displacement for a specific node')
ylabel('$$u_y \mbox{ [cm]}$$','Interpreter','latex')
xlabel('Time [s]')
set(gca,'FontName','ComputerModern');
grid on
%legend({'ROM-n','ROM-d','PROM-n','PROM-d','Approximation of ROM-d using PROM sensitivities'}, 'Location', 'southoutside','Orientation','horizontal')
legend({'ROM-n','ROM-d','PROM-n','Approximation of ROM-d using PROM sensitivities'}, 'Location', 'southoutside','Orientation','horizontal')
hold off


%% Animations _____________________________________________________________
%% ROM-n __________________________________________________________________
AnimateFieldonDeformedMesh(nodes, elementPlot,TI_NL_ROMn.Solution.u,'factor',1,'index',1:2,'filename','result_video')

%% ROM-d __________________________________________________________________
AnimateFieldonDeformedMesh(nodes, elementPlot,TI_NL_ROMd.Solution.u,'factor',100,'index',1:2,'filename','result_video')


%% OPTIMIZATION ___________________________________________________________
d = [1;0];
%pd_fhydro_PROM = @(q,qd,xi)DpROM_hydro_derivatives(q,qd,xi,tensors_hydro_PROM);
[xiStar,LrEvo] = optimization_pipeline_1(V,d,tensors_PROM, tensors_hydro_PROM,TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd,TI_sens.Solution.q,TI_sens.Solution.qd)
%xi_star = optimization_pipeline_1(Vd,d,tensors_ROMd, tensors_hydro_ROMd,TI_NL_ROMd.Solution.q,TI_NL_ROMd.Solution.qd)

