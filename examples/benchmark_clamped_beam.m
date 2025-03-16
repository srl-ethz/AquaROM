% ------------------------------------------------------------------------ 
% Accuracy comparison: clamped beam
% 
% Last modified: 12/02/2024, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc
percentageError = [];
%% GET REAL DATA __________________________________________________________

measfilename = 'gravity_left_z_real_data.dat'; 

% copy the file to a legible format 
measfilenamefilename2 = [measfilename(1:end-3), 'txt']; 
copyfile(measfilename, measfilenamefilename2);

% read table of legible file type
solReal = readtable(measfilenamefilename2, 'ReadVariableNames', 0);
solReal = table2array(solReal);

% delete the file that we created bc we don't need it anymore
delete(measfilenamefilename2);

%% PREPARE MODEL __________________________________________________________                                                   

% DATA ____________________________________________________________________
E       = 263824;       % Young's modulus [Pa] - from sim-to-real paper
rho     = 1070;         % density [kg/m^3]
nu      = 0.499;        % Poisson's ratio 

% material
elementType = 'TET4';

myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	    % set "false" for plane_strain
switch elementType
    case 'HEX8'
        myElementConstructor = @()Hex8Element(myMaterial);
        nID = 8;
    case 'HEX20'
        myElementConstructor = @()Hex20Element(myMaterial);
        nID = 20;
    case 'TET4'
        myElementConstructor = @()Tet4Element(myMaterial);
        nID = 4;
    case 'TET10'
        myElementConstructor = @()Tet10Element(myMaterial);
        nID = 10;
    case 'WED15'
        myElementConstructor = @()Wed15Element(myMaterial);
        nID = 15;
%         quadrature = struct('lin', 5, 'tri', 12);
%         myElementConstructor = @()Wed15Element(myMaterial, quadrature);
end


% MESH ____________________________________________________________________
l = 0.1;
w = .03;
t = .03;
nx = 15;
ny = 3;
nz = 3;
[nodes, elements, nset]=mesh_3Dparallelepiped(elementType,l,w,t,nx,ny,nz);

nel = size(elements,1);

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of airfoil
Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of airfoil
Lz = abs(max(nodes(:,3))-min(nodes(:,3)));  % vertical length of airfoil

% plot nominal mesh
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
%PlotMeshAxis(nodes, elementPlot, 0);
hold off

% boundary conditions: clamped end at x=0
nel = size(elements,1);
MeshNominal.set_essential_boundary_condition([nset{1}],1:3,0)

% ASSEMBLY (FOM) __________________________________________________________
BeamAssembly = Assembly(MeshNominal);
M = BeamAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( MeshNominal.nDOFs, 1);
[K,~] = BeamAssembly.tangent_stiffness_and_force(u0);

% store matrices
BeamAssembly.DATA.K = K;
BeamAssembly.DATA.M = M;

alfa = 0.912;
beta = 0.002;
D = alfa*M + beta*K; % Rayleigh damping 
BeamAssembly.DATA.D = D;
BeamAssembly.DATA.C = D;
Dc = BeamAssembly.constrain_matrix(D);

%% STATIC TEST ____________________________________________________________
% gravitational force
g = 9.81;
F = -rho*g*BeamAssembly.uniform_body_force();
[u_lin, u] = static_equilibrium(BeamAssembly, F, 'method', 'newton', ...
    'nsteps', 10, 'display', 'iter');
ULIN = reshape(u_lin,3,[]).';	% Linear response
UNL = reshape(u,3,[]).';        % Nonlinear response

maxDispSim = max(abs(UNL(:,3)))
maxDispReal = abs(solReal(end,2)-0.024)

percentageError = [percentageError;
    nID,(maxDispSim - maxDispReal)/maxDispReal*100];



%% SIMULATION PARAMETERS __________________________________________________
h = 0.01;
tmax = 1.63;

%% FOM SIMULATION _________________________________________________________
tStart = tic;
% gravitational force
g = 9.81;
Fext = @(t,q,qd) -rho*g*BeamAssembly.uniform_body_force();

% initial conditions
nUncDOFs = size(MeshNominal.EBC.unconstrainedDOFs,2);
q0 = zeros(nUncDOFs,1);
qd0 = zeros(nUncDOFs,1);
qdd0 = zeros(nUncDOFs,1);

% instantiate object for nonlinear time integration
TI_NL_FOM = ImplicitNewmark('timestep',h,'alpha',0.005,'MaxNRit',400,'MaxNRit',200,'RelTol',1e-6);
Residual_NL = @(q,qd,qdd,t)residual_nonlinear(q,qd,qdd,t,BeamAssembly,Fext);

% nonlinear Time Integration
TI_NL_FOM.Integrate(q0,qd0,qdd0,tmax,Residual_NL);
TI_NL_FOM.Solution.u = zeros(BeamAssembly.Mesh.nDOFs,size(TI_NL_FOM.Solution.q,2));

for t=1:size(TI_NL_FOM.Solution.q,2)
    TI_NL_FOM.Solution.u(:,t) = BeamAssembly.unconstrain_vector(TI_NL_FOM.Solution.q(:,t));
end
sol = TI_NL_FOM.Solution.u;
toc

fprintf('Time needed to solve the problem using FOM: %.2fsec\n',toc(tStart))
timeSim = toc(tStart);


%% PLOT ___________________________________________________________________

timePlotSim = linspace(0,tmax,size(sol,2));
idxNodeMarker = find_node(0.1,0,0.024,nodes);
ySolMarker = sol(idxNodeMarker*MeshNominal.nDOFPerNode,:);
posMarker = nodes(idxNodeMarker,3);
constDif = 0.024 - posMarker;
ySolMarkerPos = (posMarker + constDif)*ones(1,size(ySolMarker,2))+ySolMarker;

figure
plot(solReal(:,1),solReal(:,2),'DisplayName','Real data')
hold on
plot(timePlotSim, ySolMarkerPos,'DisplayName','Simulation')
xlabel('Time [s]')
ylabel('y-position motion marker')
legend('Location','southeast')
hold off

%% ROM ____________________________________________________________________
% EIGENMODES - VIBRATION MODES (VMs) ______________________________________
n_VMs = 4;  % number of vibration modes to include in the ROM

% eigenvalue problem
Kc = BeamAssembly.constrain_matrix(K);
Mc = BeamAssembly.constrain_matrix(M);
[VM,om] = eigs(Kc, Mc, n_VMs, 'SM');
[f0n,ind] = sort(sqrt(diag(om))/2/pi);
VM = VM(:,ind);
for ii = 1:n_VMs
    VM(:,ii) = VM(:,ii)/max(sqrt(sum(VM(:,ii).^2,2)));
end
VM = BeamAssembly.unconstrain_vector(VM);

% MODAL DERIVATIVES (MDs) _________________________________________________                   
[MD, MDname] = modal_derivatives(BeamAssembly, elements, VM);

% ROM TENSORS (INTERNAL FORCES) ___________________________________________                       
% define reduced order basis
V = [VM MD];     % reduced order basis (ROM)

% orthonormalize reduction basis
V = orth(V);	% ROM

% ROM-n: reduced order model for nominal mesh
toc
USEJULIA = 1;
tensors_ROMn = reduced_tensors_ROM(BeamAssembly, elements, V, USEJULIA);

% REDUCED ASSEMBLY ________________________________________________________
ROM_Assembly = ReducedAssembly(MeshNominal, V);
ROM_Assembly.DATA.M = ROM_Assembly.mass_matrix();     % reduced mass matrix 
ROM_Assembly.DATA.C = V.'*D*V;                      % reduced damping matrix 
ROM_Assembly.DATA.K = V.'*K*V;                      % reduced stiffness matrix 

%% TIME INTEGRATION _______________________________________________________

% initial condition: equilibrium
fprintf('solver')
tic
q0 = zeros(size(V,2),1);
qd0 = zeros(size(V,2),1);
qdd0 = zeros(size(V,2),1);

% instantiate object for nonlinear time integration
TI_NL_ROM = ImplicitNewmark('timestep',h,'alpha',0.005,'MaxNRit',60,'RelTol',1e-6);

% modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear(q,qd,qdd,t, ROM_Assembly, Fext)

% nonlinear Time Integration
TI_NL_ROM.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);

%%
TI_NL_ROM.Solution.u = V * TI_NL_ROM.Solution.q; % get full order solution
sol = TI_NL_ROM.Solution.u;
toc


%% PLOT ___________________________________________________________________

timePlotSim = linspace(0,tmax,size(sol,2));
idxNodeMarker = find_node(0.1,0,0.024,nodes);
ySolMarker = sol(idxNodeMarker*MeshNominal.nDOFPerNode,:);
posMarker = nodes(idxNodeMarker,3);
constDif = 0.024 - posMarker;
ySolMarkerPos = (posMarker + constDif)*ones(1,size(ySolMarker,2))+ySolMarker;

figure
plot(solReal(:,1),solReal(:,2),'DisplayName','Real data')
hold on
plot(timePlotSim, ySolMarkerPos,'DisplayName','Simulation')
xlabel('Time [s]')
ylabel('y-position motion marker')
legend('Location','southeast')
hold off








