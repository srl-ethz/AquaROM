% ------------------------------------------------------------------------ 
% Playground to test fish locomotion forces using a FOM formulation
% Used element type: TRI3.
% 
% Last modified: 09/10/2023, Mathieu Dubied, ETH Zurich
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

% nodes(:,1) = -nodes(:,1);

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
        % if  nodes(elements(el,n),1)>-Lx*0.1 && ~any(cat(2, nset{:}) == elements(el,n))
        if  nodes(elements(el,n),1)>-Lx*0.1 && ~any(cat(2, nset{:}) == elements(el,n))
            nset{end+1} = elements(el,n);
        end
    end   
end
for l=1:length(nset)
    MeshNominal.set_essential_boundary_condition([nset{l}],2,0);
end

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

%% DAMPING ________________________________________________________________
alfa = 0.912;
beta = 0.002;

% nominal
Dn = alfa*Mn + beta*Kn; % Rayleigh damping 
NominalAssembly.DATA.D = Dn;
NominalAssembly.DATA.C = Dn;
Dc = NominalAssembly.constrain_matrix(Dn);

%% VISUALISATION __________________________________________________________
% plot nominal mesh
elementPlot = elements(:,1:3); 
figure('units','normalized','position',[.2 .1 .6 .4],'name','Nominal mesh')
PlotMesh(nodes, elementPlot, 1);

%% SPINE AND TAIL OPERATIONS ______________________________________________
% get nodes of the spine
[spineNodes, spineElements, spineElementWeights, nodeIdxPosInElements] = find_spine_TRI3(elements,nodes);

% get nodes of the tail
[tailNode, tailElement, tailElementWeights] = find_tail_TRI3(elements,nodes,spineElements,nodeIdxPosInElements);

tailElementWeights = sparse(tailElementWeights);

% get normalisation factors, i.e., the inverse of the distance between each consecutive spine node
normalisationFactors = compute_normalisation_factors(nodes, elements, spineElements, nodeIdxPosInElements);

% compute force for a hypothetical position u and ud
q = NominalAssembly.constrain_vector(reshape(nodes.',[],1));
qd = zeros(size(q));
qd(2:2:end) = 4;
qdd = zeros(size(q));
mTilde = 1000;

% compute force
% fReactive = reactive_force_TRI3(NominalAssembly, spineElementWeights, nodeIdxPosInElements,normalisationFactors, mTilde, q, qd);
fTailPressure = tail_pressure_force_TRI3(NominalAssembly, tailElementWeights, nodeIdxPosInElements,normalisationFactors, mTilde, q, qd);
% disp(NominalAssembly.unconstrain_vector(fTailPressure));
disp(fTailPressure)

%% SPINE CHANGE OF MOMENTUM TENSOR ________________________________________
T = spine_momentum_tensor_TRI(NominalAssembly, spineElementWeights,nodeIdxPosInElements,normalisationFactors,mTilde);

q = reshape(nodes.',[],1);
qd = zeros(size(q));
qd(2:2:end) = 4;
qd(1:2:end) = 0;
qdd = zeros(size(q));
qdd(2:2:end) = -4;
mTilde = 3;
spineForceDirect = spine_momentum_force_TRI(NominalAssembly, spineElementWeights,nodeIdxPosInElements,normalisationFactors,mTilde,q,qd,qdd)
spineForce = ttv(ttv(ttv(T,qdd,2),q,2),q,2)+ttv(ttv(ttv(T,qd,2),qd,2),q,2)+ttv(ttv(ttv(T,qd,2),q,2),qd,2)



%% PLOT FORCE _____________________________________________________________
% PlotMeshandForce(nodes,elements,0,NominalAssembly.unconstrain_vector(fTailPressure)/30)
PlotMeshandForce(nodes,elements,0,fTailPressure/30)
%% TAIL OPERATIONS ________________________________________________________
% get tail elements and tail nodes
[tailElement, tailNodeIndexInElement, tailElementBinVec] = find_tail_TRI3(elements, nodes);

%% ROM TENSORS - HYDRODYNAMIC FORCES ______________________________________
% [skin,allfaces,skinElements, skinElementFaces] = getSkin2D(elements);
% vwater = [0.5;0.5];   % water velocity vector
% 
% rho = 997*0.01;
% c = 0.2;
% 
% fTest = hydro_force_TRI3(NominalAssembly, skinElements, skinElementFaces, vwater, rho, c, u, ud)

%% FOM TENSORS - HYDRODYNAMIC FORCES (optional) ___________________________

FOM = 1;    % FOM=1 for assembling the hydrodynamic tensors at the assembly level (FOM)

% if FOM == 1
%     tensors_hydro_FOM = unreduced_tensors_hydro_FOM(NominalAssembly, elements, skinElements, skinElementFaces, vwater, rho, c);
% end
% 
% % hydrodynamic forces
% T1 = NominalAssembly.constrain_vector(double(tensors_hydro_FOM.T1));
% Tu2 = NominalAssembly.constrain_matrix(double(tensors_hydro_FOM.Tu2));
% Tudot2 = NominalAssembly.constrain_matrix(double(tensors_hydro_FOM.Tudot2));
% Tuu3 = tensor(NominalAssembly.constrain_tensor(double(tensors_hydro_FOM.Tuu3)));
% Tuudot3 = tensor(NominalAssembly.constrain_tensor(double(tensors_hydro_FOM.Tuudot3)));
% Tudotudot3 = tensor(NominalAssembly.constrain_tensor(double(tensors_hydro_FOM.Tudotudot3)));
% 
% toc

%% ACTUATION ELEMENTS _____________________________________________________
% get skin elements
% [skin,allfaces,skinElements, skinElementFaces] = getSkin2D(elements);

% compute actuation tensors
tensors_actu = create_actuation_tensors(NominalAssembly, elements, nodes);

%% TIME INTEGRATION _______________________________________________________

% FOM-I Sn (original nonlinear force, small time step needed) _____________
% initial condition: equilibrium
fprintf('solver \n')
tic
h2=0.01;
tmax = 3; 
nUncDOFs = size(MeshNominal.EBC.unconstrainedDOFs,2);
q0 = zeros(nUncDOFs,1);
qd0 = zeros(nUncDOFs,1);
qd0(2:2:end) = 0.1;
qdd0 = zeros(nUncDOFs,1);

% external forces: tail pressure force and actuation
% k=0.01;
% B1TopMuscle = tensors_actu.B1_top;
% B2TopMuscle = tensors_actu.B2_top;
% B1BottomMuscle =tensors_actu.B1_bottom;
% B2BottomMuscle = tensors_actu.B2_bottom;

forceTest = zeros(size(reshape(nodes.',[],1)));
forceTest(1:2:end)=0.5;

F_ext = @(t,q,qd) forceTest;
%     + k/2*(1-(1+0.04*sin(t*2*pi)))*(B1TopMuscle+B2TopMuscle*q) + ...
%             k/2*(1-(1-0.04*sin(t*2*pi)))*(B1BottomMuscle+B2BottomMuscle*q);

% F_ext = @(t,q,qd) k/2*(1-(1+0.06*1*t/10))*(B1TopMuscle+B2TopMuscle*NominalAssembly.unconstrain_vector(q)) ;

% F_ext = @(t,q,qd) k/2*(1-(1+0.06*sin(t*2*pi)))*(B1TopMuscle+B2TopMuscle*NominalAssembly.unconstrain_vector(q)) + ...
%             k/2*(1-(1-0.06*sin(t*2*pi)))*(B1BottomMuscle+B2BottomMuscle*NominalAssembly.unconstrain_vector(q));
% instantiate object for nonlinear time integration
TI_NL_FOMfull = ImplicitNewmark('timestep',h2,'alpha',0.005,'MaxNRit',400,'MaxNRit',200,'RelTol',1e-6);

% modal nonlinear Residual evaluation function handle
Residual_NL = @(q,qd,qdd,t)residual_nonlinear(q,qd,qdd,t,NominalAssembly,F_ext);

% nonlinear Time Integration
TI_NL_FOMfull.Integrate(q0,qd0,qdd0,tmax,Residual_NL);
TI_NL_FOMfull.Solution.u = zeros(NominalAssembly.Mesh.nDOFs,size(TI_NL_FOMfull.Solution.q,2));

for t=1:size(TI_NL_FOMfull.Solution.q,2)
    TI_NL_FOMfull.Solution.u(:,t) = NominalAssembly.unconstrain_vector(TI_NL_FOMfull.Solution.q(:,t));
end
toc

%% PLOT
tailNode=1
plot(TI_NL_FOMfull.Solution.qd(tailNode*2-1,:))

%% CHECK TAIL PRESSURE FORCE ___________________________________________
timeStep = 40;
qTest = TI_NL_FOMfull.Solution.q(:,timeStep);
qdTest = TI_NL_FOMfull.Solution.qd(:,timeStep);

%%
elementPlot = elements(:,1:3); % plot only corners (otherwise it's a mess)
AnimateFieldonDeformedMesh(nodes, elementPlot,TI_NL_FOMfull.Solution.u, ...
    'factor',1,'index',1:2,'filename','result_video','framerate',1/h2)

%%
fTailPressure = tail_pressure_force_TRI3(NominalAssembly, tailElementWeights, nodeIdxPosInElements,normalisationFactors, mTilde, qTest, qdTest);
disp(NominalAssembly.unconstrain_vector(fTailPressure));
