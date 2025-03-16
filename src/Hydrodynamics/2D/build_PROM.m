% build_PROM
%
% Synthax:
% [V,ROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
%    build_PROM(MeshNominal,nodes,elements,mTilde,U,USEJULIA,VOLUME,FORMULATION)
%
% Description: Builds a PROM based on the nominal mesh and the shape
% variations U
%
% INPUTS: 
% (1) MeshNominal:  nominal mesh converted from Abaqus              
% (2) nodes:        nodes and their coordinates
% (3) elements:     elements described by corresponding nodes' index
% (4) mTilde:       virtual mass linear density of the fish + water
% (5) U:            shape variation matrix/vector
% (6) USEJULIA:     use of JULIA (1) for the computation of internal forces
%                   tensors - not tested, but was present as parameter in
%                   the DpROM branch
% (7) VOLUME:       integration over defected (1) or nominal volume (0)
% (8) FORMULATION:  N1/N1t/N0, Newman expansion of internal forces
%
% OUTPUTS:
% (1) V:                    ROB    
% (2) ROM_Assembly:         ROM assembly
% (3) tensors_ROM:          (reduced) tensors for the internal forces 
% (4) tailProperties:       properties of the tail pressure force
%                           (matrices, tail elements etc.)
% (5) spineProperties:      properties of the spine change in momentum
%                           (tensor, spine elements etc.)
% (6) dragProperties:       properties of the drag force
% (7) actuTop:              vectors and matrices related to the actuation
%                           muscle at the top
% (8) actuBottom:           vectors and matrices related to the actuation
%                           muscle at the bottom
%     
%
% Additional notes: -
%
% Last modified: 12/11/2023, Mathieu Dubied, ETH Zürich

function [V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
    build_PROM(MeshNominal,nodes,elements,mTilde,U,USEJULIA,VOLUME,FORMULATION)
    
    % 2D or 3D? ___________________________________________________________
    fishDim = size(nodes,2);

    % ASSEMBLY ____________________________________________________________
    NominalAssembly = Assembly(MeshNominal);
    Mn = NominalAssembly.mass_matrix();
    nNodes = size(nodes,1);
    u0 = zeros( MeshNominal.nDOFs, 1);
    [Kn,~] = NominalAssembly.tangent_stiffness_and_force(u0);
    % store matrices
    NominalAssembly.DATA.K = Kn;
    NominalAssembly.DATA.M = Mn;

    % DAMPING _____________________________________________________________ 
    alfa = 3.1;
    beta = 6.3*1e-6;
    Dn = alfa*Mn + beta*Kn; % Rayleigh damping 
    NominalAssembly.DATA.D = Dn;
    Dc = NominalAssembly.constrain_matrix(Dn);

    % ROB _________________________________________________________________
    
    % vibration modes
    n_VMs = 2;
    Kc = NominalAssembly.constrain_matrix(Kn);
    Mc = NominalAssembly.constrain_matrix(Mn);
    [VMn,om] = eigs(Kc, Mc, n_VMs, 'SM');
    [f0n,ind] = sort(sqrt(diag(om))/2/pi);
    VMn = VMn(:,ind);
    for ii = 1:n_VMs
        VMn(:,ii) = VMn(:,ii)/max(sqrt(sum(VMn(:,ii).^2,2)));
    end
    VMn = NominalAssembly.unconstrain_vector(VMn);

    % modal derivatives
    [MDn, MDname] = modal_derivatives(NominalAssembly, elements, VMn);

    % shape variation/defect sensitivities (PROM)
    [DS, names] = defect_sensitivities(NominalAssembly, elements, VMn, U, ...
    FORMULATION);
    
    % ROB formulation, case studies 1 and 2 (no rigid body mode)
    mSingle = [1 0]; % horizontal displacement
    m1 = repmat(mSingle,1,nNodes)';
    mSingle = [0 1];
    m2 = repmat(mSingle,1,nNodes)';
    V  = [m1 m2 VMn MDn DS];
    V  = orth(V);

    % % plot
    % mod = 2;
    % elementPlot = elements(:,1:3); % plot only corners (otherwise it's a mess)
    % figure('units','normalized','position',[.2 .1 .6 .8])
    % PlotMesh(nodes, elementPlot, 0);
    % v1 = reshape(VMn(:,mod), 2, []).';
    % PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', max(nodes(:,2)));
    % title(['\Phi_' num2str(mod)])

    % reduced assembly
    PROM_Assembly = DpromAssembly(MeshNominal,U,V);
    PROM_Assembly.DATA.M = PROM_Assembly.mass_matrix();  % reduced mass matrix 
    PROM_Assembly.DATA.C = V.'*Dn*V;    % reduced damping matrix, using C as needed by the residual function of Newmark integration
    PROM_Assembly.DATA.K = V.'*Kn*V;    % reduced stiffness matrix 

    % INTERNAL FORCES TENSORS _____________________________________________
    tensors_PROM = reduced_tensors_DpROM(NominalAssembly, elements, ...
    V, U, FORMULATION, VOLUME, USEJULIA)  ;
    
    % HYDRODYNAMIC FORCES _________________________________________________
    
    % find spine and tail elements
    if fishDim == 2
        [spineNodes, spineElements, spineElementWeights, nodeIdxPosInElements] = find_spine_TRI3(elements,nodes);
        [tailNode, tailElement, ~] = find_tail(elements,nodes,spineElements,nodeIdxPosInElements);
    else
        [spineNodes, spineElements, spineElementWeights, nodeIdxPosInElements] = find_spine_TET4(elements,nodes);
        [tailNode, tailElement, ~] = find_tail(elements,nodes,spineElements,nodeIdxPosInElements);
    end

    % get normalisation factors
    normalisationFactors = compute_normalisation_factors(nodes, elements, spineElements, nodeIdxPosInElements);
    wTail = normalisationFactors(tailElement);

    % tail pressure force: get matrices
    [A,B] = compute_AB_tail_pressure(nodeIdxPosInElements(tailElement,:));
    nodesTailEl = elements(tailElement,:);
    iDOFs = [nodesTailEl(1)*2-1,nodesTailEl(1)*2,...
             nodesTailEl(2)*2-1,nodesTailEl(2)*2,...
             nodesTailEl(3)*2-1,nodesTailEl(3)*2];
    VTail = V(iDOFs,:);
    UTail = U(iDOFs,:);

    R = [0 -1 0 0 0 0;
     1 0 0 0 0 0;
     0 0 0 -1 0 0;
     0 0 1 0 0 0;
     0 0 0 0 0 1;
     0 0 0 0 -1 0];     % 90 degrees rotation counterclock-wise

    % tail pressure force: group tail quantities in a struct
    tailProperties.A = A;
    tailProperties.B = B;
    tailProperties.R = R;
    tailProperties.mTilde = mTilde;
    tailProperties.w = wTail;
    tailProperties.V = VTail;
    tailProperties.U = UTail;
    tailProperties.tailNode = tailNode;
    tailProperties.tailElement = tailElement;
    tailProperties.iDOFs = iDOFs;

    % spine momentum change tensor (reduced order)
    spineTensors = compute_spine_momentum_tensors_PROM(PROM_Assembly, spineElementWeights,nodeIdxPosInElements,normalisationFactors,mTilde);
    spineProperties.tensors = spineTensors;
    spineProperties.spineNodes = spineNodes;
    spineProperties.spineElements = spineElements;

    % drag force (reduced order)
    [~,~,skinElements, skinElementFaces] = getSkin2D(elements);
    headNode = find_node_2D(0,0,nodes);
    headxDOF = 2*headNode-1;
    VHead = V(headxDOF,:);
    rho = 1000;
    kFactor = 0.05;
    tensors_drag = compute_drag_tensors_PROM(PROM_Assembly, skinElements, skinElementFaces, kFactor*rho,VHead) ;
    dragProperties.tensors = tensors_drag;
    dragProperties.skinElements = skinElements;
    dragProperties.skinElementFaces = skinElementFaces;
    

    % ACTUATION FORCES ____________________________________________________
    
    Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of the nominal fish
    Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of the nominal fish

    nel = size(elements,1);
    actuationDirection = [1;0;0];               %[1;0]-->[1;0;0] (Voigt notation)
    
    % top muscle
    topMuscle = zeros(nel,1);
    for el=1:nel
        elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2))/3;
        elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1))/3;
        if elementCenterY>0.00 &&  elementCenterX < -Lx*0.25 && elementCenterX > -Lx*0.8
            topMuscle(el) = 1;
        end    
    end
    actuTop = reduced_tensors_actuation_PROM(NominalAssembly, V, U, topMuscle, actuationDirection);
    
    % bottom muscle
    bottomMuscle = zeros(nel,1);
    for el=1:nel
        elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2))/3;
        elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1))/3;
        if elementCenterY<0.00 &&  elementCenterX < -Lx*0.25 && elementCenterX > -Lx*0.8
            bottomMuscle(el) = 1;
        end    
    end
    actuBottom = reduced_tensors_actuation_PROM(NominalAssembly, V, U, bottomMuscle, actuationDirection);


end 