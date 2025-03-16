% spine_force_FOM
%
% Synthax:
% f = spine_force_FOM(q,qd,Assembly,elements,tailProperties)
%
% Description: Computes the spine force for the FOM, expressed as an
% unconstrained vectors (size: nDOFs).
%
% INPUTS: 
% (1) q:                unconstrained displacement vector
% (2) qd:               unconstrained velocity vector
% (3) Assembly:         FOM FEM assembly
% (4) elements:         elements from the FEM assembly
% (5) spineProperties:  struct containing the mathematical property of the
%                       spine
%
% OUTPUTS:
% (1) f:                spine force, of size nDOFs x 1, where nDOFs is the 
%                       total, unconstraint number of DOFs
%     
% Last modified: 24/02/2023, Mathieu Dubied, ETH Zurich
function f = spine_force_FOM(q,qd,qdd,Assembly,elements,spineProperties)
    
    % get information from assembly
    if isempty(Assembly.Mesh.EBC)
         nUncDOFs = Assembly.Mesh.nDOFs;  % no constraints
    else
        nUncDOFs = size(Assembly.Mesh.EBC.unconstrainedDOFs,2);% some constraints 
    end
    nDOFs = Assembly.Mesh.nDOFs;
    x0 = reshape(Assembly.Mesh.nodes.',[],1);
    
    % get information about the spine
    spineNodes = spineProperties.spineNodes;
    spineElements = spineProperties.spineElements;
    nodeIdxPosInElements = spineProperties.nodeIdxPosInElements;
    matchedDorsalNodesIdx = spineProperties.dorsalNodeIdx;
    matchedDorsalNodesZPos = spineProperties.zPos;
    normalisationFactors = spineProperties.normalisationFactors;
    
    % initialise force
    f = zeros(nDOFs,1);
    
    % create force by looping over spine elements
    for i = 1:size(spineElements,1)
        
        % compute A and B matrix based on the spine element configuration
        sElementIdx = spineElements(i);
        spineNodeIndexInElement = nodeIdxPosInElements(sElementIdx,:);
        A = A_TET4(spineNodeIndexInElement(1));
        B = B_TET4(spineNodeIndexInElement(1),spineNodeIndexInElement(2));
        R = [0 -1 0 0 0 0 0 0 0 0 0 0;
             1 0 0 0 0 0 0 0 0 0 0 0;
             0 0 0 0 0 0 0 0 0 0 0 0;
             0 0 0 0 -1 0 0 0 0 0 0 0;
             0 0 0 1 0 0 0 0 0 0 0 0;
             0 0 0 0 0 0 0 0 0 0 0 0;
             0 0 0 0 0 0 0 -1 0 0 0 0;
             0 0 0 0 0 0 1 0 0 0 0 0;
             0 0 0 0 0 0 0 0 0 0 0 0;
             0 0 0 0 0 0 0 0 0 0 -1 0;
             0 0 0 0 0 0 0 0 0 1 0 0;
             0 0 0 0 0 0 0 0 0 0 0 0];     % 90 degrees rotation counterclock-wise around z axis
         
         
        % create matrix to select spine velocity and displacement from q
        nodesSpineElement = elements(sElementIdx,:);
        iDOFs = [nodesSpineElement(1)*3-2,nodesSpineElement(1)*3-1,nodesSpineElement(1)*3,...
             nodesSpineElement(2)*3-2,nodesSpineElement(2)*3-1,nodesSpineElement(2)*3,...
             nodesSpineElement(3)*3-2,nodesSpineElement(3)*3-1,nodesSpineElement(3)*3,...
             nodesSpineElement(4)*3-2,nodesSpineElement(4)*3-1,nodesSpineElement(4)*3];
         
        nodeSelMatrix = zeros(12,nDOFs);    % 12 = 4*3, specific for linear tet mesh
        for ii = 1:length(iDOFs)
            currentSpineDOF = iDOFs(ii);
            nodeSelMatrix(ii,currentSpineDOF) = 1;
        end
        
        % compute force at the spine element
        mTilde = 0.25*pi*1000*(matchedDorsalNodesZPos(sElementIdx)*2)^2;
        w = normalisationFactors(sElementIdx);
        
        qEl = nodeSelMatrix*q;
        qdEl = nodeSelMatrix*qd;
        qddEl = nodeSelMatrix*qdd;
        term1 = (A*qddEl)'*R*B*(x0(currentSpineDOF)+qEl)*R*R*B*(x0(currentSpineDOF)+qEl);
        term2 = (A*qdEl)'*R*B*qdEl*R*R*B*(x0(currentSpineDOF)+qEl);
        term3 = (A*qdEl)'*R*B*(x0(currentSpineDOF)+qEl)*R*R*B*qdEl;
        fSpineAtElement = -mTilde*w*(term1 + term2 + term3);
        
        % position element force in assembly vector
        f(iDOFs) = fSpineAtElement;
      
    end
                                    
end






