% tail_force_FOM
%
% Synthax:
% f = tail_force_FOM(q,qd,Assembly,elements,tailProperties)
%
% Description: Computes the tail force for the FOM, expressed as an
% unconstraint vectors (size: nDOFs).
%
% INPUTS: 
% (1) q:                constraint displacement vector
% (2) qd:               constraints velocity vector
% (3) Assembly:         FOM FEM assembly
% (4) elements:         elements from the FEM assembly
% (5) tailProperties:   struct containing the mathematical property of the
%                       tail and tail force
%
% OUTPUTS:
% (1) f:                tail force, of size nDOFs x 1, where nDOFs is the 
%                       total, unconstraint number of DOFs
%     
% Last modified: 24/02/2023, Mathieu Dubied, ETH Zurich
function f = tail_force_FOM(q,qd,Assembly,elements,tailProperties)
    
    % get information from assembly
    if isempty(Assembly.Mesh.EBC)
         nUncDOFs = Assembly.Mesh.nDOFs;  % no constraints
    else
        nUncDOFs = size(Assembly.Mesh.EBC.unconstrainedDOFs,2);% some constraints 
    end
    
    nDOFs = Assembly.Mesh.nDOFs;
    x0 = reshape(Assembly.Mesh.nodes.',[],1);
    
    % get information from tail
    A = tailProperties.A;
    B = tailProperties.B;
    R = tailProperties.R;
    wTail = tailProperties.w;
    tailNode = tailProperties.tailNode;
    tailElement = tailProperties.tailElement;
    tailDOFs = tailProperties.iDOFs;
    mTilde = 0.25*pi*1000*(tailProperties.z*2)^2;
    tailProperties.mTilde = mTilde;
    
    % create matrix to select tail velocity and displacement from q
    nodeSelMatrix = zeros(12,nDOFs);    % 12 = 4*3, specific for linear tet mesh
    for idx = 1:length(tailDOFs)
        currentTailDOF = tailDOFs(idx);
        nodeSelMatrix(idx,currentTailDOF) = 1;
    end
    
    % compute tail force expressed at the tail element
    tailNodeInTailElement = find(elements(tailElement,:)==tailNode);
    
    fTailAtElement = 0.5*mTilde*wTail^3*(dot(A*nodeSelMatrix*qd,R*B*...
        (x0(tailDOFs)+nodeSelMatrix*q))).^2*B*(x0(tailDOFs)+nodeSelMatrix*q);
    
    % return unconstraint force
    f = zeros(nDOFs,1);
    f(tailNode*3-2:tailNode*3) = fTailAtElement(tailNodeInTailElement*3-2:tailNodeInTailElement*3);
                      
end






