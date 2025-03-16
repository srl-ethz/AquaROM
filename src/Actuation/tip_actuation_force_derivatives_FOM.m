% tail_force_derivatives_FOM
%
% Synthax:
% der = tail_force_derivatives_FOM(q,qd,qdd,tensors)
%
% Description: This function returns the partial derivatives of the tip
% actuation force
% The partial derivatives are needed to solve the EoMs for the FOM.
%
% INPUTS: 
% (1) q:                vector of time domain displacements
% (2) t:                time
% (4) tensors:          tensors used to express the spine momentum change
%
% OUTPUTS:
% (1) der:              strucure array containing the partial derivatives 
%     
% Additional notes:
%   - q=u, qd=\dot{u}
%   
% Last modified: 18/1/2024, Mathieu Dubied, ETH Zürich
function der = tip_actuation_force_derivatives_FOM(t,q,Assembly,tailProperties)
    
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
    
    % dfdq ________________________________________________________________
    dfdq = sin(t*2*pi)*R*B; 

    
    % express at Assembly level
    d = length(tailDOFs);
    I{1} = kron(true(d,1), tailDOFs);
    J{1} = kron(tailDOFs, true(d,1));
    
    I = vertcat(I{:});
    J = vertcat(J{:});
    dfdq = vertcat(dfdq);
 
    % store results in output struct ______________________________________
    der.dfdq = sparse(I,J, dfdq, nDOFs, nDOFs);
    
end
