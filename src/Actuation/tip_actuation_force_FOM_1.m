% tip_actuation_force_FOM_1
%
% Synthax:
% f = tip_actuation_force_FOM_1(q,qd,Assembly,elements,tailProperties)
%
% Description: Computes a simple tip actuation force for the FOM, expressed 
% as an unconstraint vectors (size: nDOFs).
%
% INPUTS: 
% (1) t:                time
% (1) k:                scaling factor
% (3) Assembly:         FOM FEM assembly
% (4) tailProperties:   struct containing the mathematical property of the
%                       tail and tail force
%
% OUTPUTS:
% (1) f:                tip force, of size nDOFs x 1, where nDOFs is the 
%                       total, unconstraint number of DOFs
%     
% Last modified: 22/11/2024, Mathieu Dubied, ETH Zurich
function f = tip_actuation_force_FOM_1(t,k,Assembly,tailProperties)
    
    nDOFs = Assembly.Mesh.nDOFs;
    % get information from tail
    tailNode = tailProperties.tailNode;

    % return unconstrained force
    f = zeros(nDOFs,1);
    f(tailNode*3-1) = k*sin(t*2*pi);  % directed in y direction
                      
end






