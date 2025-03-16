% drag_force_FOM
%
% Synthax:
% f = drag_force_FOM(qd,Assembly,elements,dragProperties)
%
% Description: Computes the drag force for the FOM, expressed as an
% unconstraint vectors (size: nDOFs).
%
% INPUTS: 
% (1) qd:               unconstrained velocity vector
% (2) dragProperties:   struct containing the mathematical property of the
%                       drag
%
% OUTPUTS:
% (1) f:                drag force, of size nDOFs x 1, where nDOFs is the 
%                       total, unconstraint number of DOFs
%     
% Last modified: 18/02/2023, Mathieu Dubied, ETH Zurich
function f = drag_force_FOM(qd,dragProperties)

    % get information from drag properties
    T3 = dragProperties.tensors.T3;
    headxDOF = dragProperties.headxDOF;
    
    % compute force
    f = T3*qd(headxDOF)^2;
    
end






