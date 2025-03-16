% ROM_drag_derivatives
%
% Synthax:
% der = ROM_drag_derivatives(qd,tensors)
%
% Description: This function returns the partial derivatives of the drag
% force.
% The partial derivatives are needed to solve the sensitivity ODE as well 
% as the EoMs.
%
% INPUTS: 
% (1) qd:               vector of time domain velocities
% (2) tensors:          tensors used to express the drag force
%
% OUTPUTS:
% (1) der:              strucure array containing the partial derivatives 
%     
% Additional notes:
%   - q,qd, qdd should be understood as eta and dot{eta}, ddot{eta}.
%
% Last modified: 27/10/2023, Mathieu Dubied, ETH Zürich
function der = ROM_drag_derivatives(qd,tensors)
    % get tensors   
    T3 = tensors.Tr3;
    
    % dfdq ________________________________________________________________
    % zero
    
    % dfdqd _______________________________________________________________
    dfdqd = 2*ttv(T3,qd,3);
            

    % dfdqdd ______________________________________________________________
    % zero
 
    
    % store results in output struct ______________________________________
    der.dfdqd = double(dfdqd);
    

end
