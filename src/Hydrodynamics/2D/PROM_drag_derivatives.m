% PROM_drag_derivatives
%
% Synthax:
% der = PROM_drag_derivatives(qd,tensors)
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
% Last modified: 07/11/2023, Mathieu Dubied, ETH Zürich
function der = PROM_drag_derivatives(qd,xi,tensors)

    % get tensors   
    T3 = tensors.Tr3;
    T4 = tensors.Tr4;
    T5 = tensors.Tr5;
    
    % dfdp ________________________________________________________________
    dfdp = ttv(ttv(T4,qd,4),qd,3) ...
         + ttv(ttv(ttv(T5,xi,3),qd,4),qd,3) ...
         + ttv(ttv(ttv(T5,xi,2),qd,4),qd,3);
    
    % dfdq ________________________________________________________________
    % zero
    
    % dfdqd _______________________________________________________________
    dfdqd = 2*ttv(T3,qd,3) ...
            + 2*ttv(ttv(T4,xi,2),qd,3) ...
            + 2*ttv(ttv(ttv(T5,xi,2),xi,2),qd,3);      

    % dfdqdd ______________________________________________________________
    % zero
 

    % store results in output struct ______________________________________
    der.dfdp = double(dfdp);
    der.dfdqd = double(dfdqd);

    

end

