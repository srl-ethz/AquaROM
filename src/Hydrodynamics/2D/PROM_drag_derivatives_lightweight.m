% PROM_drag_derivatives_lightweigh
%
% Synthax:
% der = PROM_drag_derivatives_lightweigh(qd,tensors)
%
% Description: This function returns the partial derivatives of the drag
% force.
% This lightweight version only entails the partial derivatives needed in
% the algorithm proposed in the paper: dfdp for xi=0
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
function der = PROM_drag_derivatives_lightweight(qd,tensors)

    % get tensors   
    T4 = tensors.Tr4;
    
    % dfdp ________________________________________________________________
    dfdp = ttv(ttv(T4,qd,4),qd,3);
        
    % store results in output struct ______________________________________
    der.dfdp = double(dfdp);

end

