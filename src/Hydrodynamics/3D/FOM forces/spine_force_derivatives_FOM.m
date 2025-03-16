% spine_force_derivatives_FOM
%
% Synthax:
% der = spine_force_derivatives_FOM(q,qd,qdd,tensors)
%
% Description: This function returns the partial derivatives of the change
% in the fish momentum, which is understood as a force acting on the spine.
% The partial derivatives are needed to solve the EoMs for the FOM.
%
% INPUTS: 
% (1) q:                vector of time domain displacements
% (2) qd:               vector of time domain velocities
% (3) qdd:              vector of time domain acceleration
% (4) tensors:          tensors used to express the spine momentum change
%
% OUTPUTS:
% (1) der:              strucure array containing the partial derivatives 
%     
% Additional notes:
%   - q=u, qd=\dot{u}, qdd=\ddot{u}
%   
% Last modified: 27/10/2024, Mathieu Dubied, ETH Zürich
function der = spine_force_derivatives_FOM(q,qd,qdd,tensors)

    % get tensors   
    T2a = tensor(tensors.T2a);
    T3b = tensors.T3b;
    T3c = tensors.T3c;
    T4d = tensors.T4d;
    
    % dfdq ________________________________________________________________
    dfdq = ttv(T3b,qdd,2) + ttv(T3c,qdd,2) ...
         + ttv(ttv(T4d,q,4),qdd,2) ...
         + ttv(ttv(T4d,q,3),qdd,2) ...
         + ttv(ttv(T4d,qd,3),qd,2) ...
         + ttv(ttv(T4d,qd,4),qd,2);
    
    % dfdqd _______________________________________________________________
    dfdqd = ttv(T3c,qd,3) + ttv(T3c,qd,2) ...
            + ttv(ttv(T4d,q,4),qd,3) + ttv(ttv(T4d,q,4),qd,2) ...
            + ttv(T3b,qd,3) + ttv(T3b,qd,2) ...
            + ttv(ttv(T4d,qd,4),q,3) + ttv(ttv(T4d,q,3),qd,2);
            
    % dfdqdd ______________________________________________________________
    dfdqdd = T2a  + ttv(T3b,q,3) + ttv(T3c,q,3) + ttv(ttv(T4d,q,4),q,3);
 
    % store results in output struct ______________________________________
    der.dfdq = double(dfdq);
    der.dfdqd = double(dfdqd);
    der.dfdqdd = double(dfdqdd);

end
