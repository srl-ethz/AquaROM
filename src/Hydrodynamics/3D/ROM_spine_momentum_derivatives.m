% ROM_spine_momentum_derivatives
%
% Synthax:
% der = ROM_spine_momentum_derivatives(q,qd,qdd,x0,A,B)
%
% Description: This function returns the partial derivatives of the change
% in the fish momentum, which is understood as a force acting on the spine.
% The partial derivatives are needed to solve the sensitivity ODE as well 
% as the EoMs.
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
%   - q,qd, qdd should be understood as eta and dot{eta}, ddot{eta}.
%
% Last modified: 19/10/2023, Mathieu Dubied, ETH Zürich
function der = ROM_spine_momentum_derivatives(q,qd,qdd,tensors)
    % get tensors   
    Txx = tensors.Txx;
    TxV = tensors.TxV;
    TVx = tensors.TVx;
    TVV = tensors.TVV;

    if isa(Txx,'struct')    % the tensors were computed for a PROM
        Txx = Txx.f0;
        TxV = TxV.f0;
        TVx = TVx.f0;
        TVV = TVV.f0;
    end

    
    % dfdq ________________________________________________________________
    dfdq = ttv(TxV,qdd,2) + ttv(TVx,qdd,2) ...
         + ttv(ttv(TVV,q,4),qdd,2) ...
         + ttv(ttv(TVV,q,3),qdd,2) ...
         + ttv(ttv(TVV,qd,3),qd,2) ...
         + ttv(ttv(TVV,qd,4),qd,2);
    
    % dfdqd _______________________________________________________________
    dfdqd = ttv(TVx,qd,3) + ttv(TVx,qd,2) ...
            + ttv(ttv(TVV,q,4),qd,3) + ttv(ttv(TVV,q,4),qd,2) ...
            + ttv(TxV,qd,3) + ttv(TxV,qd,2) ...
            + ttv(ttv(TVV,qd,4),q,3) + ttv(ttv(TVV,q,3),qd,2);
            

    % dfdqdd ______________________________________________________________
    dfdqdd = Txx + ttv(TxV,q,3) + ttv(TVx,q,3) + ttv(ttv(TVV,q,4),q,3);
 
    
    % store results in output struct ______________________________________
    der.dfdq = double(dfdq);
    der.dfdqd = double(dfdqd);
    der.dfdqdd = double(dfdqdd);

end
