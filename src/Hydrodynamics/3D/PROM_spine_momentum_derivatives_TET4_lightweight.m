% PROM_spine_momentum_derivatives_TET4_lightweight
%
% Synthax:
% der = PROM_spine_momentum_derivatives_TET4_lightweight(q,qd,qdd,eta0,xi,tensors)
%
% Description: This function returns the partial derivatives of the change
% in the fish momentum, which is understood as a force acting on the spine.
% This lightweight version only entails the partial derivatives needed in
% the algorithm proposed in the paper: dfdp for xi=0
%
% INPUTS: 
% (1) q:                vector of time domain displacements
% (2) qd:               vector of time domain velocities
% (3) qdd:              vector of time domain acceleration
% (4) xi:               shape variation parameter vector/scalar
% (5) tensors:          tensors (struct) used to express the spine momentum
%                       change
%
% OUTPUTS:
% (1) der:              strucure array containing the partial derivatives 
%     
% Additional notes:
%   - q,qd, qdd should be understood as eta and dot{eta}, ddot{eta}.
%
% Last modified: 06/11/2023, Mathieu Dubied, ETH Zürich
function der = PROM_spine_momentum_derivatives_TET4_lightweight(q,qd,qdd,tensors)

    % get tensors   
    % Txx = tensors.Txx.f0;
    % TxV = tensors.TxV.f0;
    % TVx = tensors.TVx.f0;
    % TVV = tensors.TVV.f0;
    TUV = tensors.TUV.f0;
    TVU = tensors.TVU.f0;
    % TUU = tensors.TUU.f0;
    TxU = tensors.TxU.f0;
    TUx = tensors.TUx.f0;

    Txx3 = tensors.Txx.f1;
    TxV4 = tensors.TxV.f1;
    TVx4 = tensors.TVx.f1;
    TVV5 = tensors.TVV.f1;
    % TUV5 = tensors.TUV.f1;
    % TVU5 = tensors.TVU.f1;
    % % TUU5 = tensors.TUU.f1;
    % TxU4 = tensors.TxU.f1;
    % TUx4 = tensors.TUx.f1;


    % dfdp ________________________________________________________________
    df0dp = ttv(TxU,qdd,2) ...
            + ttv(TUx,qdd,2) ...
            + ttv(ttv(TUV,q,4),qdd,2) ...
            + ttv(ttv(TVU,q,3),qd,2) ...
            + ttv(ttv(TVU,qd,3),qd,2) ...
            + ttv(ttv(TUV,qd,4),qd,2);

    df1dp = ttv(Txx3,qdd,2) ...
            + ttv(ttv(TxV4,q,3),qdd,2) ...
            + ttv(ttv(TVx4,q,3),qdd,2) ...
            + ttv(ttv(ttv(TVV5,q,4),q,3),qdd,2) ...
            + ttv(ttv(TVx4,qd,3),qd,2) ...
            + ttv(ttv(ttv(TVV5,q,4),qd,3),qd,2) ...
            + ttv(ttv(TxV4,qd,3),qd,2) ...
            + ttv(ttv(ttv(TVV5,qd,4),q,3),qd,2);

    dfdp = df0dp + df1dp;
    
    % store results in output struct ______________________________________
    der.dfdp = double(dfdp);

end