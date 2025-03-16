% DpROM_hydro_derivatives
%
% Synthax:
% der = DpROM_hydro_derivatives(q,qd,xi,tensors_DpROM)
%
% Description: This function returns the partial derivatives of the
% sensitity ODE, evaluated at the solution displacements q and velocities qdot.
%
% INPUTS: 
% (1) q:                vector of time domain displacements
% (2) qdot:             vector of time domain velocities
% (3) tensors_DpROM:    structure array containing reduced tensors of
%                       F_hydro model
% (4) FOURTHORDER       logical value, 1 if fourth order tensors need to be
%                       computed, 0 else
% (5) secondOrderDer    logical value, 1 if second order derivatives need
%                       to be computed, 0 else
%
% OUTPUTS:
% (1) der:              strucure array containing first (and possibly
%                       second) order partial derivatives needed to solve 
%                       the sensitivity ODE. 
%     
%
% Additional notes:
%   - q and qd should be understood as eta and dot{eta}.
%   - p should be understood as xi.
%
% Last modified: 21/03/2023, Mathieu Dubied, ETH Zürich
function der = DpROM_hydro_derivatives(q,qd,xi,tensors_PROM,FOURTHORDER,secondOrderDer)
% populate tensors
Tr2 = tensors_PROM.Tr2; 
Tru2 = tensors_PROM.Tru2;
Tru3 = tensors_PROM.Tru3;
Trudot2 = tensors_PROM.Trudot2;
Trudot3 = tensors_PROM.Trudot3;
Truu3 = tensors_PROM.Truu3;
Truudot3 = tensors_PROM.Truudot3;
Trudotudot3 = tensors_PROM.Trudotudot3;
if FOURTHORDER
    Truu4 = tensors_PROM.Truu4;
    Truudot4 = tensors_PROM.Truudot4;
    Trudotudot4 = tensors_PROM.Trudotudot4;
end

% First order derivatives _________________________________________________
% evaluate partial derivatives
df1dp = Tr2;
df2dp = ttv(Tru3,q,2) + ttv(Trudot3,qd,2);

dfdp = df1dp + double(df2dp);
if FOURTHORDER
    df3dp = ttv(ttv(Truu4,q,3),q,2) + ...
            ttv(ttv(Truudot4,qd,3),q,2) + ...
            ttv(ttv(Trudotudot4,qd,3),qd,2);
    dfdp = dfdp + double(df3dp);
end

df1dq = zeros(size(Tr2,1));
df1dqd = zeros(size(Tr2,1));
df2dq = Tru2 + double(ttv(Tru3,xi,3));
df2dqd = Trudot2 + double(ttv(Trudot3,xi,3));
df3dq = double(ttv(Truu3,q,3) + ttv(Truu3,q,2) + ttv(Truudot3,qd,3));
df3dqd = double(ttv(Trudotudot3,qd,3) + ttv(Trudotudot3,qd,2) + ttv(Truudot3,q,2));
if FOURTHORDER
    df3dq = df3dq + double(ttv(ttv(Truu4,xi,4),q,3) + ...
                            ttv(ttv(Truu4,xi,4),q,2) + ...
                            ttv(ttv(Truudot4,xi,4),qd,3));  
    df3dqd = df3dqd +  double(ttv(ttv(Trudotudot4,xi,4),q,3) + ...
                                ttv(ttv(Trudotudot4,xi,4),q,2) + ...
                                ttv(ttv(Truudot4,xi,4),q,2));
end


dfdq = df1dq + df2dq + df3dq;
dfdqd = df1dqd + df2dqd + df3dqd;

% store results in output struct
der.dfdp = double(dfdp);
der.dfdq = double(dfdq);
der.dfdqd = double(dfdqd);

% Second order derivatives ________________________________________________
% only compute the non-zero derivatives. Syntax: df3dqdqd= dell(f3)/(dell dqdot dell dq)
% The script could be simplified, but it should be easier to follow written
% as it is
if secondOrderDer
    % without 4th order tensors
    df3dqdq = Truu3 + permute(Truu3,[1 3 2]);
    df3dqdqd = Truudot3;
    df3dqddq = permute(Truudot3,[1 3 2]);
    df3dqddqd = Trudotudot3 + permute(Trudotudot3,[1 3 2]);

    df2dqdp = Tru3;
    df2dpdq = permute(Tru3,[1,3,2]);
    df2dqddp = Trudot3;
    df2dpdqd = permute(Trudot3,[1,3,2]);
    
    % with 4th order tensors
    if FOURTHORDER
        df3dqdq = df3dqdq + ttv(Truu4,xi,4) + permute(ttv(Truu4,xi,4),[1 3 2]);
        df3dqdqd = df3dqdqd + ttv(Truudot4,xi,4);
        df3dqddq = df3dqddq + permute(ttv(Truudot4,xi,4),[1 3 2]);
        df3dqddqd = df3dqddqd + ttv(Trudotudot4,xi,4) + permute(ttv(Trudotudot4,xi,4),[1 3 2]);

        df3dqdp = ttv(Truu4,q,3) + ttv(Truu4,q,2) + ttv(Truudot4,qd,3);
        df3dpdq = permute(ttv(Truu4,q,3),[1 3 2]) + permute(ttv(Truu4,q,2),[1 3 2]) + permute(ttv(Truudot4,qd,3), [1 3 2]);
        df3dqddp = ttv(Trudotudot4,qd,3) + ttv(Trudotudot4,qd,2) + ttv(Truudot4,q,2);
        df3dpdqd = permute(ttv(Truudot4,q,2),[1 3 2]) + permute(ttv(Trudotudot4,qd,3),[1 3 2]) + permute(ttv(Trudotudot4,qd,2),[1 3 2]);
    end

    % assemble contributions of f1r,f2r,f3r
    dfdqdq = df3dqdq;
    dfdqdqd = df3dqdqd; 
    dfdqddq = df3dqddq;
    dfdqddqd = df3dqddqd;

    dfdqdp = df2dqdp;
    dfdpdq = df2dpdq; 
    dfdqddp = df2dqddp; 
    dfdpdqd = df2dpdqd; 

    % with 4th order tensors
    if FOURTHORDER
        dfdqdp = dfdqdp + df3dqdp;
        dfdpdq = dfdpdq + df3dpdq; 
        dfdqddp = dfdqddp + df3dqddp; 
        dfdpdqd = dfdpdqd + df3dpdqd; 
    end
    
    % store results in output struct
    der.dfdqdq = dfdqdq;
    der.dfdqdqd = dfdqdqd; 
    der.dfdqddq = dfdqddq;
    der.dfdqddqd = dfdqddqd;

    der.dfdqdp = dfdqdp;
    der.dfdpdq = dfdpdq; 
    der.dfdqddp = dfdqddp; 
    der.dfdpdqd = dfdpdqd; 

end


end
