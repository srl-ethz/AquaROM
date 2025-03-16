% FOM_hydro_derivatives
%
% Synthax:
% der = FOM_hydro_derivatives(q,qd,tensors_ROM)
%
% Description: This function returns the partial derivatives of the
% sensitity ODE, evaluated at the solution displacements q and velocities qdot.
%
% INPUTS: 
% (1) q:                vector of time domain displacements
% (2) qdot:             vector of time domain velocities
% (3) tensors_ROM:      structure array containing reduced tensors of
%                       F_hydro model
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
% Last modified: 16/04/2023, Mathieu Dubied, ETH Zürich
function der = FOM_hydro_derivatives(q,qd,Tu2,Tudot2,Tuu3,Tuudot3,Tudotudot3)


% First order derivatives _________________________________________________
% evaluate partial derivatives

df1dq = zeros(size(Tu2));
df1dqd = zeros(size(Tu2));
df2dq = Tu2 ;
df2dqd = Tudot2 ;
df3dq = double(ttv(Tuu3,q,3) + ttv(Tuu3,q,2) + ttv(Tuudot3,qd,3));
df3dqd = double(ttv(Tudotudot3,qd,3) + ttv(Tudotudot3,qd,2) + ttv(Tuudot3,q,2));

dfdq = df1dq + df2dq + df3dq;
dfdqd = df1dqd + df2dqd + df3dqd;a

% store results in output struct
der.dfdq = double(dfdq);
der.dfdqd = double(dfdqd);




end
