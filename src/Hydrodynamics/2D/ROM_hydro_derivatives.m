% ROM_hydro_derivatives
%
% Synthax:
% der =ROM_hydro_derivatives(q,qd,xi,tensors_ROM)
%
% Description: This function returns the partial derivatives of the
% sensitity ODE, evaluated at the solution displacements q and velocities qdot.
%
% INPUTS: 
% (1) q:                vector of time domain displacements
% (2) qdot:             vector of time domain velocities
% (3) tensors_ROM:      structure array containing reduced tensors of
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
%   - NOT TESTED FOR secondOrderDer = 1 
%
% Last modified: 16/04/2023, Mathieu Dubied, ETH Zürich
function der = ROM_hydro_derivatives(q,qd,xi,tensors_ROM,FOURTHORDER,secondOrderDer)
% populate tensors
 
Tru2 = tensors_ROM.Tru2;
Trudot2 = tensors_ROM.Trudot2;
Truu3 = tensors_ROM.Truu3;
Truudot3 = tensors_ROM.Truudot3;
Trudotudot3 = tensors_ROM.Trudotudot3;

% First order derivatives _________________________________________________
% evaluate partial derivatives

df1dq = zeros(size(Tru2));
df1dqd = zeros(size(Tru2));
df2dq = Tru2 ;
df2dqd = Trudot2 ;
df3dq = double(ttv(Truu3,q,3) + ttv(Truu3,q,2) + ttv(Truudot3,qd,3));
df3dqd = double(ttv(Trudotudot3,qd,3) + ttv(Trudotudot3,qd,2) + ttv(Truudot3,q,2));

dfdq = df1dq + df2dq + df3dq;
dfdqd = df1dqd + df2dqd + df3dqd;

% store results in output struct
der.dfdq = double(dfdq);
der.dfdqd = double(dfdqd);




end
