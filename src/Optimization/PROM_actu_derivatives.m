% PROM_actu_derivatives
%
% Synthax:
% der = PROM_actu_derivatives(tensors_PROM,a)
%
% Description: This function returns the partial derivatives of the
% sensitity ODE. As f_act is linear in eta and xi, the first order partial
% derivatives are constant that are not dependent on eta or xi. However,
% they depedent on the actuation value a.
%
% INPUTS: 
% (1) tensors_PROM: structure array containing reduced tensors of
%                   F_act model
% (2) a:            actuation signal
%
% OUTPUTS:
% (1) der:              strucure array containing first order partial 
%                       derivatives needed to solve the sensitivity ODE. 
%
% Last modified: 15/10/2023, Mathieu Dubied, ETH Zurich
function der = PROM_actu_derivatives(tensors_PROM,a)

% evaluate partial derivatives
dfdp = a*tensors_PROM.B3;
dfdq = a*tensors_PROM.B2;

% store results in output struct
der.dfdp = double(dfdp);
der.dfdq = double(dfdq);

end
